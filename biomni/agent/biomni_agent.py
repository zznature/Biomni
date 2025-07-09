from biomni.llm import get_llm
from biomni.utils import pretty_print, api_schema_to_langchain_tool, \
    function_to_api_schema, textify_api_dict, library_content_dict, \
    data_lake_dict, load_pkl, run_bash_script, run_r_code, run_with_timeout
from biomni.tool.tool_registry import ToolRegistry
from biomni.model.retriever import ToolRetriever
from biomni.tool.support_tools import run_python_repl

from typing import List, Dict, Any, Optional, Annotated, Sequence, TypedDict, Literal
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage, SystemMessage
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langgraph.graph import StateGraph, Graph, START, END
from langgraph.checkpoint.memory import MemorySaver

import os, re, glob, signal, inspect, numpy as np, subprocess, tempfile
from IPython.display import Image, display
from multiprocessing import Process, Queue

class AgentState(TypedDict):
    messages: List[BaseMessage]
    next_step: Optional[str]

class biomni_agent:
    def __init__(self, path = './data', llm='claude-3-7-sonnet-latest', 
                 use_tool_retriever=True, timeout_seconds=600):
        """
        Initialize the biomni agent.
        
        Args:
            path: Path to the data
            llm: LLM to use for the agent
            use_tool_retriever: If True, use a tool retriever
        """
        self.path = path
        
        if not os.path.exists(path): 
            os.makedirs(path)
            print(f'Created directory: {path}')
            ### TODO: Download the data
        else:
            print(f'Data directory already exists: {path}, loading...')

        module2api = load_pkl(f'{path}/api/all_250415.pkl')

        self.llm = get_llm(llm, stop_sequences = ['</execute>', '</solution>'])
        api = [i for i in module2api['biomni.tool.support_tools'] if i['name'] == 'run_python_repl'][0]
        self.tools = [api_schema_to_langchain_tool(api, mode = 'custom_tool', module_name = 'biomni.tool.support_tools')]
        self.module2api = module2api
        self.use_tool_retriever = use_tool_retriever
        
        if self.use_tool_retriever:
            self.tool_registry = ToolRegistry(module2api)
            self.retriever = ToolRetriever()

        # Add timeout parameter
        self.timeout_seconds = timeout_seconds  # 10 minutes default timeout
        self.configure()

    def add_tool(self, api):
        function_code = inspect.getsource(api)
        schema = function_to_api_schema(function_code, self.llm)
        ## TODO: add a new tool to the tool registry
        #self.tools += [api_schema_to_langchain_tool(schema, mode = 'custom_tool', module_name = api.__module__)]

    def add_data(self, data):
        ## TODO: add a new data to the data lake
        #self.data_lake += [data]
        pass

    def add_software(self, software):
        ## TODO: add a new software to the software library
        #self.software_library += [software]
        pass

    def _generate_system_prompt(self, tool_desc, data_lake_content, library_content_list, self_critic=False, is_retrieval=False):
        """
        Generate the system prompt based on the provided resources.
        
        Args:
            tool_desc: Dictionary of tool descriptions
            data_lake_content: List of data lake items
            library_content_list: List of libraries
            self_critic: Whether to include self-critic instructions
            is_retrieval: Whether this is for retrieval (True) or initial configuration (False)
            
        Returns:
            The generated system prompt
        """
        def format_item_with_description(name, description):
            """Format an item with its description in a readable way."""
            # Handle None or empty descriptions
            if not description:
                description = f"Data lake item: {name}"
                
            # Check if the item is already formatted (contains a colon)
            if isinstance(name, str) and ': ' in name:
                return name
                
            # Wrap long descriptions to make them more readable
            max_line_length = 80
            if len(description) > max_line_length:
                # Simple wrapping for long descriptions
                wrapped_desc = []
                words = description.split()
                current_line = ""
                
                for word in words:
                    if len(current_line) + len(word) + 1 <= max_line_length:
                        if current_line:
                            current_line += " " + word
                        else:
                            current_line = word
                    else:
                        wrapped_desc.append(current_line)
                        current_line = word
                
                if current_line:
                    wrapped_desc.append(current_line)
                
                # Join with newlines and proper indentation
                formatted_desc = f"{name}:\n  " + "\n  ".join(wrapped_desc)
                return formatted_desc
            else:
                return f"{name}: {description}"
        
        # Format the data lake content
        if isinstance(data_lake_content, list) and all(isinstance(item, str) for item in data_lake_content):
            # Simple list of strings - check if they already have descriptions
            data_lake_formatted = []
            for item in data_lake_content:
                # Check if the item already has a description (contains a colon)
                if ': ' in item:
                    data_lake_formatted.append(item)
                else:
                    description = self.data_lake_dict.get(item, f"Data lake item: {item}")
                    data_lake_formatted.append(format_item_with_description(item, description))
        else:
            # List with descriptions
            data_lake_formatted = []
            for item in data_lake_content:
                if isinstance(item, dict):
                    name = item.get('name', '')
                    description = self.data_lake_dict.get(name, f"Data lake item: {name}")
                    data_lake_formatted.append(format_item_with_description(name, description))
                else:
                    # Check if the item already has a description (contains a colon)
                    if isinstance(item, str) and ': ' in item:
                        data_lake_formatted.append(item)
                    else:
                        description = self.data_lake_dict.get(item, f"Data lake item: {item}")
                        data_lake_formatted.append(format_item_with_description(item, description))
        
        # Format the library content
        if isinstance(library_content_list, list) and all(isinstance(item, str) for item in library_content_list):
            if len(library_content_list) > 0 and isinstance(library_content_list[0], str) and ',' not in library_content_list[0]:
                # Simple list of strings
                libraries_formatted = []
                for lib in library_content_list:
                    description = self.library_content_dict.get(lib, f"Software library: {lib}")
                    libraries_formatted.append(format_item_with_description(lib, description))
            else:
                # Already formatted string
                libraries_formatted = library_content_list
        else:
            # List with descriptions
            libraries_formatted = []
            for lib in library_content_list:
                if isinstance(lib, dict):
                    name = lib.get('name', '')
                    description = self.library_content_dict.get(name, f"Software library: {name}")
                    libraries_formatted.append(format_item_with_description(name, description))
                else:
                    description = self.library_content_dict.get(lib, f"Software library: {lib}")
                    libraries_formatted.append(format_item_with_description(lib, description))
        
        # Base prompt
        prompt_modifier = """
You are a helpful biomedical assistant assigned with the task of problem-solving. 
To achieve this, you will be using an interactive coding environment equipped with a variety of tool functions, data, and softwares to assist you throughout the process.

Given a task, make a plan first. The plan should be a numbered list of steps that you will take to solve the task. Be specific and detailed.
Format your plan as a checklist with empty checkboxes like this:
1. [ ] First step
2. [ ] Second step
3. [ ] Third step

Follow the plan step by step. After completing each step, update the checklist by replacing the empty checkbox with a checkmark:
1. [✓] First step (completed)
2. [ ] Second step
3. [ ] Third step

If a step fails or needs modification, mark it with an X and explain why:
1. [✓] First step (completed)
2. [✗] Second step (failed because...)
3. [ ] Modified second step
4. [ ] Third step

Always show the updated plan after each step so the user can track progress.

At each turn, you should first provide your thinking and reasoning given the conversation history.
After that, you have two options:

1) Interact with a programming environment and receive the corresponding output within <observe></observe>. Your code should be enclosed using "<execute>" tag, for example: <execute> print("Hello World!") </execute>. IMPORTANT: You must end the code block with </execute> tag.
   - For Python code (default): <execute> print("Hello World!") </execute>
   - For R code: <execute> #!R\nlibrary(ggplot2)\nprint("Hello from R") </execute>
   - For Bash scripts and commands: <execute> #!BASH\necho "Hello from Bash"\nls -la </execute>
   - For CLI softwares, use Bash scripts.

2) When you think it is ready, directly provide a solution that adheres to the required format for the given task to the user. Your solution should be enclosed using "<solution>" tag, for example: The answer is <solution> A </solution>. IMPORTANT: You must end the solution block with </solution> tag.

You have many chances to interact with the environment to receive the observation. So you can decompose your code into multiple steps. 
Don't overcomplicate the code. Keep it simple and easy to understand.
When writing the code, please print out the steps and results in a clear and concise manner, like a research log.
When calling the existing python functions in the function dictionary, YOU MUST SAVE THE OUTPUT and PRINT OUT the result.
For example, result = understand_scRNA(XXX) print(result) 
Otherwise the system will not be able to know what has been done.

For R code, use the #!R marker at the beginning of your code block to indicate it's R code.
For Bash scripts and commands, use the #!BASH marker at the beginning of your code block. This allows for both simple commands and multi-line scripts with variables, conditionals, loops, and other Bash features.

In each response, you must include EITHER <execute> or <solution> tag. Not both at the same time. Do not respond with messages without any tags. No empty messages.
"""

        # Add self-critic instructions if needed
        if self_critic:
            prompt_modifier += """
You may or may not receive feedbacks from human. If so, address the feedbacks by following the same procedure of multiple rounds of thinking, execution, and then coming up with a new solution.
"""

        # Add environment resources
        prompt_modifier += """

Environment Resources:

- Function Dictionary:
{function_intro}
--- 
{tool_desc}
---

{import_instruction}

- Biological data lake
You can access a biological data lake at the following path: {data_lake_path}. 
{data_lake_intro}
Each item is listed with its description to help you understand its contents.
----
{data_lake_content}
----

- Software Library:
{library_intro}
Each library is listed with its description to help you understand its functionality.
----
{library_content_formatted}
----

- Note on using R packages and Bash scripts:
  - R packages: Use subprocess.run(['Rscript', '-e', 'your R code here']) in Python, or use the #!R marker in your execute block.
  - Bash scripts and commands: Use the #!BASH marker in your execute block for both simple commands and complex shell scripts with variables, loops, conditionals, etc.
        """
        
        # Set appropriate text based on whether this is initial configuration or after retrieval
        if is_retrieval:
            function_intro = "Based on your query, I've identified the following most relevant functions that you can use in your code:"
            data_lake_intro = "Based on your query, I've identified the following most relevant datasets:"
            library_intro = "Based on your query, I've identified the following most relevant libraries that you can use:"
            import_instruction = "IMPORTANT: When using any function, you MUST first import it from its module. For example:\nfrom [module_name] import [function_name]"
        else:
            function_intro = "In your code, you will need to import the function location using the following dictionary of functions:"
            data_lake_intro = "You can write code to understand the data, process and utilize it for the task. Here is the list of datasets:"
            library_intro = "The environment supports a list of libraries that can be directly used. Do not forget the import statement:"
            import_instruction = ""
        
        # Format the content consistently for both initial and retrieval cases
        library_content_formatted = '\n'.join(libraries_formatted)
        data_lake_content_formatted = '\n'.join(data_lake_formatted)
        
        # Format the prompt with the appropriate values
        formatted_prompt = prompt_modifier.format(
            function_intro=function_intro,
            tool_desc=textify_api_dict(tool_desc) if isinstance(tool_desc, dict) else tool_desc,
            import_instruction=import_instruction,
            data_lake_path=self.path + '/data_lake',
            data_lake_intro=data_lake_intro,
            data_lake_content=data_lake_content_formatted,
            library_intro=library_intro,
            library_content_formatted=library_content_formatted
        )
        
        return formatted_prompt

    def configure(self, self_critic=False, test_time_scale_round=0):
        """
        Configure the agent with the initial system prompt and workflow.
        
        Args:
            self_critic: Whether to enable self-critic mode
            test_time_scale_round: Number of rounds for test time scaling
        """
        # Store self_critic for later use
        self.self_critic = self_critic
        
        # Get data lake content
        data_lake_path = self.path + '/data_lake'
        data_lake_content = glob.glob(data_lake_path + '/*')
        data_lake_items = [x.split('/')[-1] for x in data_lake_content]
        
        # Store data_lake_dict as instance variable for use in retrieval
        self.data_lake_dict = data_lake_dict
        # Store library_content_dict directly without library_content
        self.library_content_dict = library_content_dict

        # Prepare tool descriptions
        tool_desc = {i: [x for x in j if x['name']!='run_python_repl'] for i, j in self.module2api.items()}
        
        # Prepare data lake items with descriptions
        data_lake_with_desc = []
        for item in data_lake_items:
            description = self.data_lake_dict.get(item, f"Data lake item: {item}")
            data_lake_with_desc.append({"name": item, "description": description})
        
        # Generate the system prompt for initial configuration (is_retrieval=False)
        # Use library_content_dict.keys() instead of library_content
        self.system_prompt = self._generate_system_prompt(
            tool_desc=tool_desc,
            data_lake_content=data_lake_with_desc,
            library_content_list=list(self.library_content_dict.keys()),
            self_critic=self_critic,
            is_retrieval=False
        )
        
        # Define the nodes
        def generate(state: AgentState) -> AgentState:
            
            messages = [SystemMessage(content=self.system_prompt)] + state['messages']
            response = self.llm.invoke(messages)
            
            # Parse the response
            msg = str(response.content)

            # Check for incomplete tags and fix them
            if '<execute>' in msg and '</execute>' not in msg:
                msg += '</execute>'
            if '<solution>' in msg and '</solution>' not in msg:
                msg += '</solution>'
            if '<think>' in msg and '</think>' not in msg:
                msg += '</think>'
            
            think_match = re.search(r'<think>(.*?)</think>', msg, re.DOTALL)
            execute_match = re.search(r'<execute>(.*?)</execute>', msg, re.DOTALL)
            answer_match = re.search(r'<solution>(.*?)</solution>', msg, re.DOTALL)
            
            # Add the message to the state before checking for errors
            state['messages'].append(AIMessage(content=msg.strip()))
                        
            if answer_match:
                state['next_step'] = 'end'
            elif execute_match:
                state['next_step'] = 'execute'
            elif think_match:
                state['next_step'] = 'generate'
            else:
                print('parsing error...')
                # Check if we already added an error message to avoid infinite loops
                error_count = sum(1 for m in state['messages'] if isinstance(m, AIMessage) and "There are no tags" in m.content)
                
                if error_count >= 2:
                    # If we've already tried to correct the model twice, just end the conversation
                    print('Detected repeated parsing errors, ending conversation')
                    state['next_step'] = 'end'
                    # Add a final message explaining the termination
                    state['messages'].append(AIMessage(content="Execution terminated due to repeated parsing errors. Please check your input and try again."))
                else:
                    # Try to correct it
                    state['messages'].append(AIMessage(content="There are no tags (e.g. <execute><solution>). Please follow the instruction, fix and update."))
                    state['next_step'] = 'generate'
            return state
            
        def execute(state: AgentState) -> AgentState:
            last_message = state['messages'][-1].content
            # Only add the closing tag if it's not already there
            if '<execute>' in last_message and '</execute>' not in last_message:
                last_message += '</execute>'
            
            execute_match = re.search(r'<execute>(.*?)</execute>', last_message, re.DOTALL)
            if execute_match:
                code = execute_match.group(1)
                
                # Set timeout duration (10 minutes = 600 seconds)
                timeout = self.timeout_seconds
                
                # Check if the code is R code
                if code.strip().startswith("#!R") or code.strip().startswith("# R code") or code.strip().startswith("# R script"):
                    # Remove the R marker and run as R code
                    r_code = re.sub(r'^#!R|^# R code|^# R script', '', code, 1).strip()
                    result = run_with_timeout(run_r_code, [r_code], timeout=timeout)
                # Check if the code is a Bash script or CLI command
                elif code.strip().startswith("#!BASH") or code.strip().startswith("# Bash script") or code.strip().startswith("#!CLI"):
                    # Handle both Bash scripts and CLI commands with the same function
                    if code.strip().startswith("#!CLI"):
                        # For CLI commands, extract the command and run it as a simple bash script
                        cli_command = re.sub(r'^#!CLI', '', code, 1).strip()
                        # Remove any newlines to ensure it's a single command
                        cli_command = cli_command.replace('\n', ' ')
                        result = run_with_timeout(run_bash_script, [cli_command], timeout=timeout)
                    else:
                        # For Bash scripts, remove the marker and run as a bash script
                        bash_script = re.sub(r'^#!BASH|^# Bash script', '', code, 1).strip()
                        result = run_with_timeout(run_bash_script, [bash_script], timeout=timeout)
                # Otherwise, run as Python code
                else:
                    result = run_with_timeout(run_python_repl, [code], timeout=timeout)
                
                if len(result) > 10000:
                    result = 'The output is too long to be added to context. Here are the first 10K characters...\n'  + result[:10000]
                observation = f"\n<observation>{result}</observation>"
                state['messages'].append(AIMessage(content=observation.strip()))
                
            return state

        def routing_function(state: AgentState) -> Literal['execute', 'generate', 'end']:
            next_step = state.get('next_step')
            if next_step == 'execute':
                return 'execute'
            elif next_step == 'generate':
                return 'generate'
            elif next_step == 'end':
                return 'end'
            else:
                raise ValueError(f"Unexpected next_step: {next_step}")

        def routing_function_self_critic(state: AgentState) -> Literal['generate', 'end']:
            next_step = state.get('next_step')
            if next_step == 'generate':
                return 'generate'
            elif next_step == 'end':
                return 'end'
            else:
                raise ValueError(f"Unexpected next_step: {next_step}")

        def self_critic(state: AgentState) -> AgentState:
                       
            if self.critic_count < test_time_scale_round:
                # Generate feedback based on message history
                messages = state['messages']
                feedback_prompt = """
                Here is a reminder of what is the user requested: {user_task}
                Examine the previous executions, reaosning, and solutions. 
                Critic harshly on what could be improved? 
                Be specific and constructive. 
                Think hard what are missing to solve the task.
                No question asked, just feedbacks.
                """.format(user_task=self.user_task)
                feedback = self.llm.invoke(messages + [HumanMessage(content=feedback_prompt)])
                
                # Add feedback as a new message
                state['messages'].append(HumanMessage(content=f"Wait... this is not enough to solve the task. Here are some feedbacks for improvement:\n{feedback.content}"))
                self.critic_count += 1
                state['next_step'] = 'generate'
            else:
                state['next_step'] = 'end'
            
            return state
        # Create the workflow
        workflow = StateGraph(AgentState)

        # Add nodes
        workflow.add_node('generate', generate)
        workflow.add_node('execute', execute)

        if self_critic:
            workflow.add_node('self_critic', self_critic)
            # Add conditional edges
            workflow.add_conditional_edges('generate', routing_function, 
                                                    path_map={
                                                    'execute': 'execute',
                                                    'generate': 'generate',
                                                    'end': 'self_critic'
                                                })
            workflow.add_conditional_edges('self_critic', routing_function_self_critic,
                                                    path_map={
                                                    'generate': 'generate',
                                                    'end': END
                                                })
        else:
            # Add conditional edges
            workflow.add_conditional_edges('generate', routing_function, 
                                                    path_map={
                                                    'execute': 'execute',
                                                    'generate': 'generate',
                                                    'end': END
                                                })
        workflow.add_edge('execute', 'generate')
        workflow.add_edge(START, 'generate')

        # Compile the workflow
        self.app = workflow.compile()
        self.checkpointer = MemorySaver()
        self.app.checkpointer = self.checkpointer
        #display(Image(self.app.get_graph().draw_mermaid_png()))
        
    
    def go(self, prompt):
        """
        Execute the agent with the given prompt.
        
        Args:
            prompt: The user's query
        """
        self.critic_count = 0
        self.user_task = prompt
        
        if self.use_tool_retriever:            
            # Gather all available resources
            # 1. Tools from the registry
            all_tools = self.tool_registry.tools if hasattr(self, 'tool_registry') else []
            
            # 2. Data lake items with descriptions
            data_lake_path = self.path + '/data_lake'
            data_lake_content = glob.glob(data_lake_path + '/*')
            data_lake_items = [x.split('/')[-1] for x in data_lake_content]
            
            # Create data lake descriptions for retrieval
            data_lake_descriptions = []
            for item in data_lake_items:
                description = self.data_lake_dict.get(item, f"Data lake item: {item}")
                data_lake_descriptions.append({
                    "name": item,
                    "description": description
                })
            
            # 3. Libraries with descriptions - use library_content_dict directly
            library_descriptions = []
            for lib_name, lib_desc in self.library_content_dict.items():
                library_descriptions.append({
                    "name": lib_name,
                    "description": lib_desc
                })
            
            # Use retrieval to get relevant resources
            resources = {
                'tools': all_tools,
                'data_lake': data_lake_descriptions,
                'libraries': library_descriptions
            }
            
            # Use prompt-based retrieval with the agent's LLM
            selected_resources = self.retriever.prompt_based_retrieval(prompt, resources, llm=self.llm)
            print('Using prompt-based retrieval with the agent\'s LLM')
            
            # Log the selected resources
            print('Selected tools:')
            for tool in selected_resources['tools']:
                if isinstance(tool, dict):
                    print(f"- {tool.get('name', 'Unknown')}: {tool.get('description', '')}")
                else:
                    print(f"- {getattr(tool, 'name', str(tool))}: {getattr(tool, 'description', '')}")
                
            print('\nSelected data lake items:')
            for item in selected_resources['data_lake']:
                if isinstance(item, dict):
                    name = item.get('name', 'Unknown')
                    description = self.data_lake_dict.get(name, f"Data lake item: {name}")
                    print(f"- {name}: {description}")
                elif isinstance(item, str) and ': ' in item:
                    # If the item already has a description, print it as is
                    print(f"- {item}")
                else:
                    description = self.data_lake_dict.get(item, f"Data lake item: {item}")
                    print(f"- {item}: {description}")
                
            print('\nSelected libraries:')
            for lib in selected_resources['libraries']:
                if isinstance(lib, dict):
                    print(f"- {lib.get('name', 'Unknown')}: {lib.get('description', '')}")
                else:
                    print(f"- {lib}")
            
            # Extract the names from the selected resources for the system prompt
            selected_resources_names = {
                'tools': selected_resources['tools'],
                'data_lake': [],
                'libraries': [lib['name'] if isinstance(lib, dict) else lib for lib in selected_resources['libraries']]
            }
            
            # Process data lake items to extract just the names
            for item in selected_resources['data_lake']:
                if isinstance(item, dict):
                    selected_resources_names['data_lake'].append(item['name'])
                elif isinstance(item, str) and ': ' in item:
                    # If the item already has a description, extract just the name
                    name = item.split(': ')[0]
                    selected_resources_names['data_lake'].append(name)
                else:
                    selected_resources_names['data_lake'].append(item)
            
            # Update the system prompt with the selected resources
            self.update_system_prompt_with_selected_resources(selected_resources_names)

        inputs = {
            'messages': [HumanMessage(content=prompt)],
            'next_step': None
        }
        config = {"recursion_limit": 500, "configurable": {"thread_id": 42}}
        self.log = []
                
        for s in self.app.stream(inputs, stream_mode="values", config = config):
            message = s["messages"][-1]
            out = pretty_print(message)
            self.log.append(out)
            
        return self.log, message.content
        
    def update_system_prompt_with_selected_resources(self, selected_resources):
        """Update the system prompt with the selected resources."""
        # Extract tool descriptions for the selected tools
        tool_desc = {}
        for tool in selected_resources['tools']:
            # Get the module name from the tool
            if isinstance(tool, dict):
                module_name = tool.get('module', None)
                
                # If module is not specified, try to find it in the module2api
                if not module_name and hasattr(self, 'module2api'):
                    for mod, apis in self.module2api.items():
                        for api in apis:
                            if api.get('name') == tool.get('name'):
                                module_name = mod
                                # Update the tool with the module information
                                tool['module'] = module_name
                                break
                        if module_name:
                            break
                
                # If still not found, use a default
                if not module_name:
                    module_name = "biomni.tool.scRNA_tools"  # Default to scRNA_tools as a fallback
                    tool['module'] = module_name
            else:
                module_name = getattr(tool, 'module_name', None)
                
                # If module is not specified, try to find it in the module2api
                if not module_name and hasattr(self, 'module2api'):
                    tool_name = getattr(tool, 'name', str(tool))
                    for mod, apis in self.module2api.items():
                        for api in apis:
                            if api.get('name') == tool_name:
                                module_name = mod
                                # Set the module_name attribute
                                setattr(tool, 'module_name', module_name)
                                break
                        if module_name:
                            break
                
                # If still not found, use a default
                if not module_name:
                    module_name = "biomni.tool.scRNA_tools"  # Default to scRNA_tools as a fallback
                    setattr(tool, 'module_name', module_name)
                
            if module_name not in tool_desc:
                tool_desc[module_name] = []
                
            # Add the tool to the appropriate module
            if isinstance(tool, dict):
                # Ensure the module is included in the tool description
                if 'module' not in tool:
                    tool['module'] = module_name
                tool_desc[module_name].append(tool)
            else:
                # Convert tool object to dictionary
                tool_dict = {
                    "name": getattr(tool, 'name', str(tool)),
                    "description": getattr(tool, 'description', ''),
                    "parameters": getattr(tool, 'parameters', {}),
                    "module": module_name  # Explicitly include the module
                }
                tool_desc[module_name].append(tool_dict)
        
        # Generate the system prompt with the selected resources (is_retrieval=True)
        print(f"Data lake items for system prompt: {selected_resources['data_lake']}")
        
        # Prepare data lake items with descriptions
        data_lake_with_desc = []
        for item in selected_resources['data_lake']:
            description = self.data_lake_dict.get(item, f"Data lake item: {item}")
            data_lake_with_desc.append({"name": item, "description": description})
        
        self.system_prompt = self._generate_system_prompt(
            tool_desc=tool_desc,
            data_lake_content=data_lake_with_desc,
            library_content_list=selected_resources['libraries'],
            self_critic=getattr(self, 'self_critic', False),
            is_retrieval=True
        )
        
        # Print the raw system prompt for debugging
        print("\n" + "="*20 + " RAW SYSTEM PROMPT FROM AGENT " + "="*20)
        print(self.system_prompt)
        print("="*70 + "\n")

    def result_formatting(self, output_class, task_intention):

        self.format_check_prompt = ChatPromptTemplate.from_messages(
            [
                (
                    "system",
                    (
                        "You are evaluateGPT, tasked with extract and parse the task output based on the history of an agent. "
                        "Review the entire history of messages provided. "
                        "Here is the task output requirement: \n"
                        f"'{task_intention.replace('{', '{{').replace('}', '}}')}'.\n"
                    ),
                ),
                ("placeholder", "{messages}"),
            ]
        )

        checker_llm = self.format_check_prompt | self.llm.with_structured_output(output_class)
        result = checker_llm.invoke({ "messages": [("user", str(self.log))]}).dict()
        return result