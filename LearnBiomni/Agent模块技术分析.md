# Biomni Agent 模块技术分析

## 1. 概述

Biomni Agent 模块是整个 Biomni 系统的核心，负责接收用户查询、规划执行路径、调用适当的工具、处理执行结果并生成最终答案。该模块实现了三种不同的代理模式：

1. **A1 代理**：主要代理实现，基于 LangGraph 的工作流系统
2. **ReAct 代理**：基于推理-行动-观察循环的代理实现
3. **QA LLM 代理**：简单的问答代理，直接使用 LLM 回答问题

本文档将详细分析 Agent 模块的实现方法和工具调用机制。

## 2. 核心类和组件

### 2.1 A1 代理 (`a1.py`)

A1 是 Biomni 的主要代理实现，具有以下核心组件：

#### 状态管理
```python
class AgentState(TypedDict):
    messages: list[BaseMessage]  # 消息历史
    next_step: str | None        # 下一步操作
```

#### 初始化
```python
def __init__(
    self,
    path="./data",
    llm="claude-sonnet-4-20250514",
    use_tool_retriever=True,
    timeout_seconds=600,
):
```

A1 代理在初始化时会：
1. 检查并创建数据目录
2. 下载必要的数据（如果不存在）
3. 加载工具注册表
4. 初始化 LLM
5. 设置工具检索器（如果启用）
6. 设置超时参数

### 2.2 ReAct 代理 (`react.py`)

ReAct 代理基于推理-行动-观察循环，具有以下特点：

#### 状态管理
```python
class AgentState(TypedDict):
    """The state of the agent."""
    # add_messages is a reducer that combines message sequences
    messages: Annotated[Sequence[BaseMessage], add_messages]
```

#### 初始化
```python
def __init__(
    self,
    path="./data",
    llm="claude-3-7-sonnet-latest",
    use_tool_retriever=False,
    timeout_seconds=600,
):
```

### 2.3 QA LLM 代理 (`qa_llm.py`)

最简单的代理实现，直接使用 LLM 回答问题：

```python
def __init__(self, path="./data", llm="claude-3-haiku-20240307", lab_bench_reproduce=False):
    self.path = path
    self.llm = get_llm(llm)
    # 设置提示词修饰符
    if lab_bench_reproduce:
        self.prompt_modifier = """
The following is a multiple choice question about biology.
Please answer by responding with the letter of the correct answer.

Think step by step. \n
        """
    else:
        self.prompt_modifier = ""
    self.log = []
```

## 3. 工作流程分析

### 3.1 A1 代理工作流程

A1 代理基于 LangGraph 实现了一个复杂的工作流系统，主要包含以下步骤：

#### 1. 配置阶段 (`configure`)
```python
def configure(self, self_critic=False, test_time_scale_round=0):
```

在配置阶段，A1 代理会：
1. 加载数据湖内容
2. 准备工具描述
3. 生成系统提示词
4. 设置 LangGraph 工作流

#### 2. 工作流节点

A1 代理定义了三个主要的工作流节点：

##### a. 生成节点 (`generate`)
```python
def generate(state: AgentState) -> AgentState:
```
该节点负责：
- 调用 LLM 生成响应
- 解析响应中的标签（`<execute>`, `<solution>`, `<think>`）
- 确定下一步操作

##### b. 执行节点 (`execute`)
```python
def execute(state: AgentState) -> AgentState:
```
该节点负责：
- 解析执行代码
- 根据代码类型选择执行方式（Python、R、Bash）
- 设置超时控制
- 将执行结果添加到状态

##### c. 路由函数 (`routing_function`)
```python
def routing_function(state: AgentState) -> Literal["execute", "generate", "end"]:
```
该函数决定下一步应该执行哪个节点。

#### 3. 执行阶段 (`go`)
```python
def go(self, prompt):
```

在执行阶段，A1 代理会：
1. 使用工具检索器选择相关工具（如果启用）
2. 更新系统提示词
3. 初始化状态并启动工作流
4. 记录执行日志
5. 返回执行结果

### 3.2 ReAct 代理工作流程

ReAct 代理基于 LangGraph 实现了一个简化的工作流：

#### 1. 配置阶段 (`configure`)
```python
def configure(self, plan=False, reflect=False, data_lake=False, react_code_search=False, library_access=False):
```

根据参数配置不同的提示词修饰符，可以启用：
- 计划生成
- 反思机制
- 数据湖访问
- 代码搜索
- 库访问

#### 2. 自定义 ReAct 代理创建
```python
def _create_custom_react_agent(self, llm, tools, prompt):
```

该方法创建一个自定义的 ReAct 代理，包含：
- 模型调用节点 (`call_model`)
- 工具执行节点 (`tool_node`)
- 条件边缘判断 (`should_continue`)

#### 3. 执行阶段 (`go`)
```python
def go(self, prompt):
```

执行过程与 A1 类似，但更专注于工具的循环调用。

### 3.3 QA LLM 代理工作流程

QA LLM 代理是最简单的实现：

```python
def go(self, input):
    self.log = []
    self.log.append(("user", input))
    message = self.llm.invoke(self.prompt_modifier + input)
    self.log.append(("assistant", message.content))
    return [message.content], message.content
```

直接将用户输入传递给 LLM，并返回结果。

## 4. 工具调用机制

### 4.1 工具注册和管理

#### 工具注册表 (`ToolRegistry`)
```python
self.tool_registry = ToolRegistry(module2api)
```

工具注册表负责：
- 存储所有可用工具
- 提供工具查询功能
- 管理工具 ID 和名称

#### 添加自定义工具
```python
def add_tool(self, api):
```

该方法允许动态添加自定义工具，步骤包括：
1. 获取函数源代码
2. 生成 API schema
3. 验证和增强 schema
4. 将工具添加到注册表
5. 更新 `module2api` 结构

### 4.2 工具检索机制

#### 工具检索器 (`ToolRetriever`)
```python
self.retriever = ToolRetriever()
```

工具检索器通过提示词方法选择最相关的工具：

```python
selected_resources = self.retriever.prompt_based_retrieval(prompt, resources, llm=self.llm)
```

检索过程：
1. 构建提示词描述所有可用资源
2. 让 LLM 选择最相关的资源
3. 解析 LLM 响应获取选定的工具、数据和库

### 4.3 工具执行机制

#### A1 代理中的执行
```python
def execute(state: AgentState) -> AgentState:
```

A1 代理根据代码类型选择执行方式：
- Python 代码：使用 `run_python_repl`
- R 代码：使用 `run_r_code`
- Bash 脚本：使用 `run_bash_script`

所有执行都有超时控制：
```python
result = run_with_timeout(run_python_repl, [code], timeout=timeout)
```

#### ReAct 代理中的执行
```python
def tool_node(state: AgentState):
```

ReAct 代理使用 LangChain 工具接口执行工具：
```python
tool_result = tools_by_name[tool_call["name"]].invoke(tool_call["args"])
```

### 4.4 超时控制机制

#### A1 代理的超时控制
```python
result = run_with_timeout(run_python_repl, [code], timeout=timeout)
```

使用 `run_with_timeout` 函数控制执行时间。

#### ReAct 代理的超时控制
```python
def _add_timeout_to_tools(self, tools):
```

ReAct 代理使用多进程方法为每个工具添加超时控制：
1. 创建一个单独的进程执行工具
2. 设置超时时间
3. 如果超时，强制终止进程
4. 返回超时错误消息

## 5. 提示词工程

### 5.1 系统提示词生成

```python
def _generate_system_prompt(
    self,
    tool_desc,
    data_lake_content,
    library_content_list,
    self_critic=False,
    is_retrieval=False,
    custom_tools=None,
    custom_data=None,
    custom_software=None,
):
```

系统提示词包含：
- 代理角色描述
- 可用工具列表及描述
- 数据湖内容描述
- 可用库描述
- 输出格式指导
- 思考-执行-解决方案框架

### 5.2 ReAct 代理的提示词修饰符

ReAct 代理可以根据不同参数配置不同的提示词修饰符：

```python
if (not plan) and (not reflect):
    prompt_modifier = """You are a helpful biologist and expert geneticist."""
elif plan and (not reflect):
    prompt_modifier = """You are a helpful biologist and expert geneticist.
    Given the question from the user,
    - First, come up with a high level plan...
    """
```

## 6. 自我批评机制

A1 代理实现了自我批评机制：

```python
def self_critic(state: AgentState) -> AgentState:
```

自我批评流程：
1. 分析当前状态和结果
2. 识别问题和改进点
3. 调整执行策略
4. 重新生成响应

## 7. 结果格式化

所有代理都实现了结果格式化方法：

```python
def result_formatting(self, output_class, task_intention):
```

该方法使用 LLM 将代理输出解析为指定的结构化格式。

## 8. 技术架构对比

### 8.1 三种代理模式对比

| 特性 | A1 代理 | ReAct 代理 | QA LLM 代理 |
|------|---------|------------|-------------|
| 复杂度 | 高 | 中 | 低 |
| 工作流 | 基于 LangGraph | 基于 LangGraph | 简单调用 |
| 标签解析 | `<execute>`, `<solution>`, `<think>` | 工具调用 | 无 |
| 工具检索 | 支持 | 支持 | 不支持 |
| 自我批评 | 支持 | 不支持 | 不支持 |
| 超时控制 | 支持 | 支持 | 不支持 |
| 适用场景 | 复杂任务 | 中等复杂任务 | 简单问答 |

### 8.2 工具调用方法对比

| 特性 | A1 代理 | ReAct 代理 |
|------|---------|------------|
| 调用方式 | 基于标签解析 | 基于工具对象 |
| 代码执行 | 支持多语言 (Python, R, Bash) | 主要支持 Python |
| 参数传递 | 通过代码字符串 | 通过结构化参数 |
| 错误处理 | 观察结果中处理 | 工具节点中处理 |
| 超时实现 | 函数级别 | 工具级别 |

## 9. 总结与最佳实践

### 9.1 实现亮点

1. **模块化设计**：三种代理模式满足不同复杂度需求
2. **灵活的工具管理**：动态注册和检索工具
3. **强大的执行控制**：超时机制和错误处理
4. **自适应提示词**：根据任务和可用资源调整系统提示词
5. **自我批评机制**：提高代理执行质量

### 9.2 工具调用最佳实践

1. **使用工具检索**：根据任务选择最相关的工具
2. **设置合理超时**：防止工具执行时间过长
3. **结构化输出**：使用格式化方法处理代理输出
4. **错误处理**：捕获并处理工具执行错误
5. **自定义工具**：使用 `add_tool` 方法扩展功能

### 9.3 扩展建议

1. **添加新工具**：使用 `add_tool` 方法添加自定义工具
2. **自定义数据**：使用 `add_data` 方法添加自定义数据
3. **自定义软件**：使用 `add_software` 方法添加自定义软件
4. **调整提示词**：修改系统提示词以适应特定任务
5. **配置工作流**：根据需要调整工作流节点和路由 