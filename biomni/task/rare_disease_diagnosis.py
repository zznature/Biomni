import sys, json, numpy as np, pandas as pd, requests
from pydantic import BaseModel, Field

from biomni.task.base_task import base_task

class rare_disease_diagnosis(base_task):
    def __init__(self, path = './data', eval_llm = 'claude-3-5-sonnet-20241022', num_samples = None):
        self.llm = get_llm(eval_llm)
        data_path = path + '/rare_disease_diagnosis/mygene.json'
        data = []
        with open(data_path, "r") as file:
            for line in file:
                data.append(json.loads(line))

        if num_samples is None:
            self.data = pd.DataFrame(data)
        else:
            self.data = pd.DataFrame(data)[:num_samples]

        # Ensure the data contains all necessary columns
        required_columns = ['id', 'positive_phenotypes', 'all_candidate_genes', 'omim', 'disease_name', 'orpha_id']
        for col in required_columns:
            if col not in self.data.columns:
                raise ValueError(f"Dataset is missing required column: {col}")

        self.query = []
        self.answer = []
        for _, row in self.data.iterrows():
            phenotypes = row['positive_phenotypes']
            candidate_genes = row['all_candidate_genes']
            disease_name = row['disease_name']
            omim_id = row['omim']

            self.query.append({
                "phenotypes": phenotypes,
                "candidate_genes": candidate_genes
            })
            self.answer.append({
                "disease_name": disease_name,
                "OMIM_ID": omim_id
            })

        self.task_description = """
Task: given a patient's phenotypes and a list of candidate genes, diagnose the rare disease that the patient has.
Phenotypes: {phenotype_list}
Candidate genes: {candidate_genes}

Output format: {{'disease_name': XXX, 'OMIM_ID': XXX}}
        """

        self.completion_checker = """
Given an answer and a solution, check if the answer is correct.

Answer: {answer}
Solution: {solution}

Return 'task completed' if the answer is correct, and 'task not completed' otherwise.
        """

    def __len__(self):
        return len(self.query)
        
    def get_example(self, index=None):
        if index is None:
            index = np.random.randint(len(self.query))
        
        q = self.query[index]
        a = self.answer[index]

        prompt = self.task_description.format(
            phenotype_list=', '.join(q['phenotypes']),
            candidate_genes=q['candidate_genes']
        )
            
        return {"prompt": prompt, "answer": a}

    def split(self, ratio=0.8, seed=42):
        np.random.seed(seed)
        indices = np.arange(len(self.query))
        np.random.shuffle(indices)
        split_idx = int(ratio * len(self.query))
        train_indices = indices[:split_idx]
        val_indices = indices[split_idx:]
        return train_indices, val_indices

    def reward(self, input, output):
        answer = self.get_example(input)['answer']
        return 1 if output.OMIM_ID == answer.OMIM_ID else 0

    def get_iterator(self):
        for i in range(len(self.query)):
            yield self.get_example(i)

    def output_class(self):
        from pydantic import BaseModel, Field
        from typing import Optional
        
        class DiagnosisOutput(BaseModel):
            """A diagnosis for a rare disease."""
            
            disease_name: Optional[str] = Field(
                description="The name of the diagnosed rare disease, e.g., 'Marfan Syndrome'"
            )
            OMIM_ID: Optional[str] = Field(
                description="The OMIM ID of the diagnosed disease, e.g., '154700'"
            )
        
        return DiagnosisOutput

    def evaluate(self, response, ground_truth = None):
        from sklearn.metrics import accuracy_score
        if ground_truth is None:
            ground_truth = self.answer
        predicted = response
        correct = []
        results = []
        for pred, gt in zip(predicted, ground_truth):
            # Use the LLM-based completion checker to verify each prediction
            check_prompt = self.completion_checker.format(
                answer=json.dumps(pred),
                solution=json.dumps(gt)
            )
            # Assuming an LLM API call here; replace with the actual implementation
            result = self.call_llm_to_check(check_prompt)
            correct.append(result == 'task completed')
            results.append(result)
        
        accuracy = accuracy_score([1] * len(correct), correct)
        return {
            'completion_rate': accuracy,
            'num_of_tasks_completed': sum(correct),
            'num_of_total_tasks': len(correct),
            'results': results
        }

    def call_llm_to_check(self, prompt):
        class output_format(BaseModel):
            """Parse if the task is completed or not."""
            completion_status: str = Field(
                description="""'task completed' if the answer shows that it is completed, and 'task not completed' otherwise."""
            )
        
        output_parser = self.llm.with_structured_output(output_format)
        return self.llm.invoke(prompt).completion_status