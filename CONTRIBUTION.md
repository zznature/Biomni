# Contributing to Biomni

Thank you for your interest in contributing to Biomni! We're building the infrastructure layer for biomedical AI agents, and we welcome contributions from the community. Contributors with significant contributions will be invited to co-author publications in top-tier journals and conferences.

## Getting Started

Before contributing, please ensure you:
- Have tested your changes locally
- Follow the existing code style and conventions
- Include appropriate documentation

## Types of Contributions

### 🛠️ Adding a New Tool

Tools are implemented as Python functions in `biomni/tool/XXX.py`, organized by subject area.

**Steps:**
1. **Implement and test** your function locally
2. **Choose the appropriate subject** category
3. **Create a tool description** in `biomni/tool/tool_description/XXX.py` following the existing format

   *Tip: Use this helper to auto-generate descriptions:*
   ```python
   from biomni.utils import function_to_api_schema
   from biomni.llm import get_llm
   
   llm = get_llm('claude-sonnet-4-20250514')
   desc = function_to_api_schema(function_code, llm)
   ```

4. **Create a test prompt** that uses your tool and verify the agent works correctly
5. **Submit a pull request** for review, don't forget to include your test prompt as well

### 📊 Adding New Data

**Requirements:**
- Data must not overlap with existing datasets
- Must have proper redistribution rights
- Accessible download link required

**Steps:**
1. **Verify uniqueness** - ensure no overlap with existing data
2. **Prepare download link** with verified redistribution rights
3. **Add entry** to `data_lake_dict` in `biomni/env_desc.py`
4. **Submit a pull request** with the download link

### 💻 Adding New Software

**Steps:**
1. **Test locally** to ensure no conflicts with existing environments
2. **Create installation script** as a bash file
3. **Add entry** to `library_content_dict` in `biomni/env_desc.py`
4. **Submit a pull request** including:
   - Installation bash script
   - Screenshot demonstrating no environment conflicts

### 🎯 Adding a New Benchmark

Create benchmarks in the `biomni/task/` folder.

**Required implementation:**
```python
class YourBenchmark:
    def __init__(self):
        # Initialize benchmark
        pass
    
    def __len__(self):
        # Return dataset size
        pass
    
    def get_example(self, index):
        # Return dataset item at index
        pass
    
    def evaluate(self):
        # Evaluation logic (flexible input format)
        pass
    
    def output_class(self):
        # Return expected agent output format
        pass
```

**Steps:**
1. **Create benchmark file** in `biomni/task/[benchmark_name].py`
2. **Implement required methods** as shown above
3. **Provide data download link** for associated datasets
4. **Submit a pull request**

### 🐛 Bug Fixes & Enhancements

We welcome all bug fixes and enhancements to the existing codebase!

**Guidelines:**
- Clearly describe the issue or enhancement
- Include tests when applicable
- Follow existing code patterns
- Update documentation if needed

## Submission Process

1. **Fork** the repository
2. **Create a feature branch** from `main`
3. **Make your changes** following the guidelines above
4. **Test thoroughly** in your local environment
5. **Submit a pull request** with a clear description

## Review Process

The Biomni team will review all pull requests promptly. We may request changes or provide feedback to ensure code quality and consistency.

## Questions?

If you have questions about contributing, please open an issue or reach out to the maintainers.

---

*Together, let's build the future of biomedical AI agents!*
