# Universal AI Agent Blueprint (AGENTS.md)

This document is the absolute source of truth and operational blueprint for all AI agents and autonomous coding harnesses operating within this repository. You MUST adhere strictly to these constraints, conventions, and software carpentry guidelines to ensure scientifically reproducible and maintainable code.

## 1. Core Directives & Agent Execution

* **Read Before Writing:** Always read the current state of a file, module, or README before proposing or applying edits. Never assume the structure of the codebase.
* **Objective Truth over Agreement (Anti-Sycophancy):** Do not blindly agree with the user if they propose a flawed architectural choice, insecure code, or an anti-pattern. Push back logically and propose the scientifically and computationally optimal solution.
* **No Unnecessary Apologies:** Do not waste token space apologizing. Acknowledge corrections and output the fixed code or plan immediately.
* **Directness:** Be concise, direct, and authoritative in your domain knowledge. Provide facts over filler.
* **No Destructive Actions Without Consent:** Never execute `rm -rf`, `drop table`, or equivalent destructive commands without explicit, unambiguous user confirmation.
* **Secret Protection:** Never log, print, or commit API keys, tokens, or passwords. Always respect `.env.template` and strictly ignore actual `.env` files.

## 2. Test-Driven Dialogue & Clarification

* **The Zero-Assumption Policy:** Do not write core logic or add new functionality based on ambiguous requests. If a prompt lacks specific constraints, edge-case definitions, or architectural clarity, you should pause execution and ask the user formatted, numbered questions to resolve the ambiguity.
* **Demand Ground-Truth Examples:** Before generating functional scripts, data pipelines, or algorithms, you should ask the user to provide known-correct example inputs and their exact expected outputs. Use these to write tests. Do not guess data structures.

## 3. Scientific Rigor & Software Carpentry

* **Keep It Simple, Stupid (KISS):** Avoid over-engineering. Write clear, modular, and readable code over clever, complex one-liners. Junior researchers must be able to read, audit, and understand the codebase.
* **Reproducibility First:** Scientific code must be perfectly reproducible. Ensure all random seeds are settable, data pipelines are deterministic, and outputs can be perfectly replicated on a standard Linux environment.
* **Mathematical & Chemical Formatting:** When documenting scientific processes in markdown or docstrings, utilize LaTeX markup for rendering (e.g., $$CH_4 + 2O_2 \rightarrow CO_2 + 2H_2O$$).

## 4. Environment & Dependency Management

* **Mamba Exclusivity:** This project uses `mamba` exclusively for environment management and dependency resolution. ALWAYS execute tasks, tests, and scripts from within the active mamba environment.
* **Mamba Initialization:** When opening a new shell or persistent terminal, AI agents must initialize mamba using the following snippet to ensure `mamba` is in the PATH:
  ```bash
  export MAMBA_EXE='/home/$USER/miniforge3/bin/mamba'
  export MAMBA_ROOT_PREFIX='/home/$USER/miniforge3'
  __mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
  if [ $? -eq 0 ]; then
      eval "$__mamba_setup"
  else
      alias mamba="$MAMBA_EXE"
  fi
  unset __mamba_setup
  ```
* **Adding Dependencies:** Use `mamba install <package>`. You may only use `pip install` as a last resort if absolutely no conda/mamba-forge package is available.
* **Tracking Dependencies:** If you add a new library, you MUST immediately update the `environment.yaml` and/or `pyproject.toml` files.

## 5. Code Architecture & CLI

* **Library-First Design:** All Python scripts must be structured as importable libraries. Use `__init__.py` files appropriately to create a cohesive Python module. Do not write standalone, procedural "spaghetti" scripts.
* **Execution Block:** When creating runnable scripts, you MUST include a standard execution block at the bottom of the file to encourage proper modularity, and also use the `click` library as recommended. Core logic must be decoupled from the CLI argument parsing so it can be safely imported and executed natively within Jupyter without triggering the CLI:
  ```python
  if __name__ == "__main__":
      main()
  ```
* **Command Line Interfaces (CLI):** Expose functional entry points using the `click` library. Ensure the CLI is intuitive, clearly parameterized, and self-documenting.
* **Mandatory Versioning (`--version`):** Every CLI MUST include a `--version` argument. This must print the current version of the script AND the versions of its primary computational dependencies (e.g., `pandas`, `matplotlib`, `numpy`, `scipy`) to guarantee environment traceability.
* **Preferred Libraries:** Default to `pandas` for tabular data manipulation and `matplotlib` or `seaborn` for plotting, unless specifically instructed otherwise.

## 6. Jupyter Notebook Integration & Sync

* **Jupytext Synchronization:** All Jupyter notebooks (`.ipynb`) must be paired with a corresponding Python script (`.py`) using `jupytext`.
* **Percent Format:** The paired `.py` files must strictly use the `percent` format, demarcating Jupyter cells with `# %%` comments to ensure a clear mapping between scripts and notebooks.
* **Develop Only in `.py`:** AI agents MUST perform all code authoring, refactoring, and debugging in the plain-text `.py` file. Never attempt to edit the raw JSON structure of the `.ipynb` file directly.
* **Sync Before and After Edits:** Do not assume the state of the codebase. Always verify and sync the files using `jupytext --sync <notebook.ipynb>` before making edits, and sync again after completing the edits to ensure the notebook reflects the updated `.py` script.
* **Notebook Execution:** Syncing the code structure is sufficient for updating the notebook. Ask the user if the notebook should be explicitly executed to regenerate outputs (e.g., plots or dataframe renders) after a sync.

## 7. Pre-Commit Checklist: Testing, Linting & Docs

Before proposing a commit or finalizing a task, you MUST complete this checklist:

1. **Testing:** Run the test suite using the standard Python `unittest` library. Write new tests for any newly added logic based on the user's ground-truth examples and exspected input and outputs.
2. **Linting:** Run standard linters and formatters. Always use **ruff** for linting and formatting. code formatting via `ruff format`, and linting via `ruff check`. Fix all warnings and errors.
3. **Documentation Sync (README.md):** The `README.md` must be a complete, standalone guide. Whenever you update code or CLI arguments, you MUST update the README to include:
    * **Installation Instructions:** Always include exact `mamba install` instructions to recreate the environment.
    * **Jupyter Notebook Usage:** Include instructions on how to use the Jupyter notebooks, explicitly showing how to register the `mamba`-built environment with Jupyter (e.g., using `python -m ipykernel install --user --name <env_name>`).
    * **Library & API Usage:** Instructions for using the tool as an importable Python library.
    * **Current CLI Usage:** Run the CLI with `--help` and explicitly copy/paste the exact, updated output into the README as well as providing examples. 
    * **Version Tracking:** Run the CLI with `--version` and explicitly copy/paste the output into the README.
    * **Spellcheck comments and README.md:** For any file the user has edited, spellcheck the comments and text in README.md or other files and correct, asking questions if you are unsure. 

## 8. Version Control & Git Flow

* **Branch Naming:** All AI-generated work must be done on an isolated branch following the exact format: `ai/<model>/<ticket-or-feature>` (e.g., `ai/gemini/refactor-data-loader`).
* **Workflow Constraints:** Always verify the current branch (`git branch --show-current`) before executing edits. If on `main` or a non-compliant branch, immediately create and checkout the appropriate `ai/` branch. NEVER commit directly to `main`.

## 9. Multimodal Ingestion
* **Leverage Vision Context:** When generating or debugging data visualizations or plots, proactively request the user to provide screenshot outputs so you can visually verify correctness.
* **Complex Ingestion:** Utilize native multimodal capabilities to parse user-provided architectural diagrams, PDFs, and structural images, translating them directly into the target architecture.
* **MCP Integration:** An `mcp.json` file is provisioned in the project root containing domain-specific MCP servers (Jupyter, SciAgent, BioinfoMCP, scientific-visualization). Rely on these servers—especially the Jupyter MCP—for stateful data exploration, bioinformatics workflows, and iterative problem-solving natively.

