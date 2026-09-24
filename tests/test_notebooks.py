"""Execute the teaching notebooks.

The notebooks are the part of this project students actually read, and they
are the part most likely to rot: they call the library by name, and a rename
or a changed default breaks them silently while every unit test still passes.
Running them is the only check that the documented usage is the real usage.

Each notebook is executed from its jupytext ``.py`` pair rather than from the
``.ipynb``. The ``.py`` is the source of truth -- AGENTS.md section 6 requires
all authoring to happen there -- and reading it means the test never depends
on committed outputs and never writes to the working tree.

Only the numbered teaching notebooks run. ``Scratchbook`` and ``tmp`` are
working files that change under the user's hands, and ``-Copy1`` duplicates
are JupyterLab's, gitignored; failing the suite on any of those would be
reporting on something other than the library.

These are the slowest tests in the suite, at roughly 20 seconds. Set
``MT_SKIP_NOTEBOOKS=1`` to skip them while iterating on something else.
"""

import os
import re
import unittest
from pathlib import Path

NOTEBOOK_DIR = Path(__file__).resolve().parent.parent / "notebooks"

#: Numbered teaching notebooks, excluding JupyterLab's "-Copy1" duplicates.
TEACHING_PATTERN = re.compile(r"^\d+_[A-Za-z0-9_]+\.py$")

SKIP = os.environ.get("MT_SKIP_NOTEBOOKS") == "1"


def teaching_notebooks() -> list[Path]:
    if not NOTEBOOK_DIR.is_dir():
        return []
    return sorted(path for path in NOTEBOOK_DIR.glob("*.py") if TEACHING_PATTERN.match(path.name))


#: A notebook that has not finished in this long has hung, not slowed down.
#: The slowest here runs in about 15 seconds, so this is generous by an order
#: of magnitude -- but it is the difference between a hung run costing four
#: minutes and costing fifteen.
CELL_TIMEOUT_S = 240

#: Kernels occasionally hang partway through execution: the process is alive
#: and idle, the client waits forever. Seen three times over the life of this
#: suite, on unchanged notebooks that pass on the next run. It is a jupyter
#: infrastructure flake rather than anything about the library, so one retry
#: is allowed -- but ONLY for kernel-level failures. A CellExecutionError is a
#: real bug in the notebook and is never retried.
KERNEL_FLAKE_RETRIES = 1


def execute(path: Path):
    """Run one notebook to completion, returning the executed notebook.

    Raises whatever the notebook raised, with the failing cell's traceback
    attached, which is what makes a failure here readable.
    """
    from nbclient.exceptions import CellTimeoutError, DeadKernelError

    last = None
    for attempt in range(KERNEL_FLAKE_RETRIES + 1):
        try:
            return _execute_once(path)
        except (CellTimeoutError, DeadKernelError) as exc:
            # The kernel, not the notebook. Retry once, then give up loudly.
            last = exc
            if attempt < KERNEL_FLAKE_RETRIES:
                continue
    raise AssertionError(
        f"{path.name}: the kernel hung or died on every attempt "
        f"({KERNEL_FLAKE_RETRIES + 1}). Last error: {type(last).__name__}. "
        "This is usually a jupyter flake rather than a notebook bug -- run "
        "the notebook by hand before assuming the library broke."
    ) from last


def _execute_once(path: Path):
    import jupytext
    from nbclient import NotebookClient

    notebook = jupytext.read(path)
    client = NotebookClient(
        notebook,
        timeout=CELL_TIMEOUT_S,
        # A kernel that fails to come up must fail, not hang. One run of this
        # suite sat on a kernel that never started; without this it waits
        # forever and looks like a slow test rather than a stuck one.
        startup_timeout=120,
        kernel_name=_kernel_name(),
        # Relative paths inside a notebook are written against its own
        # directory, so run it there.
        resources={"metadata": {"path": str(NOTEBOOK_DIR)}},
    )
    client.execute()
    return notebook


def _kernel_name() -> str:
    """Prefer the project's registered kernel, fall back to the plain one.

    setup.sh registers ``microbial-thermo``; a bare checkout may only have
    ``python3``, which in an activated environment is that environment's own.
    """
    try:
        from jupyter_client.kernelspec import KernelSpecManager

        available = KernelSpecManager().find_kernel_specs()
    except Exception:
        return "python3"
    return "microbial-thermo" if "microbial-thermo" in available else "python3"


class TestNotebooksArePaired(unittest.TestCase):
    """The jupytext discipline itself, which costs nothing to check."""

    def test_there_are_teaching_notebooks_to_run(self):
        """Guard against the discovery pattern silently matching nothing --
        which would make every other test here pass by doing nothing."""
        self.assertGreaterEqual(len(teaching_notebooks()), 3)

    def test_every_script_has_a_notebook_beside_it(self):
        for path in teaching_notebooks():
            with self.subTest(notebook=path.name):
                self.assertTrue(
                    path.with_suffix(".ipynb").exists(),
                    f"{path.name} has no paired .ipynb -- run jupytext --sync",
                )

    def test_every_script_is_in_percent_format(self):
        for path in teaching_notebooks():
            with self.subTest(notebook=path.name):
                text = path.read_text()
                self.assertIn("# %%", text, f"{path.name} is not in percent format")

    def test_jupyterlab_duplicates_are_not_collected(self):
        """They are the user's working copies and are gitignored."""
        names = [p.name for p in teaching_notebooks()]
        self.assertFalse([n for n in names if "-Copy" in n])

    def test_scratch_notebooks_are_not_collected(self):
        names = [p.name for p in teaching_notebooks()]
        self.assertNotIn("Scratchbook.py", names)
        self.assertNotIn("tmp.py", names)


@unittest.skipIf(SKIP, "MT_SKIP_NOTEBOOKS=1")
class TestNotebooksExecute(unittest.TestCase):
    """One test method per notebook, generated below.

    Generated rather than looped inside a single test so that a failure names
    the notebook that failed, and so one broken notebook does not hide the
    others.
    """


def _add_execution_tests():
    for path in teaching_notebooks():
        name = f"test_{path.stem}_runs_clean"

        def check(self, path=path):
            notebook = execute(path)
            errors = [
                output
                for cell in notebook.cells
                for output in cell.get("outputs", [])
                if output.get("output_type") == "error"
            ]
            self.assertEqual(
                errors,
                [],
                f"{path.name} raised: "
                + "; ".join(f"{e.get('ename')}: {e.get('evalue')}" for e in errors),
            )

        check.__name__ = name
        check.__doc__ = f"Execute {path.name} end to end."
        setattr(TestNotebooksExecute, name, check)


_add_execution_tests()


if __name__ == "__main__":
    unittest.main()
