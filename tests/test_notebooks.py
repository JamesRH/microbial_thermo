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

These are the slowest tests in the suite: about two minutes of a ~225 s run
plugged in, with notebook 13 the largest at ~35 s and 02 next at ~20 s, and
over twice that on battery. Set ``MT_SKIP_NOTEBOOKS=1`` to skip them while
iterating on something else.
"""

import contextlib
import os
import re
import threading
import time
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


#: A cell that has not finished in this long has hung, not slowed down. The
#: slowest whole notebook takes about 35 seconds, so this is generous by an
#: order of magnitude. An idle hang is now caught sooner by the stall
#: watchdog below; this remains the backstop for a cell that is busy forever.
CELL_TIMEOUT_S = 240

#: Kernels occasionally hang partway through execution: the process is alive
#: and idle, the client waits forever. Seen three times over the life of this
#: suite, on unchanged notebooks that pass on the next run. It is a jupyter
#: infrastructure flake rather than anything about the library, so one retry
#: is allowed -- but ONLY for kernel-level failures. A CellExecutionError is a
#: real bug in the notebook and is never retried.
KERNEL_FLAKE_RETRIES = 1

#: How long a running cell may go without using CPU before its kernel is
#: declared stalled and killed. A hung kernel is not slow, it is idle: one
#: sitting in ``time.sleep`` was measured at about 1 ms of CPU per second,
#: where a cell doing real work uses a whole second per second. So a stall
#: shows within seconds, and this window only has to be long enough not to
#: mistake a brief pause for one. It replaces the 240 s cell timeout as what
#: catches a hang.
STALL_S = 30.0

#: How often the watchdog samples the kernel's CPU time.
STALL_POLL_S = 5.0

#: A window counts as busy if the kernel used at least this fraction of it
#: on the CPU -- 1.5 s of CPU in 30 s. An idle kernel uses about 0.1 %.
STALL_MIN_BUSY = 0.05


def execute(path: Path, stall_s: float = STALL_S):
    """Run one notebook to completion, returning the executed notebook.

    Raises whatever the notebook raised, with the failing cell's traceback
    attached, which is what makes a failure here readable.
    """
    from nbclient.exceptions import CellTimeoutError, DeadKernelError

    last = None
    for attempt in range(KERNEL_FLAKE_RETRIES + 1):
        try:
            return _execute_once(path, stall_s=stall_s)
        except (CellTimeoutError, DeadKernelError) as exc:
            # The kernel, not the notebook. Retry once, then give up loudly.
            last = exc
            if attempt < KERNEL_FLAKE_RETRIES:
                continue
    raise AssertionError(
        f"{path.name}: the kernel hung or died on every attempt "
        f"({KERNEL_FLAKE_RETRIES + 1}). Last error: {type(last).__name__}: "
        f"{last}. This is usually a jupyter flake rather than a notebook bug "
        "-- run the notebook by hand before assuming the library broke."
    ) from last


def _execute_once(path: Path, stall_s: float = STALL_S):
    import jupytext

    return run_notebook(jupytext.read(path), stall_s=stall_s)


def run_notebook(notebook, stall_s: float = STALL_S, poll_s: float = STALL_POLL_S):
    """Execute an in-memory notebook under the stall watchdog, in place.

    Raises ``DeadKernelError`` whose message starts "Kernel stalled" when the
    watchdog killed the kernel, so a hang is not reported as a crash.
    """
    from nbclient import NotebookClient
    from nbclient.exceptions import DeadKernelError

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
    watchdog = StallWatchdog(client, stall_s=stall_s, poll_s=poll_s)
    try:
        client.execute()
    except DeadKernelError as exc:
        if watchdog.fired:
            raise DeadKernelError(f"Kernel stalled and was killed: {watchdog.fired}") from exc
        raise
    finally:
        watchdog.stop()
    return notebook


class StallWatchdog:
    """Kill a kernel that is running a cell but has stopped using the CPU.

    A thread samples the CPU time of the kernel and its children. If a cell
    has been running for ``stall_s`` and the kernel used less than
    ``min_busy`` of that on the CPU, the kernel's process tree is killed.
    nbclient polls the kernel once a second, sees it dead, and raises
    ``DeadKernelError``, which :func:`execute` already retries.

    Which hooks mark a running cell matters. ``on_cell_complete`` fires as
    soon as the execute request is *sent*, not when the cell finishes, so a
    cell runs from ``on_cell_execute`` to ``on_cell_executed`` (which fires
    after the reply). A hung client may never reach ``on_cell_executed`` --
    which is the case to catch -- so the clock starts on ``on_cell_execute``
    and only that hook resets it.
    """

    def __init__(
        self,
        client,
        stall_s: float = STALL_S,
        poll_s: float = STALL_POLL_S,
        min_busy: float = STALL_MIN_BUSY,
    ):
        self.client = client
        self.stall_s = stall_s
        self.poll_s = poll_s
        self.min_busy = min_busy
        #: Why the kernel was killed, or None if it never was.
        self.fired = None
        self._lock = threading.Lock()
        self._stopped = threading.Event()
        self._thread = None
        self._process = None
        self._cell = None  # index of the running cell, None between cells
        self._mark = (0.0, 0.0)  # (wall clock, CPU seconds) at window start
        client.on_cell_execute = self._cell_started
        client.on_cell_executed = self._cell_finished

    def stop(self):
        self._stopped.set()

    def _cell_started(self, cell, cell_index):
        import psutil

        if self._thread is None:
            # The kernel exists by the time the first cell is sent.
            self._process = psutil.Process(self.client.km.provisioner.process.pid)
            self._thread = threading.Thread(target=self._watch, daemon=True)
            self._thread.start()
        with self._lock:
            self._mark = (time.monotonic(), self._cpu_seconds())
            self._cell = cell_index

    def _cell_finished(self, cell, cell_index, execute_reply):
        with self._lock:
            self._cell = None

    def _cpu_seconds(self) -> float:
        import psutil

        try:
            processes = [self._process, *self._process.children(recursive=True)]
        except psutil.NoSuchProcess:
            return 0.0
        total = 0.0
        for process in processes:
            try:
                times = process.cpu_times()
            except psutil.NoSuchProcess:
                continue
            total += times.user + times.system + times.children_user + times.children_system
        return total

    def _watch(self):
        while not self._stopped.wait(self.poll_s):
            with self._lock:
                if self._cell is None:
                    continue
                started, cpu_at_start = self._mark
                now = time.monotonic()
                used = self._cpu_seconds() - cpu_at_start
                if used >= self.min_busy * self.stall_s:
                    # Busy: start a new window from here.
                    self._mark = (now, cpu_at_start + used)
                    continue
                if now - started < self.stall_s:
                    continue
                self.fired = f"cell {self._cell} used {used:.2f} s of CPU in {now - started:.0f} s"
            self._kill()
            return

    def _kill(self):
        import psutil

        try:
            processes = [*self._process.children(recursive=True), self._process]
        except psutil.NoSuchProcess:
            return
        for process in processes:
            with contextlib.suppress(psutil.NoSuchProcess):
                process.kill()


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


class TestNotebooksAreInSync(unittest.TestCase):
    """The ``.py`` and the ``.ipynb`` must hold the same cells.

    Nothing else here would catch drift. The execution tests deliberately run
    from the ``.py`` -- it is the source of truth, and reading it keeps the
    suite independent of committed outputs -- so a ``.ipynb`` that has fallen
    behind its script executes nothing, breaks nothing, and is silently wrong
    for every reader who opens the notebook instead of the script.

    That matters more than it used to now that the notebooks are linked from
    the README and rendered by GitHub: what a student reads is the ``.ipynb``.

    Compares cell *sources* only. Outputs, execution counters and metadata
    differ freely -- re-running two cells in JupyterLab is not drift.
    """

    #: Every pair, not just the numbered teaching set: the photo* notebooks
    #: are read by students too, and Scratchbook is paired like the rest.
    @staticmethod
    def pairs():
        return sorted(
            (path, path.with_suffix(".ipynb"))
            for path in NOTEBOOK_DIR.glob("*.py")
            if path.with_suffix(".ipynb").exists()
        )

    def test_there_are_pairs_to_check(self):
        self.assertGreaterEqual(len(self.pairs()), 10)

    def test_no_notebook_has_drifted_from_its_script(self):
        import jupytext

        for script, notebook in self.pairs():
            with self.subTest(notebook=script.name):
                from_script = [c.source.strip() for c in jupytext.read(script).cells]
                from_notebook = [c.source.strip() for c in jupytext.read(notebook).cells]
                self.assertEqual(
                    from_script,
                    from_notebook,
                    f"{script.name} and {notebook.name} differ in cell source -- "
                    "edit the .py, then run "
                    f"`jupytext --to ipynb notebooks/{script.name} "
                    f"-o notebooks/{notebook.name}` and re-execute it",
                )

    def test_every_notebook_has_a_script_beside_it(self):
        """The other direction: a .ipynb with no .py is unpaired and will not
        be executed by this suite at all."""
        orphans = [
            path.name
            for path in NOTEBOOK_DIR.glob("*.ipynb")
            if not path.with_suffix(".py").exists()
        ]
        self.assertEqual(orphans, [], f"unpaired notebooks: {orphans}")


@unittest.skipIf(SKIP, "MT_SKIP_NOTEBOOKS=1")
class TestStallWatchdog(unittest.TestCase):
    """The watchdog must kill an idle kernel and leave a busy one alone.

    Both run a throwaway one-cell notebook with a stall window of a few
    seconds, so each costs a kernel start and little more.
    """

    STALL_S = 3.0
    POLL_S = 0.5

    @staticmethod
    def one_cell(source):
        import nbformat

        return nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(source)])

    def test_a_sleeping_kernel_is_killed_and_reported_as_stalled(self):
        from nbclient.exceptions import DeadKernelError

        started = time.monotonic()
        with self.assertRaisesRegex(DeadKernelError, "Kernel stalled"):
            run_notebook(
                self.one_cell("import time\ntime.sleep(600)"),
                stall_s=self.STALL_S,
                poll_s=self.POLL_S,
            )
        # Killed in about stall + kernel start-up, far inside the cell timeout.
        self.assertLess(time.monotonic() - started, 60)

    def test_a_busy_kernel_outlasting_the_window_is_left_alone(self):
        source = (
            "import time\n"
            f"end = time.monotonic() + {3 * self.STALL_S}\n"
            "while time.monotonic() < end:\n"
            "    pass\n"
            "print('finished')"
        )
        notebook = run_notebook(self.one_cell(source), stall_s=self.STALL_S, poll_s=self.POLL_S)
        self.assertEqual(notebook.cells[0].outputs[0]["text"].strip(), "finished")


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
