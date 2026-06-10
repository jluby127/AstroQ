"""Pytest fixtures shared by the tests.

The tests in :mod:`tests.test_sample` reference example inputs (configs,
request/allocation/past CSVs) by their literal in-repo paths under
``examples/hello_world/`` and ``examples/bench/``. AstroQ writes its outputs
into ``[global] workdir`` from each config, which points back into
``examples/...``. Running the suite in-place would dirty the repo.

The autouse session fixture below sidesteps that:

1. Once per session, copy ``examples/`` into a fresh ``tmp_path`` directory.
2. ``chdir`` into that directory for the duration of the session.
3. Restore the original CWD at teardown; the tmp directory is cleaned up by
   pytest automatically.

Effect: every test reads its fixtures from the in-repo ``examples/`` (via the
copy) and writes outputs into the tmp tree. The repo working directory stays
clean. Tests that depend on artifacts created by earlier tests (e.g.
``test09_hdf5_validation`` reads the ``semester_planner.h5`` produced by
``test01_helloworld``) still work because the tmp directory persists for the
session.
"""

import os
import shutil
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session", autouse=True)
def sandbox_examples(tmp_path_factory):
    """Copy ``examples/`` to a tmp dir and ``chdir`` there for the whole session.

    Yields the sandbox path in case a test wants it explicitly.
    """
    sandbox = tmp_path_factory.mktemp("astroq_examples_sandbox")
    shutil.copytree(REPO_ROOT / "examples", sandbox / "examples")

    original_cwd = os.getcwd()
    os.chdir(sandbox)
    print(f"\n[conftest] test sandbox: {sandbox}")
    try:
        yield sandbox
    finally:
        os.chdir(original_cwd)
