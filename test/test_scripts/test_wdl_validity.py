import errno
import glob
import os
import subprocess
import pytest
from pathlib import Path


@pytest.mark.parametrize(
    "tool, subcommand, dir_names",
    (
        ["womtool", "validate", ["pipelines", "tasks", "deprecated_wdl"]],
        ["miniwdl",  "check", ["pipelines", "tasks"]],  # miniwdl validation is stricter and takes more time than womtool, so we only run it on 'pipelines' and 'tasks' since 'deprecated_wdl' will eventually be deleted.
    )
)
def test_wdl_validity(script_runner, tool: str, subcommand: str, dir_names: list):

    for dir_name in dir_names:
        retsum = 0
        wdl_dir = Path(__file__).resolve().parents[2] / dir_name
        if not wdl_dir.exists():
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), wdl_dir)
        wdl_paths = wdl_dir / "**/*.wdl"

        for wdl in glob.glob(str(wdl_paths), recursive=True):
            print(f'{wdl}:')
            ret = subprocess.run([tool, subcommand, wdl])
            retsum += ret.returncode

        assert retsum == 0
