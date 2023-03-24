import errno
import glob
import os
import subprocess
import pytest
from pathlib import Path


@pytest.mark.parametrize(
    "dir_name",
    ["pipelines", "tasks", "deprecated_wdl"],
)
def test_wdl_validity(script_runner, dir_name):
    retsum = 0

    wdl_dir = Path(__file__).resolve().parents[2] / dir_name
    if not wdl_dir.exists():
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), wdl_dir)
    wdl_paths = wdl_dir / "**/*.wdl"

    for wdl in glob.glob(str(wdl_paths), recursive=True):
        print(f'{wdl}:')
        ret = subprocess.run(["womtool", "validate", wdl])
        retsum += ret.returncode

    assert retsum == 0
