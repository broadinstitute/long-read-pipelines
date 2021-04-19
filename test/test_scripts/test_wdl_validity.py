import glob
import subprocess


def test_wdl_validity(script_runner):
    retsum = 0
    for wdl in glob.glob("wdl/**.wdl", recursive=True):
        print(f'{wdl}:')
        ret = subprocess.run(["womtool", "validate", wdl])
        retsum += ret.returncode

    assert retsum == 0
