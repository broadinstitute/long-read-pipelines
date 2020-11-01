import pysam
import pathlib
import tempfile
import shutil


def get_zmws(file):
    bf = pysam.Samfile(file, 'rb', check_sq=False)

    zmws = set()
    for read in bf:
        zmw = read.get_tag("zm")
        zmws.add(zmw)

    bf.close()

    return zmws


def test_extract_uncorrected_reads(script_runner):
    subreads_bam = "test/test_data/for_scripts/subreads.bam"
    consensus_bam = "test/test_data/for_scripts/consensus.bam"

    testdir = tempfile.mkdtemp()
    out_bam = f"{testdir}/out.bam"

    pathlib.Path(testdir).mkdir(parents=True, exist_ok=True)
    ret = script_runner.run("docker/lr-pb/extract_uncorrected_reads.py", "-o", out_bam, subreads_bam, consensus_bam)

    subreads_zmws = get_zmws(subreads_bam)
    consensus_zmws = get_zmws(consensus_bam)

    exp_zmws = subreads_zmws.difference(consensus_zmws)
    act_zmws = get_zmws(out_bam)

    print(exp_zmws)
    print(act_zmws)

    shutil.rmtree(testdir)

    assert ret.success
    assert exp_zmws == act_zmws

