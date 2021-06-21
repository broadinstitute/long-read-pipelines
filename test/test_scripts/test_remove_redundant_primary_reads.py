import pysam
import pathlib
import tempfile
import shutil


def get_primary_read_name_counts(file):
    bf = pysam.Samfile(file, 'rb', check_sq=False)

    counts = {}
    for read in bf:
        if not read.is_secondary and not read.is_supplementary:
            counts[read.query_name] = counts.get(read.query_name, 0) + 1

    bf.close()

    return counts


def test_remove_redundant_primary_reads(script_runner):
    bam = "test/test_data/for_scripts/redundant_reads.bam"
    testdir = tempfile.mkdtemp()

    prefix = f"{testdir}/filtered"

    pathlib.Path(testdir).mkdir(parents=True, exist_ok=True)
    ret = script_runner.run("docker/lr-utils/remove_redundant_primary_reads.py", "-p", prefix, bam)

    assert ret.success

    # Silence message about the .bai file not being found.
    pysam.set_verbosity(0)

    counts_orig = get_primary_read_name_counts(bam)
    counts_filt = get_primary_read_name_counts(f'{prefix}.bam')

    all_keys = set(list(counts_orig.keys()) + list(counts_filt.keys()))

    for read_name in all_keys:
        assert read_name in counts_orig
        assert read_name in counts_filt
        assert counts_filt[read_name] == 1

    shutil.rmtree(testdir)
