import pysam
import pathlib


def get_read_zmw_counts(file):
    bf = pysam.Samfile(file, 'rb', check_sq=False)

    zmw_counts = {}
    for read in bf:
        zmw = read.get_tag("zm")
        zmw_counts[zmw] = zmw_counts.get(zmw, 0) + 1

    bf.close()

    return zmw_counts


def test_shard_bam(script_runner):
    bam = "test/test_data/for_scripts/shard_bam_test_file.bam"
    testdir = "test/test_output"
    prefix = f"{testdir}/shard"
    num_shards = 2

    pathlib.Path(testdir).mkdir(parents=True, exist_ok=True)
    ret = script_runner.run("../../docker/lr-pb/shard_bam.py", "-p", prefix, "-n", str(num_shards), bam)

    assert ret.success

    files_with_zmws = {}
    zmw_counts_orig = get_read_zmw_counts(bam)
    for i in range(0, num_shards):
        zmw_counts_shard = get_read_zmw_counts(f'{prefix}{i}.bam')

        for zmw in zmw_counts_shard:
            zmw_counts_orig[zmw] = zmw_counts_orig.get(zmw, 0) - zmw_counts_shard.get(zmw, 0)

            if zmw not in files_with_zmws:
                files_with_zmws[zmw] = set()

            files_with_zmws[zmw].add(i)

    for zmw in files_with_zmws:
        # Verify that ZMWs only ever appear in one shard.
        assert len(files_with_zmws[zmw]) == 1

        # Verify that every instance of a ZMW seen in the original file are accounted for in the shards.
        assert zmw_counts_orig[zmw] == 0
