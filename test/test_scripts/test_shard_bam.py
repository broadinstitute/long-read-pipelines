def test_shard_bam(script_runner):
    ret = script_runner.run("../../docker/lr-pb/shard_bam.py", "-p", "shard", "-n", "2", "test/data/tiny.bam")

    exp = """Reading index (test/data/tiny.bam.pbi). This may take a few minutes...
Writing 2 shards using 2 threads...
  - wrote 18 reads to shard0.bam
  - wrote 27 reads to shard1.bam
Sharded 45/45 reads across 2 shards.
"""

    assert ret.success
    assert ret.stdout == exp
