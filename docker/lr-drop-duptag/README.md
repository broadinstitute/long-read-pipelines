# lr-drop-duptag

`lr_drop_duptag` streams a BAM and removes duplicate lowercase `mx` auxiliary
tags from each alignment record.

```bash
lr_drop_duptag \
  --input in.bam \
  --output-bam out.bam \
  --mismatches mismatches.txt \
  --threads 8 \
  --progress-seconds 60
```

For records with multiple `mx` tags, the first `mx` value seen in the record is
kept and all later `mx` tags are dropped. The read name is written to the
mismatch report when any later `mx` value differs from the first, or when a
record has three or more `mx` tags.

Progress is written to stderr as periodic log lines. Set `--progress-seconds 0`
to disable periodic progress updates.
