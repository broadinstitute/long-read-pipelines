# Long read pipelines
This repository contains pipelines for processing of long read data from PacBio and/or Oxford Nanopore platforms.  The pipelines are written in [WDL 1.0](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#introduction) intended for use with Google Cloud Platform via the scientific workflow engine, [Cromwell](https://github.com/broadinstitute/cromwell).  Processing is designed to be reasonably consistent between both long read platforms, and use platform-specific options or tasks where necessary.

## Available pipelines:
The following pipelines are implemented:

| Pipeline                        | Description                                                                                              |
|---------------------------------|----------------------------------------------------------------------------------------------------------|
| LRWholeGenomeSingleSample.wdl   | Error correction, alignment, and variant discovery on >= 1 {SMRT,flow}cell of data from the same sample. |
| LRTranscriptomeSingleSample.wdl | Error correction and splice-aware alignment on >= 1 {SMRT,flow}cell of data from the same sample.        |
