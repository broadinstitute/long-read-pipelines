[![Generic badge](https://img.shields.io/badge/version-2.1.20-blue.svg)](https://shields.io/)
![CI](https://github.com/broadinstitute/long-read-pipelines/workflows/CI/badge.svg?branch=master&event=push)

# Long read pipelines
This repository contains pipelines for processing of long read data from PacBio and/or Oxford Nanopore platforms.  The pipelines are written in [WDL 1.0](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#introduction) intended for use with Google Cloud Platform via the scientific workflow engine, [Cromwell](https://github.com/broadinstitute/cromwell).  Processing is designed to be reasonably consistent between both long read platforms, and use platform-specific options or tasks where necessary.

High level workflows can be found in the wdl/ directory.
