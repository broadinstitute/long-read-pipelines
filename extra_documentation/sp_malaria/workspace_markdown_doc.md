# _P. falciparum_ Short Read Whole Genome Workspace
This is the workspace for short read whole genome variant discovery and analysis in _Plasmodium falciparum_.   This workspace can call variants in a single-sample, joint call cohorts of samples, and perform various tertiary analyses (e.g. drug resistance screening, rapid diagnostic test evasion screening, etc.).

While the current focus of this workspace is _P. falciparum_, but the processing steps here are generalized and can be adapted to other _Plasmodium_ species.

## Variant Calling Pipeline

As part of this workspace there are workflows to call variants on both single samples, and for joint calling across cohorts of samples.

The main variant calling pipeline has has the following high-level structure:

![LRMA SP Malaria Variant Calling](https://github.com/broadinstitute/long-read-pipelines/raw/jts_kvg_sp_malaria/extra_documentation/sp_malaria/lrma_sr_malaria_pipeline_diagram_high_level.png)

## Data

### Datasets

The following datasets are currently in this workspace:
- [PF7](https://www.malariagen.net/apps/pf7/)
- The MalariaGEN [crosses](https://www.malariagen.net/parasite/p-falciparum-genetic-crosses)
- 2022 data collected in Senegal
- 2019 data collected in Senegal

### Data Structure

The data processing is broken down into three levels (similar to other LRMA projects) in the following Terra data tables:
*  Sample (flowcell data)
*  Sample Set (sample data / single-sample calling)
*  Sample Set Set (cohort data for joint calling)

_Sample / Flowcell_ data consists of reads from a single flowcell.  The sample from which these reads have been processed may or may not be represented in other flowcells.

_Sample Set_ data consists of all data from a specific sample.  This may include data from multiple flowcells that belong to the same "participant" (i.e. same strain / clone).

_Sample Set Set / Cohort_ data consists of data from multiple samples.  

