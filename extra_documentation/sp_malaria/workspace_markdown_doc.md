## LRMA Special Project: Malaria

This is the Long Read Methods and Applications (LRMA) Special Project Workspace for short-read malaria data processing.  This workspace is a rapidly-evolving prototype and is not yet ready for heavy usage.  The current focus of this workspace is _P. falciparum_, but the processing steps here are generalized and can be adapted to other _Plasmodium_ species.

## Variant Calling Pipeline

As part of this workspace there are workflows to call variants on both single samples, and for joint calling across cohorts of samples.  The 

The main variant calling pipeline has has the following high-level structure:

![LRMA SP Malaria Variant Calling](https://github.com/broadinstitute/long-read-pipelines/raw/jts_kvg_sp_malaria/extra_documentation/sp_malaria/lrma_sr_malaria_pipeline_diagram.png)

## Data

### Datasets

This workspace has been seeded with the [PF7](https://www.malariagen.net/apps/pf7/) samples, and the [crosses](https://www.malariagen.net/parasite/p-falciparum-genetic-crosses).  Over time further samples will be added.

### Data Structure

The data processing is broken down into three levels (similar to other LRMA projects) in the following Terra data tables:
*  Sample (flowcell data)
*  Sample Set (sample data / single-sample calling)
*  Sample Set Set (cohort data for joint calling)

_Sample / Flowcell_ data consists of reads from a single flowcell.  The sample from which these reads have been processed may or may not be represented in other flowcells.

_Sample Set_ data consists of all data from a specific sample.  This may include data from multiple flowcells.

_Sample Set Set / Cohort_ data consists of data from multiple samples.  


## Task List
The following is a preliminary task list for this workspace.  It may or may not be up-to-date:

| Level | Category                     | To-do items                                         | Assignment | Status |
|-------|------------------------------|-----------------------------------------------------|------------|--------|
| 0     | Download and prepare data    | Download Pf7                                        | Kiran      | Done    |
|       |                              | Experimental crosses                                | Jonn      | Done   |
|       |                              | Additional Broad datasets                           | Wes/Jonn      |        |
| 1     | Sample-level processing      | Alignment with bwa-mem2 to Pf3D7                    | Jonn       | Done       |
|       |                              | Mark duplicates                                     | Jonn       | Done       |
|       |                              | BQSR                                                | Jonn       | Done       |
| 2     | Participant-level processing | Sample merging                                      | Kiran      | Done   |
|       |                              | Call with DeepVariant                               | Kiran      | Done   |
|       |                              | Call with HaplotypeCaller                           | Jonn       | Done       |
|       |                              | Apply VQSR                               | Jonn      | Done       |
|       |                              | Annotate with SnpEff                                | Kiran      | Done       |
| 3     | Cohort-level processing      | Apply VQSR                               | Jonn      | Done       |
|       |                              | Joint call DeepVariant calls with GLNexus           | Kiran      | Done   |
|       |                              | Joint call HaplotypeCaller calls with GenotypeGVCFs | Jonn       | Done       |
|       |                              | Joint call HaplotypeCaller calls with GVS | Jonn       |        |
|       |                              | Estimate MOI with DEploid                           | Kiran      |        |
|       |                              | **Call CNVs with gCNV**                                 | Jonn       |        |
|       |                              | Convert SNV/indel calls to Hail MatrixTables        | Kiran      | Done   |
|       |                              | Convert SNV/indel calls to Zarr Store        | Kiran      | Done   |
| 4     | Evaluation                   | Comparison to long read references                  | Jonn       |        |
|       |                              | Compute Mendelian violation rate                    | Jonn       |        |
|       |                              | Compare results to Pf3k/Pf7 manuscripts             | Jonn       |        |
| 5     | Biology                      | Drug resistance markers                              | Kiran      | Done   |
|       |                              | HRP2/HRP3 deletions                                 | Jonn       |        |
|       |                              | Population genetics                                 | Jonn       |        |
|       |                              | **Mapping (lat/long)**                                             | Jonn / Bridget      | In Progress       |
