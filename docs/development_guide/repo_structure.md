# Repository Structure

The repository file and directory is as follows:

* LICENSE: The license for the repository.
* README.md: This document, which provides an overview of the repository.
* VERSION: The version number of the repository.
* cloudbuild.yaml: A Cloud Build configuration file that defines how the repository is built.
* docker: A directory containing Dockerfiles for building the pipelines.
* docs: A directory containing documentation for the pipelines.
* requirements.txt: A file listing the Python dependencies for the pipelines.
* resources: A directory containing resources used by the pipelines.
* scripts: A directory containing scripts used to build and run the pipelines.
* site: A directory containing static files for the website.
* test: A directory containing unit tests for the pipelines.
* tox.ini: A configuration file for the tox test runner.
* wdl: A directory containing WDL files for the pipelines.

## WDL Directory Structure:
The WDL directory is further divided into subdirectories. The subdirectories are as follows:
* tasks: A directory containing WDL tasks for the pipelines.
* pipelines: A directory containing WDL workflows for the pipelines.
* structs: A directory containing WDL structs for the pipelines.

### Tasks Directory Structure:
The task directory has an additional subdirectory to organize wdl tasks by analysis type. The subdirectories are as follows:

* alignment: A directory containing WDL tasks for aligning reads to a reference genome.
* annotation: A directory containing WDL tasks for annotating reads.
* assembly: A directory containing WDL tasks for assembling reads.
* epigenomics: A directory containing WDL tasks for epigenomics analysis.
* preprocessing: A directory containing WDL tasks for preprocessing reads.
* qc: A directory containing WDL tasks for quality control of reads.
* transcriptomics: A directory containing WDL tasks for transcriptomics analysis.
* utility: A directory containing WDL tasks for utility functions.
* variantcalling: A directory containing WDL tasks for calling variants from reads.
* visualization: A directory containing WDL tasks for visualizing data.

### Pipelines Directory Structure:
The pipelines directory has two additional subdirectories to organize wdl workflows, first by platform then by analysis type. 
The first level subdirectories are as follows:
* Illumina: A directory containing WDL workflows for processing Illumina data.
* PacBio: A directory containing WDL workflows for processing PacBio data.
* ONT: A directory containing WDL workflows for processing Oxford Nanopore data.
* TechAgnostic: A directory containing WDL workflows for processing data from PacBio and ONT platforms.

The second level subdirectories are as follows:
* alignment: A directory containing WDL workflows for aligning reads to a reference genome.
* annotation: A directory containing WDL workflows for annotating reads.
* assembly: A directory containing WDL workflows for assembling reads.
* epigenomics: A directory containing WDL workflows for epigenomics analysis.
* multianalysis: A directory containing WDL workflows for performing multiple analyses on reads.
* preprocessing: A directory containing WDL workflows for preprocessing reads.
* utility: A directory containing WDL workflows for utility functions.
* variantcalling: A directory containing WDL workflows for calling variants from reads.
