# Repository Structure

```angular2html
├── cloudbuild.yaml
├── LICENSE
├── README.md
├── mkdocs.yml
├── requirements.txt
├── tox.ini
├── VERSION
├── docker
│   └── ...
├── docs
│   └── ...
├── resources
│   └── ...
├── scripts
│   └── ...
├── terra
│   └── ...
├── test
│   └── ...
└── wdl
    └── ...
```


The repository file and directory is as follows:

* LICENSE: The license for the repository.
* README.md: This document, which provides an overview of the repository.
* VERSION: The version number of the repository.
* cloudbuild.yaml: A Cloud Build configuration file that defines how the repository is built.
* docker: Contains Dockerfiles for building docker images used by pipelines.
* docs: Contains documentation for the pipelines and a developer's guide.
* requirements.txt: A file listing the Python dependencies for the pipelines.
* resources: A directory containing resources used by the pipelines.
* mkdocs.yml: A configuration file for the mkdocs documentation generator.
* scripts: Contains scripts used by the repository (e.g. webpage creation).
* test: Contains tests for the pipelines.
* tox.ini: A configuration file for the tox test runner.
* wdl: Contains WDL files.

## WDL Directory Structure

```angular2html
└── wdl
    └── pipelines
    │  └── ...
    └── structs
    │  └── ...
    └── tasks
         └── ...
```
The WDL directory is further divided into subdirectories. The subdirectories are as follows:

* tasks: Contains WDL files with a list of tasks to be imported and used by pipeline WDLs.
* pipelines: Contains WDL files with workflow blocks.
* structs: Contains WDL structs for the pipelines.

### Tasks Directory Structure  
The task directory has an additional subdirectory to organize wdl tasks by analysis type. The subdirectories are as follows:  

```angular2html
└── wdl
    └── tasks
    │  └── alignment
    │  │  └── ...
    │  └── annotation
    │  │  └── ...
    │  └── assembly
    │  │  └── ...
    │  └── epigenomics
    │  │  └── ...
    │  └── preprocessing
    │  │  └── ...
    │  └── qc
    │  │  └── ...
    │  └── transcriptomics
    │  │  └── ...
    │  └── utility
    │  │  └── ...
    │  └── variantcalling
    │  │  └── ...
    │  └── visualization
    │  │  └── ...
```

### Pipelines Directory Structure
The pipelines directory has two additional subdirectories to organize wdl workflows, first by platform then by analysis type.  

The first level subdirectories are as follows:  

```angular2html
└── wdl
    └── pipelines
    │  └── Illumina
    │  │  └── ...
    │  └── PacBio
    │  │  └── ...
    │  └── ONT
    │  │  └── ...
    │  └── TechAgnostic
    │  │  └── ...
```

The second level subdirectories are as follows:  

```angular2html
└── wdl
    └── pipelines
    │  └── Illumina
    │  │  └── alignment
    │  │  │  └── ...
    │  │  └── annotation
    │  │  │  └── ...
    │  │  └── assembly
    │  │  │  └── ...
    │  │  └── epigenomics
    │  │  │  └── ...
    │  │  └── multianalysis
    │  │  │  └── ...
    │  │  └── preprocessing
    │  │  │  └── ...
    │  │  └── utility
    │  │  │  └── ...
    │  │  └── variantcalling
    │  │  │  └── ...
```
