# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is the Broad Institute's long-read genomics pipelines repository, containing WDL (Workflow Description Language) workflows for processing PacBio and Oxford Nanopore sequencing data. The pipelines are designed for Google Cloud Platform using Cromwell as the workflow engine.

## Commands

### Testing
- Run all tests: `python scripts/ci/run_test_suite.py`
- Run specific test module: `pytest test/test_scripts/test_shard_bam.py`
- Run WDL validation: `scripts/wdl/validate.wdls.sh`
- Run tests with tox: `tox`

### Building
- Build all Docker images: `cd docker && make`
- Build specific Docker image: `cd docker/lr-<tool> && make`
- Create WDL dependencies: `cd wdl && zip -r lr_wdls.zip *`

### Development Workflow
- Validate WDL syntax before committing: `womtool validate <file.wdl>`
- Test Docker builds locally before pushing changes
- Use CI test suite to validate pipeline functionality

## Architecture

### Directory Structure
- `wdl/` - Core WDL workflows and tasks
  - `pipelines/` - High-level workflows organized by platform (ONT, PacBio, TechAgnostic)
  - `tasks/` - Reusable task definitions organized by analysis type
  - `structs/` - Common data structures (RuntimeAttr, DataTypeParameters)
- `docker/` - Dockerfiles for pipeline tools, organized as `lr-<tool>/` or `sr-<tool>/`
- `test/` - Test configurations and test data
- `scripts/` - Utility scripts for CI/CD, monitoring, and workflow management

### WDL Architecture
- **Platform-specific pipelines**: ONT/, PacBio/, TechAgnostic/ organize workflows by sequencing technology
- **Analysis-type organization**: Assembly/, VariantCalling/, Preprocessing/, etc. within each platform
- **Task modularity**: Tasks are separated from workflows for reusability and maintainability
- **Runtime configuration**: Uses RuntimeAttr struct from structs/Structs.wdl for consistent resource management

### Docker Organization
- Each tool has its own Docker subdirectory with Dockerfile and Makefile
- Naming convention: `lr-<tool>` for long-read tools, `sr-<tool>` for short-read tools
- Docker images are built and pushed to Google Container Registry at `us.gcr.io/broad-dsp-lrma/`

## Key Patterns

### WDL Task Structure
Tasks follow this section order: meta, parameter_meta, input, command, output, runtime. Runtime sections use the RuntimeAttr struct pattern for consistent resource allocation.

### Docker Build Process
Each Docker directory contains a Makefile that handles building and pushing images. The main docker/Makefile coordinates building multiple images.

### Testing Framework
The CI system uses a comprehensive test suite that validates WDL syntax, runs integration tests against real data, and compares outputs against known-good results stored in Google Cloud Storage.

### Platform-Agnostic Design
TechAgnostic pipelines provide common functionality that works across sequencing platforms, while platform-specific directories contain technology-specific optimizations.