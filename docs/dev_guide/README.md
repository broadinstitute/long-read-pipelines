# Development Guide

## Overview:
This repository contains pipelines for processing of long read data from PacBio and/or 
Oxford Nanopore platforms. The pipelines are written in WDL 1.0 intended for use with 
Google Cloud Platform via the scientific workflow engine, Cromwell. Processing is 
designed to be reasonably consistent between both long read platforms, and use 
platform-specific options or tasks where necessary.

To use the pipelines, you will need access to Cromwell server. The pipelines are 
designed to be run on Google Cloud Platform, but can be run on other platforms with some modifications. 


## Repository Structure:

The repository includes files such as the LICENSE and README.md files, which provide legal and informational 
overviews of the repository. Other directories such as docker, docs, resources, scripts, 
site, test, and wdl contain various files and directories that are important for 
building and testing the software in the repository. 
See [Repository Structure](./repo_structure.md) documentation for further details and [Contributing Guidelines](#contributing-guidelines) for 
information on how to contribute to the repository.

## Workflow Scripts:
All workflow scripts are located in the wdl directory and are written in WDL 1.0.


## Docker Containers:

The WDL workflows in this repository are designed to be run using Docker containers. 
This provides a number of advantages, including the ability to run the workflows on a 
variety of platforms and the ability to isolate the workflows from each other.
Many of these workflows use specialized containers that are built from the Dockerfiles in the 
docker directory. The docker directory contains Dockerfiles and other scripts for 
building containers for several tools. The docker containers are pushed to 
Google Container Registry (GCR), for internal use. External audiences interested in running workflows
using these containers should build and if needed push them to repository they have access to.

## Testing:

Scripts for running tests are located in the test directory. The test directory contains
scripts for running tests using Tox. Tox is a Python-based test automation tool that 
can be used to run tests in a variety of environments. The scripts are run through
GitHub Actions, which are configured in the .github/workflows directory. The GitHub Actions
are triggered by pushes to the repository and pull requests. 

## Contributing Guidelines:

Please adhere to the following best practices if contributing to this repository:

1. Read Style Guide: Before making any changes to the code, it's important to read the style guide. The style guide contains information on how to write code that is consistent with the rest of the codebase. See [WDL Style Guide](./wdl_style_guide.md) and [Docker Style Guide](./docker_style_guide.md) for more information.
2. Create a new branch: When making contributions to a repository, it's important to create a new branch for each change you make. The name of the branch should begin with your initials followed by an underscore and a short description of the change. For example, if Jack Sully is making a change to the README file, the name might be `js_update_readme`.
3. Keep commits small and focused: When making changes to the code, it's important to keep your commits small and focused on a specific task. This makes it easier for others to review your changes and also makes it easier to roll back changes if necessary.
4. Write clear commit messages: When committing changes to the repository, it's important to write clear and concise commit messages that describe what changes were made. This helps others understand the changes you made and why you made them.
5. Test your changes: Before submitting your changes, make sure to test them thoroughly to ensure they work as intended. This helps reduce the chance of introducing bugs or issues into the codebase.
6. Submit a pull request: Once you have made your changes and tested them, submit a pull request to the main repository. Make sure to include a clear description of the changes you made and why you made them. This makes it easier for others to review and merge your changes into the main repository.
7. Add reviewers: Once you have submitted your pull request, add reviewers to the pull request. This will notify them that you have submitted a pull request and they should review it. It's important to add at least one reviewer to your pull request.
8. Merging pull requests: Once your pull request has been reviewed and approved, it can be merged into the main repository. It's important to merge pull requests using the "Squash and merge" option. This will squash all commits in the pull request into a single commit, which makes it easier to track changes in the repository.
