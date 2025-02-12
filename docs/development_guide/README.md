# Development Guide Overview

This development guide provides information on the structure of the repository, testing infrastructure,
style guides, and contributing guidelines. The hope is that this guide will help developers 
create and maintain code that is consistent with the rest of the repository.

## Repository Structure

The repository includes files such as the LICENSE and README.md files, which provide legal and informational 
overviews of the repository. Other directories such as docker, docs, resources, scripts, 
site, test, and wdl contain various files and directories that are important for 
building and testing the software in the repository.  
See [Repository Structure](./repo_structure.md) documentation for further 
details and [Contributing Guidelines](#contributing-guidelines) for 
information on how to contribute to the repository.

## Workflow Scripts

All workflow scripts are located in the `/wdl` directory and are written in WDL 1.0 
and intended for use with Google Cloud Platform via the scientific workflow engine, Cromwell.
The WDL scripts are divided into three subdirectories: `tasks`, `structs`, and `pipelines`;
then further divided by sequencing platform and analysis type.
See [Repository Structure](./wdl_style_guide.md) for more information on directory structure.

## Docker Containers

The WDL workflows in this repository are designed to be run using Docker containers. 
This provides a number of advantages, including the ability to run the workflows on a 
variety of platforms and the ability to isolate the workflows from each other.
Many of these workflows use specialized containers that are built from the Dockerfiles in the 
docker directory. The docker directory contains Dockerfiles and other scripts for 
building containers for several tools. The docker containers are pushed to 
Google Container Registry (GCR) called `us.gcr.io/broad-dsp-lrma`, for internal use. External audiences interested in running workflows
using these containers should build and if needed push them to repository they have access to.

## Testing

Scripts for running tests are located in the test directory. The test directory contains
scripts for running tests using Tox. Tox is a Python-based test automation tool that 
can be used to run tests in a variety of environments. The scripts are run through
GitHub Actions, which are configured in the .github/workflows directory. The GitHub Actions
are triggered by pushes to the repository and pull requests. 

## Workflow Deployment

The workflows in this repository are deployed to Terra using Dockstore. Dockstore is a
platform for sharing Docker containers and workflows. The workflows are registered using the
`dockstore.yml` file in the root directory of the repository. The `dockstore.yml` file
contains information about the workflows, including the location of the WDL files and if
available the location of example input JSON files. The workflows published in Dockstore 
are automatically updated when changes are made to the repository. If you would like to
add a new workflow to this repository and have it published in Dockstore, please
update the `dockstore.yml` file in your feature branch. This should be enough for 
Dockstore github app to automatically add your workflow branch version to
the Dockstore repository.


## Contributing Guidelines

Please adhere to the following best practices if contributing to this repository:

1. **Read Style Guide**: Before making any changes to the code, it's important to read the style guide. The style guide contains information on how to write code that is consistent with the rest of the codebase. See [WDL Style Guide](./wdl_style_guide.md) and [Docker Style Guide](./docker_style_guide.md) for more information.
2. **Create a new branch**: When making contributions to a repository, it's important to create a new branch for each change you make. The name of the branch should begin with your initials followed by an underscore and a short description of the change. For example, if Janet Sully is making a change to the README file, the name might be `js_update_readme`.
3. **Keep commits small and focused**: When making changes to the code, it's important to keep your commits small and focused on a specific task. This makes it easier for others to review your changes and also makes it easier to roll back changes if necessary.
4. **Write clear commit messages**: When committing changes to the repository, it's important to write clear and concise commit messages that describe what changes were made. This helps others understand the changes you made and why you made them.
5. **Test your changes**: Before submitting your changes, make sure to test them thoroughly to ensure they work as intended. This helps reduce the chance of introducing bugs or issues into the codebase.
6. **Submit a pull request**: Once you have made your changes and tested them, submit a pull request to the main repository. Make sure to include a clear description of the changes you made and why you made them. This makes it easier for others to review and merge your changes into the main repository.
7. **Add reviewers**: Once you have submitted your pull request, add reviewers to the pull request. This will notify them that you have submitted a pull request and they should review it. It's important to add at least one reviewer to your pull request.
8. **Merging pull requests**: Once your pull request has been reviewed and approved, it can be merged into the main repository. It's important to merge pull requests using the "Squash and merge" option. This will squash all commits in the pull request into a single commit, which makes it easier to track changes in the repository.

### External Contributions

External contributions are welcome. Please follow the guidelines and submit a pull request.

#### Review Criteria

To maintain quality and consistency, contributions are reviewed based on:

- **Adherence to coding standards**: The code should follow the project's coding standards, such as naming conventions, indentation, and code style.
- **Passing tests**: All existing tests should pass after the changes are applied. New tests should be added to cover the new functionality.
- **Documentation updates**: If the contribution introduces new features or changes existing workflows, the documentation should be updated accordingly.
- **Overall quality**: The contribution should be of overall high quality, with well-structured code, clear comments, and a minimal impact on other parts of the codebase.

#### Review Process

The review process is designed to ensure that contributions are of high quality and consistent with the rest of the codebase. The review process consists of the following steps:
1. **Code Review**: Reviewers provide feedback for improvements or issue identification.
2. **Testing**: Ensure code functions as intended without introducing bugs.
3. **Merging**: Approved code becomes part of the main repository.
4. **Release**: Automatic generation of releases upon merging to the main branch.


