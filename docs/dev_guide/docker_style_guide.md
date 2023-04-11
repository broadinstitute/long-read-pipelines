# Docker Style Guide
This document will provide a guide for writing Dockerfiles for the pipelines.  
The guide will provide a list of best practices for writing Dockerfiles, and will also 
provide a list of common mistakes to avoid. The guide is for those wanting to contribute
Dockerfiles to the pipelines repository.


## Placement of Dockerfiles
Dockerfiles should be placed in the `docker` directory of the repository. The `docker` directory
contains a subdirectory for each Docker related resource. The subdirectory name should 
contian an abreviation data type the docker image is for followed by the name of the tool.
For example, the docker image for the `bwa` aligner would be placed in the `docker/lr-bwa` directory, 

## Dockerfile Guidelines: 
_This section should outline the guidelines for creating Dockerfiles, including the format, structure, and best practices for creating efficient and maintainable Docker images._



## Image Naming and Tagging Guidelines: 
_This section should cover guidelines for naming and tagging Docker images, including how to name images based on their functionality, versioning, and tagging best practices.

## Docker Compose Guidelines: 
_This section should provide guidelines for creating Docker Compose files, including best practices for structuring and organizing services, environment variables, and volumes.

## Repository Organization: 
_This section should outline the guidelines for organizing Docker-related scripts and files within the repository, including the directory structure, naming conventions, and how to maintain version history.

## Security Guidelines: 
_This section should cover security best practices for Docker images and containers, including how to minimize vulnerabilities, handle secrets and sensitive data, and how to use Docker security features such as namespaces and SELinux.

## Contribution Guidelines: 
_This section should outline the guidelines for contributing Docker-related scripts and files to the repository, including how to submit pull requests, how to review and test changes, and how to handle conflicts and merge requests.

## Testing and CI/CD Guidelines: 
_This section should provide guidelines for testing Docker images and containers, including best practices for creating automated tests and integrating with CI/CD pipelines to ensure consistent builds and deployments.
