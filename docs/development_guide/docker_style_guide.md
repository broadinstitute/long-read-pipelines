# Docker Style Guide
This document will provide a guide for writing Dockerfiles for the pipelines.  
The guide will provide a list of best practices for writing Dockerfiles, and will also 
provide a list of common mistakes to avoid. The guide is for those wanting to contribute
Dockerfiles to the pipelines repository.


## Dockerfiles Organization
All docker related resources should be placed in the `docker` directory of the repository. The `docker` directory
contains a subdirectory for each Dockerfile and its related resources. The subdirectory name should 
start with an abbreviation of data type the docker tool will process followed by the name of the tool.
For example, the docker image for the `bwa` aligner that will process long reads would be placed in the `docker/lr-bwa` directory. 

## Docker Subdirectory Folder
Each Docker subdirectory should contain the following files and folders:  

- `Dockerfile`: The Dockerfile for the Docker image.
- `Makefile`: A Makefile for building the Docker image.

Optionaly the subdirectory may contain the following files and folders:  

- `README.md`: A README file for the Docker image.
- `enironment.yml`: A conda environment file for installing dependencies.
- Any resource files (e.g. python script) needed to build the Docker image.

Example Directory Tree:  

```Text
docker
|__lr-bwa
|  |__Dockerfile
|  |__Makefile
|  |__README.md
|__lr-pb
|__lr-ont
```


## Dockerfile Guidelines 
This section outlines the guidelines for creating Dockerfiles, including the format, structure, and best practices for creating efficient and maintainable Docker images.
Docker Docs provides a valuable resource for learning about Dockerfiles. The following 
links provide a good starting point for creating Dockerfiles using general best practices: [Docker Best Practices](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/)
In addition to the Docker best practices, use the following guidelines when creating Dockerfiles for the pipelines repository.

- When appropriate be sure to add a comment proceeding Docker instructions to explain its purpose.
```Dockerfile
# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-pb/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt install -y curl git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash
```

- Specify a `MAINTAINER` for the Docker image.
```Dockerfile
FROM continuumio/miniconda3

MAINTAINER Barbra Mills
```

- Specify version numbers for all packages installed in the Docker image.
```Dockerfile
RUN conda install python=3.6.9
RUN conda create -n venv python=3.6.9 
```


## Image Naming and Tagging Guidelines  
This section outlines the guidelines for naming and tagging Docker images, including the format, structure, and best practices for creating consistent and descriptive image names and tags.

* Use descriptive names: Choose a name that clearly identifies the image and its purpose. Avoid using generic names like "docker-image" or "latest".
* Name should match directory name: When possible the name of the Docker image should match the name of the docker subdirectory it is located in.
* Use lowercase letters: Docker image names should be in lowercase letters.
* Use semantic versioning: Follow the semantic versioning pattern (major.minor.patch) to ensure consistency and compatibility between different versions of the image.
* Avoid special characters: Avoid using special characters in the image name or tag, as it may cause issues with some systems or platforms.
 

## Testing and CI/CD Guidelines 
_This section should provide guidelines for testing Docker images and containers, including best practices for creating automated tests and integrating with CI/CD pipelines to ensure consistent builds and deployments.
TBD: This section is still under development._
