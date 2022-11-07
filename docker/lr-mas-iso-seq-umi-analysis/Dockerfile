# Start with a good base python3 image:
FROM ubuntu:20.04
MAINTAINER Jonn Smith <jonn@broadinstitute.org>

# Make sure we don't need to interact with any package installations:
ARG DEBIAN_FRONTEND=noninteractive

# Set the working directory to /
WORKDIR /

################################################################################
# Install system packages:

# Make sure we can actually talk to package repos:
RUN apt-get update && apt-get install -y apt-utils && apt-get -y upgrade && add-apt-repository ppa:deadsnakes/ppa && apt-get -y update
RUN apt-get -y install python3.7

# Install our scripts:
RUN mkdir /scripts
COPY scripts/* /scripts/
RUN chmod +x /scripts/*.py /scripts/*.sh

################################################################################
# Final runtime configuration:
# Let's start at the root:
WORKDIR /

