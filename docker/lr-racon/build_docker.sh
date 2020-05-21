#!/bin/bash

set -eu

me=$(basename "$0")
show_help() {
cat << EOF

Build GPU-enabled docker on Ubuntu 18.04

Syntax
  sudo bash ${me} [Options]
  Options:
    -h or --help:
      display help and exit
    -f or --df:
      absolute path to dockerfile (note its residing dir will be the building context)
    -i or --image:
      image name
    -t or --tag:
      tag to build
    -x or --test:
      test command to be executed after successful build (surround by double quote)

NOTE:
  1. Please run the following first and reboot your VM, if you haven't already.
     Otherwise the build will fail.
         sudo apt install ubuntu-drivers-common
         sudo ubuntu-drivers autoinstall
  2. We assume the script is executed as superuser.

EOF
}

if [[ "$#" -eq 0 ]]; then
    show_help && exit 0;
fi

while [ $# -ge 1 ]; do
    case $1 in
        -h|-\?|--help)
            show_help
            exit 0
            ;;
        -f|--df)
            if [ $# -ge 2 ]; then
                DOCKER_FILE="$2"
                shift 2
            else
                throw_error "--df requires a non-empty argument"
            fi
            ;;
        --df=?*)
            DOCKER_FILE=${1#*=} # remove everything up to = and assign rest to df
            shift
            ;;
        -i|--image)
            if [ $# -ge 2 ]; then
                IMAGE="$2"
                shift 2
            else
                throw_error "--image requires a non-empty argument"
            fi
            ;;
        --image=?*)
            IMAGE=${1#*=} # remove everything up to = and assign rest to image
            shift
            ;;
        -t|--tag)
            if [ $# -ge 2 ]; then
                TAG="$2"
                shift 2
            else
                throw_error "--tag requires a non-empty argument"
            fi
            ;;
        --tag=?*)
            TAG=${1#*=} # remove everything up to = and assign rest to tag
            shift
            ;;
        -x|--test)
            if [ $# -ge 2 ]; then
                TEST_CMD="$2"
                shift 2
            else
                throw_error "--test requires a non-empty argument"
            fi
            ;;
        --test=?*)
            TEST_CMD=${1#*=} # remove everything up to = and assign rest to test
            shift
            ;;
        -?*)  # unsupported option
            throw_error "Unknown option \"$1\". use --help for syntax"
            ;;
        *)  # not an option, a positional argument. break out
            break
            ;;
    esac
done

if [[ -z "${DOCKER_FILE}" ]] || [[ -z "${IMAGE}" ]] || [[ -z "${TAG}" ]] ; then
    echo "Some or all of the requested arguments is not set properly. Please check."

    echo "DOCKER_FILE=${DOCKER_FILE}"
    echo "IMAGE=${IMAGE}"
    echo "TAG=${TAG}"
    exit 1;
fi

if [[ ! -f "${DOCKER_FILE}" ]] ; then
    echo "Provided docker file not readable. Please check:"
    echo "DOCKER_FILE=${DOCKER_FILE}"
    exit 1;
fi

BUID_CONTEXT=$(dirname "${DOCKER_FILE}")
if [[ -z ${BUID_CONTEXT+x} ]]; then
    echo "Build context is not found"
    exit 1
fi

echo -e '\033[0;93m################################################################################\033[0m'
echo -e '\033[0;93mINSTALL DOCKER ENGINE\033[0m'
apt-get remove docker docker-engine docker.io containerd runc || true
apt-get -qqy update
apt-get -qqy install \
             apt-transport-https \
             ca-certificates \
             curl \
             gnupg-agent \
             software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add -
add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
                    $(lsb_release -cs) \
                    stable"
apt-get -qqy update
apt-get -qqy install docker-ce docker-ce-cli containerd.io

echo -e '\033[0;93m################################################################################\033[0m'
echo -e '\033[0;93mINSTALL NVIDIA CONTAINER TOOKIT AND RUNTIME\033[0m'
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | tee /etc/apt/sources.list.d/nvidia-docker.list
apt-get -qqy update
apt-get -qqy install nvidia-container-toolkit nvidia-container-runtime

echo "make sure there's GPU support"
docker run --help | grep -i gpus

echo -e '\033[0;93m################################################################################\033[0m'
echo -e '\033[0;93mDOCKER ENGINE SETUP\033[0m'
mkdir -p /etc/systemd/system/docker.service.d
tee /etc/systemd/system/docker.service.d/override.conf <<EOF
[Service]
ExecStart=
ExecStart=/usr/bin/dockerd --host=fd:// --add-runtime=nvidia=/usr/bin/nvidia-container-runtime
EOF
systemctl daemon-reload
systemctl restart docker

echo -e '\033[0;93m################################################################################\033[0m'
echo -e '\033[0;93mBUILD DOCKER\033[0m'
docker build -f "${DOCKER_FILE}" -t "${IMAGE}":"${TAG}" "${BUID_CONTEXT}"

echo -e '\033[0;93m################################################################################\033[0m'
if [[ -z ${TEST_CMD+x} ]]; then
    echo -e '\033[0;93mNO TEST COMMAND FOUND. DONE.\033[0m'
else
    echo -e '\033[0;93mEXECUTING PROVIDED TEST COMMAND\033[0m'
    docker run --gpus all "${IMAGE}":"${TAG}" bash -c "${TEST_CMD}"
fi
