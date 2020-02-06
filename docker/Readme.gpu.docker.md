![image](https://cloud.githubusercontent.com/assets/3028125/12213714/5b208976-b632-11e5-8406-38d379ec46aa.png)

If you are building a docker (i.e. official docker doesn't exist) of tools that comes with GPU accelerations (nice!), then it takes some non-trivial effort to build that. Here we briefly explain how to save you some time.

First, you need a machine to with (NVDA) GPU to test run your docker (don't test with your production runs). So a cloud VM instance will be enough, if you don't have a physical one close. 

Below we assume you have such a VM with Ubuntu 18.04 (hence all the `sudo`s).

## Install CUDA driver

You can run the following to see if which CUDA drivers are available for your GPU:

```bash
sudo apt install ubuntu-drivers-common
sudo ubuntu-drivers devices
```

Then select automatically or manually

```bash
sudo ubuntu-drivers autoinstall # auto
sudo apt install nvidia-<VERSION> # manual
```

## Reboot, and Install Docker engine
Yes, reboot after the driver installation.
Then running `nvidia-smi` should confirm that driver is installed for your GPU(s).

Now we need to install a Docker engine no earlier than version 19.03 (so you don't need to go through the extra steps using nvidia-docker2 packages.)

```bash
sudo apt-get remove docker docker-engine docker.io containerd runc
sudo apt-get update
sudo apt-get -qqy install \
                  apt-transport-https \
                  ca-certificates \
                  curl \
                  gnupg-agent \
                  software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
                        $(lsb_release -cs) \
                        stable"
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

Or run

```bash
apt-cache madison docker-ce
```
to select the appropriate version, e.g.

```bash
sudo apt-get install docker-ce=5:19.03.5~3-0~ubuntu-bionic docker-ce-cli=5:19.03.5~3-0~ubuntu-bionic containerd.io
```

## Install Nvidia container tookit and runtime

```bash
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | sudo tee /etc/apt/sources.list.d/nvidia-docker.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit

sudo apt-get install nvidia-container-runtime
# make sure there's GPU support
docker run --help | grep -i gpus
```

## Docker Engine setup

To register the Nvidia runtime, here we select one approach out of [many](https://github.com/NVIDIA/nvidia-container-runtime#docker-engine-setup)

```bash
sudo mkdir -p /etc/systemd/system/docker.service.d
sudo tee /etc/systemd/system/docker.service.d/override.conf <<EOF
[Service]
ExecStart=
ExecStart=/usr/bin/dockerd --host=fd:// --add-runtime=nvidia=/usr/bin/nvidia-container-runtime
EOF
sudo systemctl daemon-reload
sudo systemctl restart docker
```

## Test run

You should be good to go after this

```bash
sudo docker pull nvidia/cuda:9.0-base
sudo docker run --gpus all nvidia/cuda:9.0-base nvidia-smi
```

Or better yet, if you tool has tests distributed with it, you should run such tests after you've successfully built the image:

```bash
sudo docker run --gpus all <COMMAND_TO_RUN_TESTS>
```