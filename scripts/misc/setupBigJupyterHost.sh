#!/usr/bin/env bash

# This script is designed to be run on a newly-created cloud VM to install all the rerequistes for running a jupyter
# notebook / jupyter lab server.
#
# Author: Jonn Smith

################################################################################
# Docker setup:
sudo apt-get update
sudo apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

sudo apt-key fingerprint 0EBFCD88

sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"

sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io

sudo docker run hello-world

un=$USER

sudo groupadd docker
sudo usermod -aG docker ${un}
newgrp docker
docker run hello-world

################################################################################
# Screenrc file:
rm -f ~/.screenrc
echo "hardstatus alwayslastline                 # Status line is always last" >> ~/.screenrc
echo "" >> ~/.screenrc
echo "#Actual status line:" >> ~/.screenrc
echo "hardstatus string '%{= kG}[ %{G}%H - \$STY %{g}][%= %{=kw}%?%-Lw%?%{r}(%{W}%n*%f%t%?(%u)%?%{r})%{w}%?%+Lw%?%?%= %{g}][ %{B}%Y-%m-%d %{W}%c %{g}]'" >> ~/.screenrc
echo "" >> ~/.screenrc
echo "defscrollback 100000                       # Change scrollback to 30000 lines" >> ~/.screenrc
echo "" >> ~/.screenrc
echo "autodetach on                             # Autodetach session on hangup instead of terminating screen completely" >> ~/.screenrc
echo "startup_message off                       # Turn off the splash screen" >> ~/.screenrc
echo "vbell off                                 # turn off visual bell " >> ~/.screenrc
echo "termcapinfo xterm ti@:te@                 # allow scrolling in xterm" >> ~/.screenrc
echo "" >> ~/.screenrc
echo "#Aliases / keybindings / commands:" >> ~/.screenrc
echo "bind X remove                             # Set ctrl-a Z to close a window inside a screen" >> ~/.screenrc
echo "" >> ~/.screenrc
echo "pow_detach_msg \"Goodbye, \$USER\"" >> ~/.screenrc
echo "" >> ~/.screenrc
echo "#Allow us to use shift to select screens 10 -> 19" >> ~/.screenrc
echo "bind  ! select 11" >> ~/.screenrc
echo "bind  @ select 12" >> ~/.screenrc
echo "bind \# select 13" >> ~/.screenrc
echo "bind  $ select 14" >> ~/.screenrc
echo "bind  % select 15" >> ~/.screenrc
echo "bind \^ select 16" >> ~/.screenrc
echo "bind  & select 17" >> ~/.screenrc
echo "bind  * select 18" >> ~/.screenrc
echo "bind  ( select 19" >> ~/.screenrc
echo "bind  ) select 10" >> ~/.screenrc
echo "" >> ~/.screenrc

