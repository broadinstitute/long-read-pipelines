#!/bin/bash
apt-get update --fix-missing
export DEBIAN_FRONTEND=noninteractive
apt-get -qqy install python3-pip
pip install nanoplot
apt-get -qqy clean
rm -rf /tmp/* \
       /var/tmp/* \
       /var/cache/apt/* \
       /var/lib/apt/lists/* \
       /usr/share/man/?? \
       /usr/share/man/??_*
