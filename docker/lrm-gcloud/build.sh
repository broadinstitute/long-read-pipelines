#!/bin/bash

PACKAGES="apt-transport-https ca-certificates gnupg curl"

############### install gcloud
set -euo pipefail

apt-get update --fix-missing
export DEBIAN_FRONTEND=noninteractive
apt-get -qqy install --no-install-recommends ${PACKAGES}
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" > /etc/apt/sources.list.d/google-cloud-sdk.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
apt-get -qqy update
apt-get -qqy install --no-install-recommends google-cloud-cli

mkdir /install
tar -C /usr/lib/google-cloud-sdk -czf /install/install.tgz .
cat <<EOF >/install/install.sh
#!/bin/bash
tar -C /usr/local -xzf /install/install.tgz
ldconfig
gcloud config set core/disable_usage_reporting true
gcloud config set component_manager/disable_update_check true
gcloud config set metrics/environment github_docker_image
rm -rf /install
EOF
chmod a+x /install/install.sh

apt-get -qqy purge ${PACKAGES} google-cloud-cli
apt-get -qqy autoremove
apt-get -qqy clean
rm -rf /tmp/* \
       /var/tmp/* \
       /var/cache/apt/* \
       /var/lib/apt/lists/* \
       /usr/share/man/?? \
       /usr/share/man/??_*

rm build.sh
