#!/bin/bash
############### build samtools and bcftools from source
SAMTOOLS_VERSION=1.13
BCFTOOLS_VERSION=1.13
PREREQ_PACKAGES="zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev"
EXTRA_PACKAGES="autoconf automake make gcc bzip2 ca-certificates wget"

set -euo pipefail

apt-get update --fix-missing
export DEBIAN_FRONTEND=noninteractive
apt-get -qqy install --no-install-recommends ${PREREQ_PACKAGES} ${EXTRA_PACKAGES}

wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd samtools-${SAMTOOLS_VERSION}
./configure --without-curses --enable-libcurl
make -s all all-htslib
make install install-htslib
cd ..
rm -rf samtools-${SAMTOOLS_VERSION}* # removes the tar file and the directory into which it was extracted

wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2
cd bcftools-${BCFTOOLS_VERSION}
./configure --without-curses
make -s all all-htslib
make install
cd ..
rm -rf bcftools-${BCFTOOLS_VERSION}* # removes the tar file and the directory into which it was extracted

mkdir /install
tar -C / -czf /install/install.tgz usr/local

cat <<EOF >/install/install.sh
#!/bin/bash
tar -C / -xzf /install/install.tgz
apt-get update --fix-missing
export DEBIAN_FRONTEND=noninteractive
apt-get -qqy install --no-install-recommends ${PREREQ_PACKAGES}
apt-get -qqy clean
rm -rf /tmp/* \
       /var/tmp/* \
       /var/cache/apt/* \
       /var/lib/apt/lists/* \
       /usr/share/man/?? \
       /usr/share/man/??_*
rm -rf /install
EOF

chmod a+x /install/install.sh
rm build.sh
