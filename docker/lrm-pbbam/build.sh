#!/bin/bash

############### build pbbam from source
BOOST_VERSION="1_81_0"
BOOST_DIR=$(echo ${BOOST_VERSION} | sed 's/_/./g')
GTEST_VERSION="v1.14.0"
PBBAM_VERSION="v2.4.0"
PACKAGES="g++ make cmake meson ninja-build git wget bzip2 ca-certificates python3"

set -euo pipefail

apt-get update --fix-missing
export DEBIAN_FRONTEND=noninteractive
apt-get -qqy install --no-install-recommends ${PACKAGES}
ln -s /usr/bin/python3 /usr/bin/python

wget --no-verbose https://boostorg.jfrog.io/artifactory/main/release/${BOOST_DIR}/source/boost_${BOOST_VERSION}.tar.bz2
tar xjf boost_${BOOST_VERSION}.tar.bz2
cd boost_${BOOST_VERSION}
./bootstrap.sh --prefix=/usr/local
./b2 --prefix=/usr/local install --without-python
cd ..
rm -rf boost_${BOOST_VERSION}.tar.bz2 boost_${BOOST_VERSION}

git clone https://github.com/google/googletest.git -b ${GTEST_VERSION}
cd googletest
mkdir build
cd build
cmake .. -DBUILD_GMOCK=OFF
make install
cd ../..
rm -rf googletest

ln -s /usr/local/lib /usr/local/lib64 # meson puts libs in /usr/local/lib64, which isn't searched by ld
git clone https://github.com/PacificBiosciences/pbbam.git -b ${PBBAM_VERSION}
cd pbbam
meson setup build
cd build
ninja install
cd ../..
rm -rf pbbam

apt-get -qqy purge ${PACKAGES}
apt-get -qqy autoremove
apt-get -qqy clean
rm -rf /tmp/* \
       /var/tmp/* \
       /var/cache/apt/* \
       /var/lib/apt/lists/* \
       /usr/share/man/?? \
       /usr/share/man/??_*

rm -rf /usr/local/include/boost
mkdir -p /install
tar -C / -czf /install/install.tgz usr/local
rm -rf /usr/local
cat <<EOF >/install/install.sh
#!/bin/bash
tar -C / -xzf /install/install.tgz
ldconfig
rm -rf /install
EOF
chmod a+x /install/install.sh
rm build.sh
