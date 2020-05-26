#!/bin/bash

set -euxo pipefail

batch_size=10,000,000,000
racon_batch_size=$(echo $batch_size | sed 's/,//g')
export RACON_PARAMS="-u -t 32 -c 4 --cudaaligner-batches 50 --split ${racon_batch_size}"

docker_id=$(sudo docker images | grep -F 'racon' | awk '{print $3}')
export dir="<my_data_dir>"

echo "POLISHING FATHER"
export READS="/data/<father_reads>"
export ASM="/data/<father_asm>"
export OVP="/data/<father_paf>"
export OUTPUT="<father_posh.fa>"
sudo docker run \
    --gpus all \
    -v ${dir}:/data \
    "${docker_id}" \
    bash polish.racon.wrapper.sh "${READS}" "${OVP}" "${ASM}" "${RACON_PARAMS}" \
    > "${OUTPUT}"

echo "POLISHING MOTHER"
export READS="/data/<mother_reads>"
export ASM="/data/<mother_asm>"
export OVP="/data/<mother_paf>"
export OUTPUT="<mother_posh.fa>"
sudo docker run \
    --gpus all \
    -v ${dir}:/data \
    "${docker_id}" \
    bash polish.racon.wrapper.sh "${READS}" "${OVP}" "${ASM}" "${RACON_PARAMS}" \
    > "${OUTPUT}"
