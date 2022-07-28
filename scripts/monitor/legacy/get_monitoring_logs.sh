#!/usr/bin/env bash

d=$( echo "${1}" | sed 's@/ *$@@' )

TMP_FILE=monitoring_logs.txt
OUT_DIR=monitoring_logs

if [[ $# -ne 1 ]] ; then
  echo "ERROR: you must give a google cloud directory." 1>&2
  exit 1
elif [[ ! ${d} =~ ^gs:// ]] ; then
  echo "ERROR: path is not a google bucket: ${d}" 1>&2
  exit 1
fi

if [[ -d ${OUT_DIR} ]] ; then
  echo "ERROR: output directory already exists: ${OUT_DIR}" 1>&2
  exit 1
fi

mkdir -p ${OUT_DIR}

echo "Getting monitoring.log file list from workflow run: ${d} ..."
gsutil -m ls ${d}/**/monitoring.log > ${TMP_FILE} 

numfiles=$(wc -l ${TMP_FILE} | awk '{print $1}')
echo "Discovered ${numfiles} monitoring logs."

echo -n "Copying files ..."

i=0
ps_limit=50
while read url ; do 
    runningJobs=$(jobs | wc -l | xargs)
    while [ $runningJobs -gt $ps_limit ] ; do
    	sleep 1
		runningJobs=$(jobs | wc -l | xargs)
	done
    nn=$(echo ${url} | sed -e "s@${d}/@@" -e 's@shard-\([0-9]\)/@shard-0\1/@g' -e 's@shard-\([0-9][0-9]\)/@shard-0\1/@g' -e 's@/@_@g')
    gsutil -m cp $url ${OUT_DIR}/${nn} &> /dev/null & 
   	echo -n "."
	let i=i+1
done < ${TMP_FILE}
echo

echo "Waiting for remaining copying processes to finish..."
wait

echo "Output stored in: ${OUT_DIR}"

