#!/bin/bash

# specificaly avoiding these as not all file tyes may not be there, and that could cause problems with gsutl
# set -eu
# set -o pipefail

# file_types=( .vcf .vcf.gz .tbi .html .txt .fasta .fasta.gz .fa .fa.gz .gfa .gfa.gz )

################################################################################
show_help() {
cat << EOF
Grant READER access to the provided email address, to certain file types under a prefix of a bucket.

Syntax
  grant_access_to_collab.sh [Options]

  Options:
    -h or --help:
      display help and exit
    -u or --user [USER]:
      email address to grant access to
    -b --bucket [BUCKET]
      the GCS bucket, without the leading gs://
    -p --prefix [PREFIX]:
      prefix under the bucket (i.e. what's commonly referred to a sub-directories)
    -e --ext [FILE_TYPE_LIST_FILE]
      a file holding type of files (as in extension) to be granted access to
EOF
}
################################################################################
throw_error() {
    echo "$1" >&2
    exit 1
}
################################################################################
if [ $# -lt 1 ]; then
  show_help
  exit 1
fi
################################################################################
while [ $# -ge 1 ]; do
    case $1 in
        -h|-\?|--help)
            show_help
            exit 0
            ;;
        -u|--user)
            if [ $# -ge 2 ]; then
                user_email="$2"
                shift 2
            else
                throw_error "--user requires a non-empty argument"
            fi
            ;;
        --user=?*)
            user_email=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        -b|--bucket)
            if [ $# -ge 2 ]; then
                gcs_bucket="$2"
                shift 2
            else
                throw_error "--bucket requires a non-empty argument"
            fi
            ;;
        --bucket=?*)
            gcs_bucket=${1#*=} # remove everything up to = and assign rest to init
            shift
            ;;
        -p|--prefix)
            if [ $# -ge 2 ]; then
                prefix="$2"
                shift 2
            else
                throw_error "--prefix requires a non-empty argument"
            fi
            ;;
        --prefix=?*)
            prefix=${1#*=} # remove everything up to = and assign rest to save
            shift
            ;;
        -e|--ext)
            if [ $# -ge 2 ]; then
                pattern_file="$2"
                shift 2
            else
                throw_error "--bucket requires a non-empty argument"
            fi
            ;;
        --ext=?*)
            pattern_file=${1#*=} # remove everything up to = and assign rest to init
            shift
            ;;
        --)   # explicit call to end of all options
            shift
            break
            ;;
        -?*)  # unsupported option
            throw_error "Unknown option \"$1\". use --help for syntax"
            ;;
        *)  # not an option, a positional argument. break out
            break
            ;;
    esac
done

################################################################################
access_level='READ'

if [[ -z ${user_email+x} ]]; then
    echo "user_email isn't set or is set to empty" && exit 1
elif [[ -z ${gcs_bucket+x} ]]; then
    echo "gcs_bucket isn't set or is set to empty" && exit 1
elif [[ -z ${prefix+x} ]]; then
    echo "prefix isn't set or is set to empty" && exit 1
elif [[ -z ${access_level+x} ]]; then
    echo "access_level isn't set or is set to empty" && exit 1
elif [[ -z ${pattern_file+x} ]]; then
    echo "access_level isn't set or is set to empty" && exit 1
fi
################################################################################
IFS=$'\n' read -d '' -r -a file_types < "${pattern_file}"
s=$(printf '%s ' "${file_types[@]}")
echo "File types to share: ${s}"
################################################################################
gcs_bucket="${gcs_bucket%/}"
prefix="${prefix%/}"
################################################################################
for t in "${file_types[@]}";
do
    echo -e "\nGranting user ${user_email} with ${access_level}-level access to file type ${t} under gs://${gcs_bucket}/${prefix}/"
    gsutil -mq acl \
        ch -u "${user_email}":"${access_level}" \
        -R "gs://${gcs_bucket}/${prefix}/**${t}"
    sleep 10 # avoid being blocked
done

echo "DONE DONE DONE DONE"
