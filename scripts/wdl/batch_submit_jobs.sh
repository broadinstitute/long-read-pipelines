#!/bin/bash

command -v cromshell >/dev/null 2>&1 || { echo >&2 "I require cromshell but it's not installed.  Aborting."; exit 1; }

me=`basename "$0"`
show_help() {
cat << EOF
Run workflow on a given list of input jsons.

Syntax
  ${me} [Options]
  Options:
    -h or --help:
      display help and exit
    -w or --wdl:
      WDL containing the workflow definition
    -i or --input:
      a file containing the relative or absolute path to input json files
    -s or --option:
      workflow options json file
    -p or --dependency:
      workflow dependency zip file
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
        -w|--wdl)
            if [ $# -ge 2 ]; then
                WORKFLOW_WDL="$2"
                shift 2
            else
                throw_error "--wdl requires a non-empty argument"
            fi
            ;;
        --wdl=?*)
            WORKFLOW_WDL=${1#*=} # remove everything up to = and assign rest to wdl
            shift
            ;;
        -i|--input)
            if [ $# -ge 2 ]; then
                INPUT_JSON_LIST="$2"
                shift 2
            else
                throw_error "--input requires a non-empty argument"
            fi
            ;;
        --input=?*)
            INPUT_JSON_LIST=${1#*=} # remove everything up to = and assign rest to input
            shift
            ;;
        -s|--option)
            if [ $# -ge 2 ]; then
                WORKFLOW_OPTION_JSON="$2"
                shift 2
            else
                throw_error "--option requires a non-empty argument"
            fi
            ;;
        --option=?*)
            WORKFLOW_OPTION_JSON=${1#*=} # remove everything up to = and assign rest to option
            shift
            ;;
        -p|--dependency)
            if [ $# -ge 2 ]; then
                WORKFLOW_WDL_DEPENDENCY="$2"
                shift 2
            else
                throw_error "--dependency requires a non-empty argument"
            fi
            ;;
        --dependency=?*)
            WORKFLOW_WDL_DEPENDENCY=${1#*=} # remove everything up to = and assign rest to dependency
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


if [[ -z "${WORKFLOW_WDL}" ]] || [[ -z "${WORKFLOW_OPTION_JSON}" ]] || [[ -z "${WORKFLOW_WDL_DEPENDENCY}" ]] || [[ -z "${INPUT_JSON_LIST}" ]] ; then
    echo "Some or all of the requested arguments is not set properly. Please check."

    echo "WORKFLOW_WDL=${WORKFLOW_WDL}"
    echo "WORKFLOW_OPTION_JSON=${WORKFLOW_OPTION_JSON}"
    echo "WORKFLOW_WDL_DEPENDENCY=${WORKFLOW_WDL_DEPENDENCY}"

    echo "INPUT_JSON_LIST=${INPUT_JSON_LIST}"
    exit 1;
fi

if [[ ! -f "${WORKFLOW_WDL}" ]] || [[ ! -f "${WORKFLOW_OPTION_JSON}" ]] || [[ ! -f "${WORKFLOW_WDL_DEPENDENCY}" ]] || [[ ! -f "${INPUT_JSON_LIST}" ]] ; then
    echo "Some or all of the files are not reada. Please check."

    echo "WORKFLOW_WDL=${WORKFLOW_WDL}"
    echo "WORKFLOW_OPTION_JSON=${WORKFLOW_OPTION_JSON}"
    echo "WORKFLOW_WDL_DEPENDENCY=${WORKFLOW_WDL_DEPENDENCY}"

    echo "INPUT_JSON_LIST=${INPUT_JSON_LIST}"
    exit 1;
fi

keep_submitting_and_alias() {

    nick_name=$1

    shift 1;

    echo "=================================================================="
    echo "$@"

    # keep trying until got submited
    echo "Cannot communicate with Cromwell server" > temp.stderr
    while  grep -qF 'Cannot communicate with Cromwell server' temp.stderr || \
           grep -qF 'job was not properly submitted' temp.stderr
    do
        cromshell \
            submit "$@" \
            > temp.stdout \
            2> temp.stderr
        sleep 10
    done

    # alias cryptic workflow ID with user provided nick name
    workflow_id=$(grep -F 'id' temp.stdout | grep -F 'status' | grep -F 'Submitted' | jq '.id' | sed 's/"//g')
    echo "${nick_name}: ${workflow_id}"
    cromshell \
        alias \
        "${workflow_id}" \
        "${nick_name}"
    rm temp.stdout temp.stderr
    echo "=================================================================="
}
export -f keep_submitting_and_alias

for input_json in `cat ${INPUT_JSON_LIST}`; do
    filename=$(basename -- "${input_json}")
    timestamp=$(date '+%Y%m%d%H%M%S')
    nick="${filename%.*}.${timestamp}"

    keep_submitting_and_alias \
        "${nick}" \
        "${WORKFLOW_WDL}" \
        "${input_json}" \
        "${WORKFLOW_OPTION_JSON}" \
        "${WORKFLOW_WDL_DEPENDENCY}"
done
