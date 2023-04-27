version 1.0

task CheckFileUpdatedDateGCP {

    meta {
        description: "Checks if file was updated within a specified time (default: 1 day)"
        note: "Specific to GCP"
    }


    input {
        Array[String] file_paths
        Int days_back = 1
        String image_to_use
    }

    #FILE_DATE description: get file info | grep the 'Update Time' row | parse the date info | reformat date info

    command <<<
        set -eu pipefail

        GS_BUCKET_PATHS=("~{sep='" "' file_paths}")
        MINIMUM_DATE=`date -d '~{days_back} day ago' +%Y-%m-%d`
        return_code=0
        echo -e "FileName\tFileUpdatedDate\tMinimumUpdateDate\tValid"

        for GS_FILE in ${GS_BUCKET_PATHS[@]};
        do
          FILE_DATE=`gsutil  stat ${GS_FILE} | grep "Update time" | awk -F '[,][ .]+' '{print $2}' | awk '{print "date -d\""$1FS$2FS$3"\" +%Y-%m-%d"}'| bash`

          if [[ $FILE_DATE >= $MINIMUM_DATE ]] ; then
            echo -e "$GS_FILE\t$FILE_DATE\t$MINIMUM_DATE\tTrue"
          else
            echo "ERROR: Date for $GS_FILE" >&2
            echo "$FILE_DATE is less than $MINIMUM_DATE" >&2

            echo -e "$GS_FILE\t$FILE_DATE\t$MINIMUM_DATE\tFalse"
            return_code=1
          fi
        done

        if [ $return_code == 1 ]; then
            exit 1
        fi

    >>>
    runtime {
        docker: image_to_use
    }
    output {
        File file_date_result = stdout()
    }
}
