output=$(echo $1 | sed s/.bed//).sorted.bed
#echo $output
bedtools sort -i $1 > $output
igvtools index $output
