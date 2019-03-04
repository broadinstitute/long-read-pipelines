truth=$1
calls=$2

gangstr20_intervals=$(cat ~/p1/results/gangstr20_confidence_region_intersection.bed | wc -l)

truthset_size=$(bcftools stats ${truth}  |& grep 'number of records' | grep SN | cut -f 4-)
callset_size=$(bcftools stats ${calls}  |& grep 'number of records' | grep SN | cut -f 4-)
total_intervals=$gangstr20_intervals

false_positives=$(bedtools window -w 10 -a ${calls} -b ${truth} -v | wc -l)
false_negatives=$(bedtools window -w 10 -a ${truth} -b ${calls} -v | wc -l)
#true_positives=$(python -c "print(${callset_size} - ${false_positives} - ${false_negatives})")  #$(bedtools window -w 10 -a ${calls} -b ${truth} -u | wc -l)
true_positives=$(bedtools window -w 10 -a ${calls} -b ${truth} -u | wc -l)
true_negatives=$(python -c "print(${gangstr20_intervals} - ${false_positives} - ${false_negatives} - ${true_positives})")  #$(bedtools window -w 10 -a ${calls} -b ${truth} -u | wc -l)
echo Truthset Size: $truthset_size
echo Callset Size: $callset_size
echo Total Intervals: $total_intervals
echo
echo True Positives: $true_positives
echo False Positives: $false_positives
echo False Negatives: $false_negatives
echo True Negatives: $true_negatives
echo

echo Precision = "tp/(tp + fp) = " $(python -c "print(${true_positives} / float(${true_positives} + ${false_positives}))")
echo Recall = "tp/(tp + fn) = " $(python -c "print(${true_positives} / float(${true_positives} + ${false_negatives}))")


#bedtools window -w 10 -a ${calls} -b ${truth} -u | head -n 20