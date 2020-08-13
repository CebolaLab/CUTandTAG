
genome_bed=$1
spike_bed=$2
scale=$3
report=$4
genome_len=$5
min_len=$6
max_len=$7

temp=`basename $genome_bed .bed`
name=`dirname $genome_bed`

if [ -f $spike_bed ] && [ -s $spike_bed ]; then 
    spike_count=$(wc -l $spike_bed | awk '{print $1}')
    scale_factor=$(echo "$scale/$spike_count" | bc -l)
else
    spike_count=0
    scale_factor=$scale
fi

echo 'scale_factor =' $scale_factor 


if [ $report == "bg" ] || [ $report == "bga" ]; then
    output=$temp.bedgraph
    echo "track type=bedGraph name=$temp" > $output
fi
#Select fragments within the length range, assumes fragment length is in column 4 of the bed file
#cat $genome_bed | awk -v min=$min_len -v max=$max_len '{if ($4 >= min && $4 <= max) print}' > $$.temp.bed

bedtools genomecov -$report -scale $scale_factor -i $genome_bed -g $genome_len >> $output
