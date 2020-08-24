#You should have pre-downloaded and indexed versions of the reference human genome (hg19 or other). Set the bt2idx variable to reflect the directory where these files are saved, e.g:
bt2idx=/rds/general/user/"$(whoami)"/projects/cebolalab_liver_regulomes/live/hg19.masked/
base=sampleID #EDIT THIS to your sample ID which prefixes _1.paired.fastq.gz and _2.paired.fastq.gz

#The alignment can be run with pre-installed versions of bowtie2 and samtools, assuming that these are in your $PATH or active conda environment:
(bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x $bt2idx/hg19.masked -1 "$base"_1.paired.fastq.gz -2 "$base"_2.paired.fastq.gz) 2> "$base".bowtie2 | samtools view -bS - > "$base"_aligned_reads.bam
#If you have an average fragment size >700 then edit the -X 700 option to reflect your average fragment size 

#Alternatively, you can direct the program to the bowtie2 and samtools versions downloaded from the Cebola Lab CUTandTAG conda directories. You can adapt the below paths to reflect where the CUTandTAG directory has been saved:
CUTandTAGbin=/rdsgpfs/general/user/"$(whoami)"/home/anaconda3/envs/CUTandTAG/bin

($CUTandTAGbin/bowtie2 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -l 10 -X 700 -x $bt2idx/hg19.masked -1 "$base"_1.paired.fastq.gz -2 "$base"_2.paired.fastq.gz) 2> "$base".bowtie2 | $CUTandTAGbin/samtools view -bS - > $aligndir/"$base"_aligned_reads.bam

#To run bowtie2 on multiple processors/cores, use -p x, where x is the number of processors/cores
