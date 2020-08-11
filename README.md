# CUT&Tag analysis pipeline
#### [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)
##### To run on the Imperial College [HPC](https://wiki.imperial.ac.uk/display/HPC/High+Performance+Computing)

The pipeline covers the following general steps:

- Alignment 
- Peak Calling
- Motif Finding
- Cut Matrix Generation
- Motif Footprinting Steps

### Alignment

Two alignments will be run to align the human DNA and carry-over E.coli DNA which will be used later for sample calibration. The recommended steps (see [here](https://www.protocols.io/view/bench-top-cut-amp-tag-bcuhiwt6?step=69)) are to *not* run adapter trimming, but to instead run the alignment with the below parameters, which should result in accurate read alignment:

1. Align reads to the reference **human** genome (hg19)
2. Align reads to the reference **E.coli** genome (strain K12, substrain MG1655)

The reference genomes and Bowtie2 indexes are available at the following paths:

Human genome, UCSC hg19 (obtained [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/))

`/rds/general/user/hm1412/projects/cebolalab_liver_regulomes/live/reference-genomes/hg19.masked`

E.coli, K-12 strain, MG1655 substrain (obtained [here](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3?report=fasta))

`/rds/general/user/hm1412/projects/cebolalab_liver_regulomes/live/reference-genomes/E.coli` 

The two alignments use different arguments, shown below:

##### Align human reads:

Align paired-end reads to hg19 using bowtie2 and the following arguments:

`--end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -l 10 -X 700`

##### Align E.coli reads:

E. coli carry-over fragments are algined to the NCBI Ecoli genome (Escherichia coli str. K12 substr. MG1655 U00096.3) using:

`--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant  --phred33 -I 10 -X 700`

