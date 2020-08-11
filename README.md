# CUT&Tag analysis pipeline
#### [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)
##### To run on the Imperial College [HPC](https://wiki.imperial.ac.uk/display/HPC/High+Performance+Computing)


- Trimming
- Alignment
- Peak Calling
- Motif Finding
- Cut Matrix Generation
- Motif Footprinting Steps

#### Trimming (and QC)

 


#### Alignment

Two alignments will be run to align the human DNA and carry-over E.coli DNA which will be used later for sample calibration.

1. Align reads to the reference **human** genome (hg19)
2. Align reads to the reference **E.coli** genome (strain K12, substrain MG1655)

##### Align human reads:

Align paired-end reads to hg19 using bowtie2 and the following arguments:

`--end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -l 10 -X 700`

The E.coli reference genome was downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3?report=fasta) for the K-12 substr. MG1655. Bowtie2 indexes were created and are located at the path:

`ls /rds/general/user/hm1412/projects/cebolalab_liver_regulomes/live/reference-genomes`

##### Align E.coli reads:

E. coli carry-over fragments are algined to the NCBI Ecoli genome (Escherichia coli str. K12 substr. MG1655 U00096.3) using:

`--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant  --phred33 -I 10 -X 700'

