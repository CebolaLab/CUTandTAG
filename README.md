# CUT&Tag analysis pipeline
#### [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)

The pipeline covers the following general steps:

- Alignment 
- Peak Calling
- Motif Finding
- Cut Matrix Generation
- Motif Footprinting Steps

The programs required have been installed and are available in the CebolaLab CUTandTAG anaconda environment. If you are running the pipeline from the Imperial College HPC, you can copy the anaconda working environment from the Cebola Lab using the following command:

`cp /rdsgpfs/general/project/cebolalab_liver_regulomes/live/CUTandTAG/CUTandTAG.shared/anaconda-env/CUTandTAG /rdsgpfs/general/user/"$(whoami)"/home/anaconda3/envs/`

The environment is also downloadable from this github and should be moved to your own `anaconda/envs/` folder.

### Alignment

Two alignments will be run to align the human DNA and carry-over E.coli DNA which will be used later for sample calibration. The [recommended steps](https://www.protocols.io/view/bench-top-cut-amp-tag-bcuhiwt6?step=69)) are to skip adapter trimming and instead run the alignment with the below parameters, which should result in accurate read alignment:

1. Align reads to the reference **human** genome (hg19)
2. Align reads to the reference **E.coli** genome (strain K12, substrain MG1655)

The UCSC hg19 ***masked*** reference genome was [downloaded](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/) and [indexed using bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome). 

The E.coli reference genome for the strain K-12 strain, MG1655 substrain was obtained [here](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3?report=fasta) and was also [indexed](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome) using bowtie2. For those working on the Imperial HPC, the reference genome and Bowtie2 indexes are available at the following path:

For those with access to the Imperial College HPC and the Cebola Lab project space, the reference genomes and index files are available at this path:

`/rds/general/user/"$whoamI"/projects/cebolalab_liver_regulomes/live/reference-genomes/` 

The two alignments use different arguments, shown below:

##### Align human reads:

Align paired-end reads to hg19 using bowtie2 and the following arguments:

`--end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -l 10 -X 700`

##### Align E.coli reads:

E. coli carry-over fragments are algined to the NCBI Ecoli genome (Escherichia coli str. K12 substr. MG1655 U00096.3) using:

`--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant  --phred33 -I 10 -X 700`

