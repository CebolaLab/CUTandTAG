# CUT&Tag analysis pipeline
#### [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)

The pipeline covers the following general steps:

- Alignment 
- Peak Calling
- Motif Finding
- Cut Matrix Generation
- Motif Footprinting Steps

All required programs required have been installed and are available in the CebolaLab [CUTandTAG anaconda environment](https://github.com/CebolaLab/CUTandTAG/tree/master/anaconda-env). If you are using anaconda, you can copy this into your own anaconda environments. On the Imperial College HPC, for example, this is located at `/rdsgpfs/general/user/"$(whoami)"/home/anaconda3/envs/`. The pipeline can be run without installing anaconda by directing the necessary scripts to the `bin` of the downloaded CUTandTAG environment (download the CUTandTAG directory and save it e.g. to your home directory). This will be described in the following steps.

### Alignment

Two alignments will be run to align the human DNA and carry-over E.coli DNA which will be used later for sample calibration. The alignment parameters are run according to CUT&Tag authors (see the [pipeline](https://www.protocols.io/view/cut-amp-tag-home-bd26i8he?step=50)). The authors recommend to skip adapter trimming and to run the alignments using bowtie2 with the below parameters, which should result in accurate read alignment. Two alignments are carried out:

1. Align reads to the reference **human** genome (hg19)
2. Align reads to the reference **E.coli** genome (strain K12, substrain MG1655)

**Human genome:** The UCSC hg19 ***masked*** reference genome was [downloaded](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/) and [indexed using bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome). 

**E.coli genome:** The E.coli reference genome for the strain K-12 strain, MG1655 substrain was [downloaded](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3?report=fasta) from UCSC and was also [indexed](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome) using bowtie2. 

For those with access to the Imperial College HPC and the Cebola Lab project space, the reference genomes and index files are available at this path:

`/rds/general/user/"$whoamI"/projects/cebolalab_liver_regulomes/live/reference-genomes/` 

The alignments are carried out using bowtie2 with the below arguments. An example script is available [here](https://github.com/CebolaLab/CUTandTAG/blob/master/alignment.sh).

##### Align human reads:

`--end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -l 10 -X 700`

##### Align E.coli reads:

`--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant  --phred33 -I 10 -X 700`

##### Sort bam file

The output bam files must be sorted, here by coordinate. `"$base"` again refers to the your filename/sample ID: 

`picard SortSam I="$base".bam O="$base"-sorted.bam SO=coordinate CREATE_INDEX=TRUE`

### Peak Calling

Peak calling will be carried out using both macs2 and SEACR. First, aligned fragments are extracted from the aligned bam file and are output in bed format using `bedtools bamtobed`. For calibration using the E.coli reads, the bed files require the length of the fragment to be added (as described in the Henikoff lab calibration [script](https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.csh)). Assuming the output file is in the format `"$base"_aligned_reads.bam` (where "$base" is your sample identifier), the following code may be used:

`bedtools bamtobed -bedpe -i "$base"-sorted.bam | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > "$base".bed`

#### Calibration 

If you are working with multiple samples, e.g. a sample and a control, they should be calibrated to be comparable. The calibration is carried out using the Henikoff lab calibration [script](https://github.com/Henikoff/Cut-and-Run/blob/mast\
er/spike_in_calibration.csh)). The calibration effectively scales the mapped counts according to the total number of E.coli reads. The ratio of primary genome to E.coli genome is expected to be the same for all samples. The script is included in this repository (see above). 

**Step 1:** Filter fragments to be within a minimum and maximum length



**Steps 2:** Run the calibration on the filtered file

To run it, seven arguments are required:

`spike_calibrate.csh genome.bed spike_genome.bed scale output(bg|bga|d) genome_chr_lens_file  min_len max_len`

- **genome.bed** the converted bed file (hg19 alignment), `"$base".bed`
- **spike_genome.bed** the converted bed file (E.coli alignment)
- **scale** 
- **output(bg|bga|d)** bg = BedGraph, bga = BedGraph including regions with 0 coverage, d = depth at each genome position with 1-based coordinates. 
- **genome_chr_lens_file**
- **min_len** minumum fragment length, `min=$(cut -f 11 SRR8383480_aligned_reads.bed | sort | uniq | head -1)`
- **max_len** maximum fragment length, `max=$(cut -f 11 SRR8383480_aligned_reads.bed | sort | uniq | tail -1)'
 
As described in the original CUT&Tag paper ([Kaya-Okur et al. 2019](https://www.nature.com/articles/s41467-019-09982-5#data-availability)), macs2 will be used with the following parameters:

`macs2 callpeak -t input_BED -f BEDPE -p 1e-5 --keep-dup all -n output_prefix` 