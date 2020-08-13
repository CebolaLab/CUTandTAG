# CUT&Tag analysis pipeline
#### [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)

The pipeline covers the following general steps:

- Alignment 
- Visualisation & Calibration 
- Peak Calling
- Motif Finding
- Cut Matrix Generation
- Motif Footprinting Steps

All required programs required have been installed and are available in the CebolaLab [CUTandTAG anaconda environment](https://github.com/CebolaLab/CUTandTAG/tree/master/anaconda-env). If you are using anaconda, you can copy this into your own anaconda environments. On the Imperial College HPC, for example, this is located at `/rdsgpfs/general/user/"$(whoami)"/home/anaconda3/envs/`. The pipeline can be run without installing anaconda by directing the necessary scripts to the `bin` of the downloaded CUTandTAG environment (download the CUTandTAG directory and save it e.g. to your home directory). This will be described in the following steps.

## Alignment

Two alignments will be run to align the human DNA and carry-over E.coli DNA later to calibrate the samples. The alignment parameters are run according to CUT&Tag authors (see the [pipeline](https://www.protocols.io/view/cut-amp-tag-home-bd26i8he?step=50)). The authors recommend to skip adapter trimming and to run the alignments using bowtie2 with the below parameters, which should result in accurate read alignment. Two alignments are carried out:

1. Align reads to the reference **human** masked genome (hg19) (download [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/))
2. Align reads to the reference **E.coli** genome (strain K12, substrain MG1655) (downloaded [here](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3?report=fasta)).

Both reference genomes should be [indexed](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome) using bowtie2. For those with access to the Imperial College HPC and the Cebola Lab project space, the reference genomes and index files are available at this path:

`/rds/general/user/"$(whoami)"/projects/cebolalab_liver_regulomes/live/reference-genomes/` 

The alignments are carried out using bowtie2 with the below arguments. An example script is available [here](https://github.com/CebolaLab/CUTandTAG/blob/master/alignment.sh).

##### Parameters to align human reads:

`--end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -l 10 -X 700`

##### Parameters to align E.coli reads:

`--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant  --phred33 -I 10 -X 700`

## Visualisation

A bedGraph file will be generated to visualise the data. In order to compare multiple samples, the samples will be calibrated using the carry-over E.coli DNA (the experiment 'spike-in'). In theory, the ratio of primary DNA to E.coli DNA is expected to be the same for each sample. As such, the calibration involves the division of the mapped count divided by the total number of reads aligned to the E.coli genome. The proportion of total DNA reads aligned to the E.coli genome is reported in the `"$(base)".Ecoli.bowtie2` output file. The general steps cover:

1. Sort aligned bam file by read name (queryname)
2. Convert bam to bed
3. Normalise samples using E.coli carry-over DNA
4. Visualise output bedGraph files
5. Convert bedGraph to bigWig

#### Sort bam file

The output bam files must be sorted by **queryname** in order to generate the BEDPE format in the next step. `"$base"` again refers to the your filename/sample ID: 

`picard SortSam I="$base".bam O="$base"-sorted.bam SO=queryorder CREATE_INDEX=TRUE`

The sorted bam file is converted to bed format using `bedtools bamtobed`. For calibration using the E.coli reads, the bed files require the length of the fragment to be added (as described in the Henikoff lab calibration [script](https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.csh)). Assuming the output file is in the format `"$base".bam` (where "$base" is your sample identifier), the following code may be used:

`bedtools bamtobed -bedpe -i "$base"-sorted.bam | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > "$base".bed`

The E.coli alignment file should also be converted to bed format (no need to sort, since the information used will be the number of reads i.e. the number of lines in the bed file):

`bedtools bamtobed -bedpe -i "$base"-E.coli.bam > "$base"-E.coli.bed`

#### Calibration

If you are working with multiple samples, e.g. a sample and a control, they should be calibrated to be comparable. The calibration is carried out using the Henikoff lab calibration [script](https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.csh)). The calibration effectively scales the mapped counts according to the total number of E.coli reads. The ratio of primary genome to E.coli genome is expected to be the same for all samples. The script is included in this repository (see above).

**Step 1:** Filter fragments to be within a minimum and maximum length

**Steps 2:** Run the calibration on the filtered file

To run the calibration, the query-sorted bed file must be sorted by coordinate:

`bedtools sort -i "$base".bed > "$base"-sorted.bed`

Seven arguments are required to run the calibration script (here converted to bash shell from the Henikoff lab C-shell [script](https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.csh)).

`spike_calibrate.sh genome.bed spike_genome.bed scale output(bg|bga|d) genome_chr_lens_file  min_len max_len`

- **genome.bed** the converted bed file (hg19 alignment), `"$base"-sorted.bed`
- **spike_genome.bed** the converted bed file (E.coli alignment)
- **scale** an arbitrary large number used to scale, e.g. 10000
- **output(bg|bga|d)** bg = BedGraph, bga = BedGraph including regions with 0 coverage, d = depth at each genome position with 1-based coordinates.
- **genome_chr_lens_file**
- **min_len** minumum fragment length, `min=$(cut -f 11 SRR8383480_aligned_reads.bed | sort | uniq | head -1)`
- **max_len** maximum fragment length, `max=$(cut -f 11 SRR8383480_aligned_reads.bed | sort | uniq | tail -1)'

The normalised bedGraphs can now be visualised, for example on the UCSC browser:

![bedGraphUCSC](UCSC-bedgraph.png?raw=true)

#### Convert to bigWig

To convert the bedGraph file to bigWig format, use the following command:

`source bedGraphToBigWig "$base".bedgraph hg19.chrom.sizes "$base".bigWig`

The bedGraphToBigWig binary, as downloaded from UCSC tools, is available in this repository.

## Peak Calling

Peak calling will be carried out using both [macs2](https://github.com/macs3-project/MACS) and SEACR. First, for macs2:

The `"$base"-sorted.bed` file is currently in the bedtools BEDPE format, which is not compatible with macs2. To convert to the macs2 BEDPE format, run the following:

`cut -f 1,2,6 "$base"-sorted.bed | sort -k1,1n -k2,2n -k3,3n > "$base"-sorted-macs.bed` 
 
As described in the original CUT&Tag paper ([Kaya-Okur et al. 2019](https://www.nature.com/articles/s41467-019-09982-5#data-availability)), macs2 will be used with the following parameters:

`macs2 callpeak -t "$base"-sorted-macs.bed -f BEDPE -p 1e-5 --keep-dup all -n output_prefix` 



