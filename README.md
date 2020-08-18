# CUT&Tag
Step-by-step analysis pipeline for CUT&Tag data
#### [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)

The CUT&Tag protocol is available as the [CUT&Tag@home](https://www.protocols.io/view/cut-amp-tag-home-bd26i8he?step=50) and [Bench top CUT&Tag V.3](https://www.protocols.io/view/bench-top-cut-amp-tag-bcuhiwt6). The following pipeline is adapted from the authors recommended analysis protocol, available [here](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-bjk2kkye). 

The following pipeline describes each analysis step:

- [Pre-alignment quality scontrol (QC)](pre-alignment-qc) 
- [Alignment](#alignment)
- [Post-alignment QC](#post-alignment-qc)
- [Visualisation & Calibration](#visualisation) 
- [Peak Calling](#peak-calling)
- [Motif Finding](#motif-finding)
- Cut Matrix Generation
- Motif Footprinting Steps

All required programs required have been installed and are available in the CebolaLab [CUTandTAG anaconda environment](https://github.com/CebolaLab/CUTandTAG/tree/master/anaconda-env). If you are using anaconda, you can copy this into your own anaconda environments (on the Imperial College HPC, for example, this is located at `/rdsgpfs/general/user/"$(whoami)"/home/anaconda3/envs/`) and then `source activate CUTandTAG`. The pipeline can also be run without installing anaconda by directing the necessary scripts to the `bin` of the downloaded CUTandTAG environment (download the CUTandTAG directory and save it e.g. to your home directory).

For the following analysis, you can save your sample file name as `base` and the example scripts will access this variable using `<sample>`. 

## Pre-alignment QC

A common tool used to assess the quality of raw sequence reads is [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). The report, which comes in html format, can be used as guidance to the general quality of the data. Some failed or 'warning' results are usually not too concerning and can be improved using various processing tools. One metric which may be observed to fail is the Per base sequence content. As [described by the authors](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-bjk2kkye?step=6) this can result from a sequence preference of the Tn5 transposase or the 10bp periodicity in the length distribution:

![FastQC](Figures/fastQC-per-base-quality.png)

If using the [Benchop CUT&Tag protocol](https://www.protocols.io/view/bench-top-cut-amp-tag-bcuhiwt6/abstract), sequence length should be 25bp and therefore there is not expected to be any contamination of adapters (which results when the sequence lenght is longer than the DNA fragment and so the sequencing extends into the adapter sequence). Therefore, adapter trimming, which is a common step in NGS pipelines, is not recommended unless the user has opted to sequence longer reads (>25bp).

## Alignment

Two alignments will be run to align the human DNA and carry-over E.coli DNA later used to standardise the samples. The alignment parameters are run as recommended by the CUT&Tag authors (see their [pipeline](https://www.protocols.io/view/cut-amp-tag-home-bd26i8he?step=50)). The authors recommend to skip adapter trimming and to run the alignments using bowtie2 with the below parameters, which should result in accurate read alignment. Two alignments are carried out:

1. Align reads to the reference **human** masked genome (hg19) (download [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/))
2. Align reads to the reference **E.coli** genome (strain K12, substrain MG1655) (downloaded [here](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3?report=fasta)).

Both reference genomes should be [indexed](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome) using bowtie2. For those with access to the Imperial College HPC and the Cebola Lab project space, the reference genomes and index files are available at this path:

`/rds/general/user/"$(whoami)"/projects/cebolalab_liver_regulomes/live/reference-genomes/` 

The alignments are carried out using bowtie2 with the below arguments. An example script is available [here](https://github.com/CebolaLab/CUTandTAG/blob/master/alignment.sh).

##### Parameters to align human reads:

`--end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -l 10 -X 700`

##### Parameters to align E.coli reads:

`--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant  --phred33 -I 10 -X 700`

## Post-alignment QC

- Remove mitochondrial reads
- Remove duplicate reads
- Remove reads with a mapping quality <30 (includes non-uniquely mapped reads)

### Remove mitochondrial reads

To assess the total % of mitochondrial reads, run `samtools idxstats` and `samtools flagstats` on your indexed bam file. For example:

`samtools index <sample>.bam`

`samtools flagstat <sample>.bam > <sample>.flagstat`

`samtools idxstats <sample>.bam > <sample>.idxstats`

The total number of reads can be seen the the first line of the flagstat file, `head -n1 <sample>.flagstat`, while the total number of mitochondrial reads can be seen using `grep chrM <sample>.idxstats` (the 2nd column is the length of the chromosome and the 3rd column is the total number of aligned reads). Flagstat and idxstats files can be generated after each QC step to assess the total number of reads removed.  Mitochondrial reads should be removed as follows:

`samtools view -h <sample>.bam | grep -v chrM | samtools sort -O bam -o <sample>.rmChrM.bam -T .`

### Mark and remove duplicate reads

`picard MarkDuplicates QUIET=true INPUT=<sample>.rmChrM.bam OUTPUT=<sample>.marked.bam METRICS_FILE=<sample>.sorted.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.`

`samtools view -h -b -F 1024 <sample>.markerd.bam > <sample>.rmDup.bam`

### Remove mapping quality <30

Reads which are uniquely mapped are assigned a high alignment quality score and one genomic position. If reads can map to more than one location, Bowtie2 reports one position and assigns a low quality score. The proportion of uniquely mapped reads can be assessed. In general, >70% uniquely mapped reads is expected, while <50% may be a cause for concern 
[(Bailey et al. 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/pdf/pcbi.1003326.pdf). A low % of uniquely mapped reads may, however, be observed when carrying out a CUT&Tag experiment for a protein which is expected to bind repetitive DNA. Alternatively, this may result from short reads, excessive PCR amplification or problems with the PCR [(Bailey et al. 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/pdf/pcbi.1003326.pdf).

`samtools view -q 10 -b <sample>.rmDup.bam > <sample>.filtered.bam`

### QC report

A hmtl QC report can be output using [qualimap](http://qualimap.bioinfo.cipf.es/). 

`qualimap bamqc -bam <sample>.filtered.bam -gd hg19 -outdir . -outfile <sample>.qualimap.report -outformat html`

## Visualisation

The aligned data (output from the example [script](https://github.com/CebolaLab/CUTandTAG/blob/master/alignment.sh)) will be in `bam` format. This will be converted to bedGraph format in order to visualise the data. For multiple samples to be compared, the samples will first be calibrated using the carry-over E.coli DNA (the effective 'spike-in'). In theory, the ratio of primary DNA to E.coli DNA is expected to be the same for each sample. As such, the calibration divides the mapped read count by the total number of reads aligned to the E.coli genome. The proportion of total DNA reads aligned to the E.coli genome is reported in the `<sample>.Ecoli.bowtie2` output file. The general steps cover:

1. Sort aligned bam file by read name (queryname)
2. Convert bam to bed
3. Normalise samples using E.coli carry-over DNA
4. Visualise output bedGraph files
5. Convert bedGraph to bigWig

### Sort bam file

The output bam files must be sorted by **queryname** in order to generate the BEDPE format in the next step. `<sample>` again refers to the your filename/sample ID: 

`picard SortSam I=<sample>.bam O=<sample>-sorted.bam SO=queryname CREATE_INDEX=TRUE`

The sorted bam file is converted to bed format using `bedtools bamtobed`. For calibration using the E.coli reads, the bed files require the length of the fragment to be added (as described in the Henikoff lab calibration [script](https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.csh)).

`bedtools bamtobed -bedpe -i <sample>-sorted.bam | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > <sample>.bed`

The E.coli alignment file should also be converted to bed format (no need to sort, since the information used will be the number of reads i.e. the number of lines in the bed file):

`bedtools bamtobed -bedpe -i <sample>-E.coli.bam > <sample>-E.coli.bed`

### Calibration

If you are working with multiple samples, e.g. a sample and a control, they should be standardized in order to be comparable. The calibration is carried out using the Henikoff lab calibration [script](https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.csh). The calibration effectively scales the mapped counts according to the total number of E.coli reads. The ratio of primary genome to E.coli genome is expected to be the same for all samples. The script is included in this repository (see above).

**Step 1:** Filter fragments to be within a minimum and maximum length

**Steps 2:** Run the calibration on the filtered file

To run the calibration, the query-sorted bed file must be sorted by coordinate:

`bedtools sort -i <sample>.bed > <sample>-sorted.bed`

Seven arguments are required to run the calibration script (here converted to bash shell from the Henikoff lab C-shell [script](https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.csh)).

`spike_calibrate.sh genome.bed spike_genome.bed scale output(bg|bga|d) genome_chr_lens_file  min_len max_len`

- **genome.bed** the converted bed file (hg19 alignment), `<sample>-sorted.bed`
- **spike_genome.bed** the converted bed file (E.coli alignment)
- **scale** an arbitrary large number used to scale, e.g. 10000
- **output(bg|bga|d)** bg = BedGraph, bga = BedGraph including regions with 0 coverage, d = depth at each genome position with 1-based coordinates.
- **genome_chr_lens_file**
- **min_len** minumum fragment length, `min=$(cut -f 11 SRR8383480_aligned_reads.bed | sort | uniq | head -1)`
- **max_len** maximum fragment length, `max=$(cut -f 11 SRR8383480_aligned_reads.bed | sort | uniq | tail -1)'

### Visualise bedGraph

The normalised bedGraphs can now be visualised, for example on the UCSC browser:

![bedGraphUCSC](UCSC-bedgraph.PNG?raw=TRUE)

### Visualise heatplots and profiles

You can also generate a heatplot to visualise the distribution of your chromatin mark / transcription factor relative to transcription start sites. This will use deeptools (included in the CUTandTAG conda bin). Gene coordinates for the reference genome hg19 were downloaded from UCSC as Gencode V34lift37 (Basic table and bed format). They are saved in this repository as `hg19-gene-coordinates.bed`.

`computeMatrix scale-regions -S <sample>.bigWig -R hg19-gene-coordinates.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --missingDataAsZero --skipZeros -o matrix.mat.gz`

`plotHeatmap -m matrix.mat.gz  -out ExampleHeatmap1.png ` 

The parameters of `plotHeatmap` can be [adjusted](https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html) to produce various heatplots: 

`plotProfile -m matrix.mat.gz -out ExampleProfile1.png `

<img src="https://github.com/CebolaLab/CUTandTAG/blob/master/Figures/ExampleHeatmap1.png" width="200"><img src="https://github.com/CebolaLab/CUTandTAG/blob/master/Figures/ExampleHeatmap2.png" width="200"><img src="https://github.com/CebolaLab/CUTandTAG/blob/master/Figures/ExampleHeatmap5.png" width="200"><img src="https://github.com/CebolaLab/CUTandTAG/blob/master/Figures/ExampleHeatmap4.png" width="200">

`deeptools` can also be used to plot profiles. Different [parameters](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html) can again be used to achieve different aesthetics:

<img src="https://github.com/CebolaLab/CUTandTAG/blob/master/Figures/ExampleProfile1.png" width="200"><img src="https://github.com/CebolaLab/CUTandTAG/blob/master/Figures/ExampleProfile2.png" width="200"><img src="https://github.com/CebolaLab/CUTandTAG/blob/master/Figures/ExampleProfile3.png" width="200">

### Convert to bigWig

To convert the bedGraph file to bigWig format, use the following command:

`source bedGraphToBigWig <sample>.bedgraph hg19.chrom.sizes <sample>.bigWig`

The bedGraphToBigWig binary, as downloaded from UCSC tools, is available in this repository.

## Peak Calling

Peak calling will be carried out using both [macs2](https://github.com/macs3-project/MACS) and SEACR. First, for macs2:

The `<sample>-sorted.bed` file is currently in the bedtools BEDPE format, which is not compatible with macs2. To convert to the macs2 BEDPE format, run the following:

`cut -f 1,2,6 <sample>-sorted.bed | sort -k1,1n -k2,2n -k3,3n > <sample>-sorted-macs.bed` 
 
As described in the original CUT&Tag paper ([Kaya-Okur et al. 2019](https://www.nature.com/articles/s41467-019-09982-5#data-availability)), macs2 will be used with the following parameters:

`macs2 callpeak -t <sample>-sorted-macs.bed -f BEDPE -p 1e-5 --keep-dup all -n output_prefix` 

## Motif Finding



