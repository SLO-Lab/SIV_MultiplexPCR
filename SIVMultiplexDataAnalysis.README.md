# SIV Multiplex PCR Sequencing May 2019

## Introduction

In order to normalize the sequencing data, we modified the Zequencer pipeline from the DHO lab to work with our samples. Pipeline available at: https://bitbucket.org/dholab/
This tool requires a variety of scripts to work properly. From the Zequencer2017 README -  

	"We have several studies where reproducible viral sequencing is important. Last year I put together a kludgy read mapping and variant calling tool, called Zequencer, that I implemented in Geneious Pro. This is convenient for users not comfortable in the command line. But it is difficult to run in a reproducible way. Zequencer 2017 is a command-line tools that masks most of the complexity of viral sequencing from users while nonetheless generating results that can be reproduced exactly. In the initial implementation, Zequencer 2017 contains three separate but related tools:

	1. prepareNcbiReference - takes an NCBI accession number and creates required indices for novoalign, samtools, and snpEff." 

We then modified two other scripts in the Zequencer for our use: 

2. normalizeShortAmpliconCoverage_qtrim - given a group of amplicon sequences from a sample, normalizes to a fixed number of reads per amplicon to create equal coverage throughout the genome. We added an additional quality trim step after merging of R1 and R2 FASTQ files.

3. mapReadsCallVariants_SIV - maps reads to a reference genome with Novoalign, calls intrahost variants with VarScan, and annotates the impact of these variants using snpEff. Minimum variant threshold is 1%, reads are mapped to all possible locations to account for LTR mapping. 

Copied from the Zequencer 2017 README:
	## Dependencies
	
	The Zequencer folder contains binaries for necessary processing applications compiled for OSX 10.12:
	
	+ bbmap_36.86
	+ novoalign
	+ novoindex
	+ samtools
	+ seqtk
	+ snpEff
	+ VarScan 2.4.3
	
	Zequencer is written in Python3 and assumes all modules present in the Anaconda Python 3.5 distribution are present. Download the OSX version of Anaconda from here and install with the graphical installer:
	
	```
	https://repo.continuum.io/archive/Anaconda3-4.2.0-MacOSX-x86_64.pkg
	```
	
	Several Python modules not included in the Anaconda distribution are needed for Zequencer 2017. After Anaconda is installed in a user home directory, install these modules with:
	
	```
	~/anaconda/bin/pip3.5 install pyfasta biopython
	```
	
	## prepare NCBI reference
	
	This script accepts an NCBI accession number as a parameter and creates the necessary indices for the mapReadsCallVariants script to operate correctly. Reference indices must be downloaded with this script before they can be used as a reference in mapReadsCallVariants. For example, to download and index the ZIKV French Polynesian strain with an accession number of KJ776791, run this command:
	
	```
	~/anaconda/bin/python3.5 \
	/Users/dho/bitbucket/xdhofs2/18582/Zequencer/prepareNcbiReference.py \
	'KJ776791'
	```  
	
	For those unaccustomed to running Python scripts, the first line is the location of the Python interpreter. The second line is the path to the file containing the prepareNcbiReference.py script that is run. Shown here is the path to the script on my computer. Your path will look different depending on where Zequencer is stored on your computer. The third line is the accession number of the reference genome to index.
	
	After you run the script, you should see output that looks like:
	
	```
	Saved
	KJ776791
	# novoindex (3.6) - Universal k-mer index constructor.
	# (C) 2008 - 2011 NovoCraft Technologies Sdn Bhd
	# novoindex /Users/dho/bitbucket/xdhofs2/18582/Zequencer/ref/KJ776791/KJ776791.fa.nix /Users/dho/bitbucket/xdhofs2/18582/Zequencer/ref/KJ776791/KJ776791.fa 
	# Creating 4 indexing threads.
	# Building with 9-mer and step of 1 bp.
	# novoindex construction dT = 0.0s
	# Index memory size   0.001Gbyte.
	# Done.
		Protein check:	KJ776791	OK: 1	Not found: 0	Errors: 0	Error percentage: 0.0%
	
	```
	
	In your Zequencer 2017 folder, there should now be a folder within the 'ref' directory that contains five files. You can check that these files are all present in the right directory by listing the contents of the directory using the UNIX ls command like this:
	
	```
	ls /Users/dho/bitbucket/xdhofs2/18582/Zequencer/ref/KJ776791
	```
	
	Once again, the path will vary depending on where Zequencer is installed on your computer. The results of the ls statement should look something like this:
	
	```
	KJ776791.fa		KJ776791.fa.nix		snpEffectPredictor.bin
	KJ776791.fa.fai		KJ776791.gbk
	```
	
	The .fa file is a FASTA formatted version of the genome.
	The .fai file is an indexed version of the FASTA file.
	The .fa.nix file is the index needed by Novoalign.
	The .gbk file is a Genbank file that includes the protein annotations that are needed by snpEff.
	The snpEffPredictor.bin file is snpEff's index that is required for it to determine the impact of variants in mapped reads.
	
	## normalize short amplicon coverage
	
	This script should be used with viral sequencing data collected using the highly multiplexed PCR methodology described in http://biorxiv.org/content/early/2017/01/09/098913. Briefly, two pools of non-overlapping short PCR amplicons tiling the genome are prepared and sequenced simultaneously. For ZIKV, this can mean dividing ~35 amplicons into two individual pools for PCR. Sequence coverage from each amplicon is inconsistent. Some amplicons will be represented by many more reads than others. Since the number of sequences in an Illumina miSeq experiment is typically much higher than the number of templates in the input sample, reads from each amplicon can be downsampled to a fixed number of reads without any effective loss in ability to resolve variants. Operationally, downsampling to 1000 reads per amplicon will likely provide high confidence data for detecting variants that are present in at least 5% of a sample. It is very important to manually inspect variant calls, however, since some variants detected very near the end of an amplicon may be artifacts.
	
	The script performs a number of steps:
	+ runs bbmap's bbmerge tool to create a single FASTQ read from each pair of R1 and R2 reads, joining the reads in the overlapping region in the middle to maximize quality
	+ run bbmap using each of the short amplicon sequences as a reference. Temporarily save the FASTQ reads that map to each amplicon
	+ run seqtk to downsample the mapped FASTQ reads from each amplicon to the user-specified number of reads
	+ concatenate downsampled, mapped reads to a single output file where each amplicon is represented by approximately the same number of sequences
	
	The script requires several parameters to run successfully:
	
	+ Path to a FASTA file containing predicted sequences for each of the short amplicons. These will be used individually as reference sequences for read mapping. Of the reads that map to each amplicon, a user-selected number will be retained in the downsampled output FASTQ file.
	+ Path to FASTQ file (R1) from paired-read Illumina miSeq dataset. Should be GZIP compressed (have a .gz extension)
	+ Path to FASTQ file (R2) from paired-read Illumina miSeq dataset. Should be GZIP compressed (have a .gz extension)
	+ Number of reads to retain from each amplicon after downsampling
	
	The syntax for running the script is:
	
	```
	~/anaconda/bin/python3.5 \
	~/PycharmProjects/xdhofs2/18582/Zequencer/normalizeShortAmpliconCoverage.py \
	'/Network/Servers/xdhofs2.pathology.wisc.edu/Volumes/odin/netUsers-xdhofs2/dho/PycharmProjects/xdhofs2/18582/Zequencer/ref/short_amplicon_ref/Zika-PR-amplicons.fasta' \
	'/Volumes/odin/netUsers-xdhofs2/dho/Desktop/18582-ensembl-vep/fastq/Original-FASTQ-reads/replicate-A/PR-ABC59-Harvestvirus-V3C2-replicateA_S21_L001_R1_001.fastq.gz' \
	'/Volumes/odin/netUsers-xdhofs2/dho/Desktop/18582-ensembl-vep/fastq/Original-FASTQ-reads/replicate-A/PR-ABC59-Harvestvirus-V3C2-replicateA_S21_L001_R2_001.fastq.gz' \
	'1000'
	```
	
	The first line is the path to the python interpreter.
	The second line is the path to the script.
	The third line is the path to the FASTA file with predicted sequences for each of the short amplicons.
	The fourth line is the R1 FASTQ file.
	The fifth line is the R2 FASTQ file.
	The sixth line is the number of reads to downsample out of each amplicon.
	
	There are two output files:
	+ merged FASTQ file from the whole R1 and R2 FASTQ datasets
	+ downsampled FASTQ file with approximately the same number of reads from each amplicon
	
	## map reads and call variants
	
	This script provides a reproducible way to map reads, call variants, and determine the functional impact of variants on open reading frames. Defaults are set for several important parameters and exact commands are logged for future reproducibility. 
	
	The script performs a number of steps:
	+ uses Novoalign (which is the recommended read mapper from Kristian Andersen) to map reads to a reference sequence, outputting a SAM file
	+ samtools sorts the SAM file and converts to BAM format
	+ samtools mpileup formatted-variants generated
	+ variants identified with VarScan
	+ variants analyzed for impact on the open reading frame using snpEff
	
	The minimum syntax for running the script is:
	```
	~/anaconda/bin/python3.5 \
	~/PycharmProjects/xdhofs2/18582/Zequencer/mapReadsCallVariants.py \
	'KU501215' \
	--fastq '/Volumes/odin/netUsers-xdhofs2/dho/Desktop/18582-ensembl-vep/fastq/Original-FASTQ-reads/replicate-A/PR-ABC59-Harvestvirus-V3C2-replicateA_S21_L001_R1_001.1000.reads.per.amplicon.fastq.gz'
	```
	
	First line is the path the python executable. Second line is the path to the script. Third line is the accession number of the reference sequence. Remember that this reference sequence needs to be indexed with the prepareNcbiReference script before it can be used in this script. The fourth line is the path to the FASTQ reads to map. Note that you **need** to include the:
	
	```
	--fastq
	```
	
	prefix before the path to the FASTQ file. What is you have two FASTQ files? No problem. Your command line should look like this:
	
	```
	~/anaconda/bin/python3.5 \
	~/PycharmProjects/xdhofs2/18582/Zequencer/mapReadsCallVariants.py \
	'KU501215' \
	--fastq '/Volumes/odin/netUsers-xdhofs2/dho/Desktop/18582-ensembl-vep/fastq/Original-FASTQ-reads/replicate-A/PR-ABC59-Harvestvirus-V3C2-replicateA_S21_L001_R1_001.fastq.gz' \
	'/Volumes/odin/netUsers-xdhofs2/dho/Desktop/18582-ensembl-vep/fastq/Original-FASTQ-reads/replicate-A/PR-ABC59-Harvestvirus-V3C2-replicateA_S21_L001_R2_001.fastq.gz'
	```
	
	Note how there are just two FASTQ files separated by a space. 
	
	You can also add additional options to the command line to modify the output, though the default settings should be sufficient for most analysis:
	
	```
	--number_reads [default = all]
	```
	
	Useful for datasets that were not previously downsampled. Can be used to limit the number of reads analyzed from the FASTQ file. 
	
	```
	--min_var_percent [default = '0.05']
	```
	
	When running VarScan, only report variants found in more than this number of reads. Default is set to 5%.
	
	```
	--qtrim [default = 't']
	--min_read_length [default = '100']
	```
	
	Low quality reads are removed using bbmap's quality filter using its default settings. To disable quality trimming, set to 'f'. Reads shorter than a specified length can also be removed prior to mapping to the reference genome. This prevents short reads from poisoning mappings.
	
	```
	--debug [default = 'f']
	```
	
	Debug mode saves intermediate files that can be useful for diagnosing problems with the script. These files are typically deleted at the end of each run to save disk space. Recall that the files can always be regenerated if needed.
	
	A command using these advanced parameters would look like this:
	
	```
	~/anaconda/bin/python3.5 \
	~/PycharmProjects/xdhofs2/18582/Zequencer/mapReadsCallVariants.py \
	'KU501215' \
	--fastq '/Volumes/odin/netUsers-xdhofs2/dho/Desktop/18582-ensembl-vep/fastq/Original-FASTQ-reads/replicate-A/PR-ABC59-Harvestvirus-V3C2-replicateA_S21_L001_R1_001.fastq.gz' \
	'/Volumes/odin/netUsers-xdhofs2/dho/Desktop/18582-ensembl-vep/fastq/Original-FASTQ-reads/replicate-A/PR-ABC59-Harvestvirus-V3C2-replicateA_S21_L001_R2_001.fastq.gz' \
	--number_reads '500' \
	--min_var_percent '0.10' \
	--qtrim 'f' \
	--min_read_length '150' \
	--debug 't'
	```
	
	This script generates three files when run in normal mode (debugging turned off):
	+ sorted BAM file
	+ impact-annotated VCF file
	+ FASTA reference file
	
	These files can be visualized in Geneious Pro or other software. In Geneious, load the FASTA file first, then the VCF file, and finally the BAM file. 
	
	## troubleshooting
	
	To test the portability of this code, I asked Shelby to transfer it to her computer and follow these instructions. There were a few issues:
	
	#### Unsupported major.minor version java error
	snpEff requires a recent version of Java. OSX 10.12 has an acceptable version installed by default. Shelby is running OSX 10.10, which doesn't. As of now, download and install the java JDK from:
	
	http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
	
When analyzing the biological data requested by reviewers, we used an updated pipeline designed more specifically to work with biological samples, available at https://github.com/gagekmoreno/SARS_CoV-2_Zequencer, version 2, but with the NCBI reference M33262, 2000000 reads per sample, 0.01 minimum variant percentage, and minimum read length 100. 
	
