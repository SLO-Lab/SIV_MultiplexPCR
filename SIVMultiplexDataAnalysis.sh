#! /bin/bash

# Shell Script that calls many other scripts for data analysis of SIV Multiplex Sequencing 
# insert your file names in quotes with a space between them. Make sure the -R1/-R2 and whatever sample ID is not part of the name! Make sure your files are unzipped. 
# It was recommended to us to run the even and odd amplicons separately and then combine the data, but they can be run together. 

declare -a files=( "0-m239-1e3")

#here it will loop through everything
for i in "${files[@]}"
do 
	#This tells you where you are in the process 
	echo $i

	~/anaconda/bin/python3.5 \
	/Volumes/DikDik/SIVmultiplexPaper_MEs/Zequencer_SIV/normalizeShortAmpliconCoverage_qtrim.py \
	'/Volumes/DikDik/SIVmultiplexPaper_MEs/Zequencer_SIV/ref/short_amplicon_ref/SIVmac239-amplicons_evens.fasta' \
	'/Volumes/DikDik/SIVmultiplexPaper_MEs/ME6_evenamplicons/'$i'-R1.fastq' \
	'/Volumes/DikDik/SIVmultiplexPaper_MEs/ME6_evenamplicons/'$i'-R2.fastq' \
	'2000'
	
	#this wait tells the computer not to start doing this top section with the next sample before it completely finishes with the first 
	wait 
	
	#commands for mapping
	~/anaconda/bin/python3.5 \
	/Volumes/DikDik/SIVmultiplexPaper_MEs/Zequencer_SIV/mapReadsCallVariants_SIV.py  \
	'M33262' \
	--fastq '/Volumes/DikDik/SIVmultiplexPaper_MEs/ME6_evenamplicons/'$i'-R1.2000.reads.per.amplicon.fastq.gz'

	
done

wait 

# here we make a CSV file of all the variants called in the VCF file

Rscript /Volumes/DikDik/SIVmultiplexPaper_MEs/Zequencer_SIV/VCFtoCSV.R ME6_evenamplicons /Volumes/DikDik/SIVmultiplexPaper_MEs/ ME6evenampliconfreqs.csv

# here we make a CSV file of all the nucleotide counts for each position in the genome 

Rscript /Volumes/DikDik/SIVmultiplexPaper_MEs/Zequencer_SIV/NucleotideCountTables.R  /Volumes/DikDik/SIVmultiplexPaper_MEs/ME6_evenamplicons