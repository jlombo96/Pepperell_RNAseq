#!/usr/bin/env bash
# Author Jonathan Lombardino, Burton Lab, University of Wisconsin Madison
# map reads from a folder of RNA-seq BAM files to an input interval file (GFF or BED).
# Usage: bash mapcoverage.bash interval.bed genomefile.genome

intervals=$1
name=${intervals%.*}
output=""$name".coverage.bed"
output_genomecov=""$name".genomecov.bed"

### Loop through every sorted BAM file and count reads that overlap with input interval file (BED/VCF/GFF) ### 
for i in *.bam;
do 
	bam_name=${i%.*} #get the name of the file
	bedtools coverage -a $intervals -b $i -s -bed | sed 's/$/\t'$bam_name'/'  >> $output #count overlap of reads in bam file and add column denoting sample, require strandedness, split reads are two intervals
	bedtools genomecov -ibam $i |  sed 's/$/\t'$bam_name'/' >> $output_genomecov #count read-depth at a per base level. 
done






