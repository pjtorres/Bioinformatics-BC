#!/bin/bash
set -eu -o pipefail

date

#Created: 2020-04-15
# Author: Pedro J. Torres
# Last update: Never

trancript_aligned=""
genome_aligned=""
prefix=""

while  test  $# -gt 0;  do
	case "$1"  in
		-h|--help|-help)
			echo -e "\n\nfilter_genome_STAR.sh\n"
      echo -e " This script filteres out reads based on whether it aligned uniquely to the trancriptome or human genome. "
      echo -e "\n"

			echo -e "\nparameters:"
      echo -e "\t-t trancript\tTrancriptome bam file from STAR"
      echo -e "\t-g genome\tGenome coordinate bam file from STAR"
      echo -e "\t-o prefix\tOutput prefix name from star"

      echo -e "\n\n"

      exit 0
      ;;
		-t)
			shift
			if test $# -gt 0; then
				trancript_aligned="$1"
			else
				echo -e "\n\nERROR! <filename>.starAligned.toTranscriptome.out.bam file not specified.\n\n"
				exit 1
			fi
			shift
			;;
      -g)
  			shift
  			if test $# -gt 0; then
  				genome_aligned="$1"
  			else
  				echo -e "\n\nERROR! <filename>.starAligned.sortedByCoord.out.bam file not specified.\n\n"
  				exit 1
  			fi
  			shift
  			;;
        -o)
          shift
          if test $# -gt 0; then
            prefix="$1"
          else
            echo -e "\n\nERROR! Prefix name used in STAR not specified.\n\n"
            exit 1
          fi
          shift
          ;;
            *)
        			echo -e "\n\nERROR!  nothing requested!\n"
        			exit 1
        			;;
        	esac
done

############################ check parameters


echo -e "\n\n\n"
echo -e "trancript_aligned=[$trancript_aligned]"
echo -e "genome_aligned=[$genome_aligned]"
echo -e "prefix=[$prefix]"

echo -e "\n\n\n"

if [   -z "$trancript_aligned"  -o  -z "$prefix" ]; then
	echo -e "ERROR!  not all inputs specified!\n"
	exit 1
fi

if [ ! -s "$trancript_aligned" ]; then
  echo -e "\n\nERROR!  <filename>.starAligned.toTranscriptome.out.bam does not exist!\n"
  exit 1
fi # 

if [ ! -s "$genome_aligned" ]; then
  echo -e "\n\nERROR!  <filename>.starAligned.sortedByCoord.out.bam does not exist!\n"
  exit 1
fi 


### execute

# Convert bam to sam file while making sure that:
#1. we do not include reads that are unmmaped or mates that are unmmaped (-F 12)
#2. include reads mapped in proper pair (-f 2)
echo -e "Converting bam to sam and inclduing reads mapped in proper pair."
samtools view -F 12 -f 2 $prefix.starAligned.toTranscriptome.out.bam > $prefix.mapped.toTranscriptome.out.sam &&
samtools view -F 12 -f 2 $prefix.starAligned.sortedByCoord.out.bam  > $prefix.mapped.sortedByCoord.out.sam &&

# Convert to bed file for simple fitlering later 
echo -e "Converting sam to bed."
awk '{print $1, 0,100,$0}' FS="\t"  OFS="\t" $prefix.mapped.toTranscriptome.out.sam > $prefix.awk.mapped.toTranscriptome.out.sam.bed &&
awk '{print $1, 0,100,$0}' FS="\t"  OFS="\t" $prefix.mapped.sortedByCoord.out.sam >         $prefix.awk.mapped.sortedByCoord.out.sam.bed &&

# Sort above files 
echo -e "Sorting bed files."
sort -t $'\t' -k1,1  -k2,2n  -k3,3n  $prefix.awk.mapped.toTranscriptome.out.sam.bed  > $prefix.sorted.awk.mapped.toTranscriptome.out.sam.bed &&
sort -t $'\t' -k1,1  -k2,2n  -k3,3n  $prefix.awk.mapped.sortedByCoord.out.sam.bed    > $prefix.sorted.awk.mapped.sortedByCoord.out.sam.bed  &&

#rm $prefix.awk.mapped.toTranscriptome.out.sam.bed &&
#rm $prefix.awk.mapped.sortedByCoord.out.sam.bed &&

# filter out the Coordinate gile to either only include reads that alinged to known human
# trancripts or that are unqiue to only the human genome : 
# -v Only report those entries in A that have no overlap in B
# -wa Write the original entry in A for each overlap
echo -e "Using bedtools."
bedtools intersect -v -a $prefix.sorted.awk.mapped.sortedByCoord.out.sam.bed  -b $prefix.sorted.awk.mapped.toTranscriptome.out.sam.bed > $prefix.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed &&
bedtools intersect -v -sorted -a $prefix.sorted.awk.mapped.sortedByCoord.out.sam.bed -b  $prefix.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed  > $prefix.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed

# covert filtered bed files back to sam
echo -e "Converting bed back to sam file"
awk '{ print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}'  FS="\t"  OFS="\t" $prefix.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed  > $prefix.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam && 
awk '{ print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}'  FS="\t"  OFS="\t" $prefix.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed  > $prefix.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam &&

# get coordinate bam file headers
echo -e "Get bam headers"
samtools view -H $prefix.starAligned.sortedByCoord.out.bam > $prefix.headers.txt

# make new sam file wiht the header
echo -e " Concatinate bam file headers with the new filtered sam file"
cat $prefix.headers.txt $prefix.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam > $prefix.header.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam && 
cat $prefix.headers.txt $prefix.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam  > $prefix.header.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam &&

# convert to bam for later
echo -e " Convert sam to bam "
samtools view -S  -b $prefix.header.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam > $prefix.header.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam &&
samtools view -S  -b $prefix.header.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam  > $prefix.header.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam &&

# mkdir
mkdir -p $prefix.IN.transcript
mkdir -p $prefix.NA.trancript.IN.genome

#sorting and idexing 
echo -e " Sorting and indexing bam file"
samtools sort $prefix.header.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam -o $prefix.sorted.header.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam &&
samtools sort $prefix.header.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam -o $prefix.sorted.header.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam &&
samtools index $prefix.sorted.header.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam &&
samtools index $prefix.sorted.header.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam &&

mv $prefix.sorted.header.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam $prefix.NA.trancript.IN.genome
mv $prefix.sorted.header.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam $prefix.IN.transcript

mv $prefix.sorted.header.NA.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam.bai $prefix.NA.trancript.IN.genome
mv $prefix.sorted.header.IN.trancriptome.sorted.awk.mapped.sortedByCoord.out.sam.bed.sam.bam.bai $prefix.IN.transcript

rm $prefix*bed*
