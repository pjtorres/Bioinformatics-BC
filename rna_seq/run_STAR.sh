#!/bin/bash
set -eu -o pipefail

date

# Created: 2020-03-13
# Author: Pedro J. Torres
# Last update: Never

infq1=""
infq2=""

outpref=""
starindex=""
numthreads="1"
gtf=""

while  test  $# -gt 0;  do
	case "$1"  in
		-h|--help|-help)
			echo -e "\n\nrun_STAR.sh\n"
      echo -e " This script runs STAR algorithm. Make sure it is set in your environment and have already made an index. "
      echo -e "\n"

			echo -e "\nparameters:"
      echo -e "\t-x ix\tdirectory of star index"
      echo -e "\t-g gtf\tfilename of star gtf file"
      echo -e "\t-1 fq1\tfilename of fastq file for mate 1 (for paired-end reads)."
			echo -e "\t-2 fq2\tfilename of fastq file for mate 2 (for paired-end reads)."
      echo -e "\t-o prefix\tprefix name of output files"

      echo -e "\t-p\tnumber of threads for multithreading [default: 1]."

      echo -e "\n\n"

      exit 0
      ;;
		-1)
			shift
			if test $# -gt 0; then
				infq1="$1"
			else
				echo -e "\n\nERROR! fq1 file not specified.\n\n"
				exit 1
			fi
			shift
			;;
      -2)
  			shift
  			if test $# -gt 0; then
  				infq2="$1"
  			else
  				echo -e "\n\nERROR! fq2 file not specified.\n\n"
  				exit 1
  			fi
  			shift
  			;;
  		-x)
  			shift
  			if test $# -gt 0; then
  				starindex="$1"
  			else
  				echo -e "\n\nERROR! STAR index directory location not specified.\n\n"
  				exit 1
  			fi
  			shift
  			;;
      -o)
    			shift
    			if test $# -gt 0; then
    				outpref="$1"
    			else
    				echo -e "\n\nERROR! output prefix not specified.\n\n"
    				exit 1
    			fi
    			shift
    			;;
      -p)
      	   shift
      			if test $# -gt 0; then
      				numthreads="$1"
      			else
      				echo -e "\n\nERROR! number of threads not specified.\n\n"
      				exit 1
      			fi
      			shift
      			;;
        -g)
        			shift
        			if test $# -gt 0; then
        				gtf="$1"
        			else
        				echo -e "\n\nERROR!gtf filename not specified.\n\n"
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
echo -e "infq1=[$infq1]"
echo -e "infq2=[$infq2]"
echo -e "starindex=[$starindex]"
echo -e "numthreads=[$numthreads]"
echo -e "gtf=[$gtf]"
echo -e "outpref=[$outpref]"

echo -e "\n\n\n"

if [   -z "$numthreads"  -o  -z "$starindex" ]; then
	echo -e "ERROR!  not all inputs specified!\n"
	exit 1
fi

if [ ! -s "$infq1" ]; then
  echo -e "\n\nERROR!  Foward read does not exist!\n"
  exit 1
fi # 

if [ ! -s "$infq2" ]; then
  echo -e "\n\nERROR!  Reverse read does not exist!\n"
  exit 1
fi 

if [ ! -s "$starindex" ]; then
  echo -e "\n\nERROR!  starindex does not exist!\n"
  exit 1
fi

if [ ! -s "$gtf" ]; then
  echo -e "\n\nERROR!  gtf does not exist!\n"
  exit 1
fi

if [ ! -z "$infq1"  -a  -z "$outpref" ];  then
  echo -e "\n\nERROR!  must provide output prefix name does not exist!\n"
  exit 1
fi

### prepare input arguemnts
echo -e "prepare input arguments."

inputreadarg=""

inputreadarg=" $infq1 $infq2  "

### execute


STAR                              \
--genomeDir $starindex            \
--runThreadN $numthreads          \
--runMode alignReads               \
--readFilesIn $inputreadarg         \
--outFileNamePrefix $outpref.star   \
--outSAMtype BAM SortedByCoordinate  \
 --outFilterType BySJout             \
 --outFilterMultimapNmax 20          \
 --outFilterMismatchNmax 999          \
 --outFilterMismatchNoverReadLmax 0.04 \
 --alignIntronMin 20                   \
 --alignIntronMax 1000000              \
 --alignMatesGapMax 1000000             \
 --alignSJoverhangMin 8                   \
 --alignSJDBoverhangMin 1                 \
 --outSAMprimaryFlag OneBestScore         \
 --outSAMattributes All                  \
 --outSAMattrIHstart 0                    \
 --quantMode TranscriptomeSAM GeneCounts  \
 --outReadsUnmapped Fastx                   \
--sjdbGTFfile $gtf                          \
--outSAMunmapped Within
