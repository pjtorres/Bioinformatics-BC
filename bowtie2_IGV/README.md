# Looking at genomic polymorphisms

## Building bowtie index

`bowtie2-build <genome>.fa index/<genome>`
## always go to inspect and make sure everything was built properly

`bowtie2-inspect -s index/<genome>`

# bowtie2 --end-to-end 
## Bowtie -  for fastq
`bowtie2 --end-to-end -x index/<genome> -1 <read_R1>.fq -2 <read_R2>.fq  -S <read>.sam`

## for fna
`bowtie2 --end-to-end   -x index/<genome>  -f  <read>.fna  -S <read>.sam`

# bowtie2 --local 

`bowtie2 --local -N 1 -x index/<genome> -f <read>.fna -S <read>.sam` 

# Now finish up and prepping for IGV

## Convert to bam
`samtools view -S -b <read>.sam > <read>.bam`

## sort samtools
`samtools sort <read>.bam >  <read>.sorted.bam`

## full script
`python pull_COVID19_reads_to_bowtie.py -l bioinfo-prod/centrifuge_taxa_output -s 1B9122304251 -i index/COVID_refseq`

## index 
`samtools index <read>.bam`

You will now have two important files needed for IGV that is the <read>.sorted.bam and <read>.sorted.bam.bai

## igv
now you can always download the fna whole genome file and then load to IGV by 1. click on load genome, 2. load genome from file and click on the genome. Then you can download the associated gff file and load that by clicking on file and then load from file. This will result in the genome size on top and the ORFs on the bottom. Remember to use the same genome used for bowtie2 aligning as what is in the gff file.  Once you have loaded your reference genome, click on the file tab on top of IV and load from file. Go to your folder which should have your <read>.bam and <read>.bam.bai files in it and click on the bam. You should now be able to visualize your read coverage. 
  
## What is you would like to get a file showing the bp and coverage at that position and then graph it yourself. This is handy when you have a lot of reads and want to quantify the coverage a bit better. 

## bedtools is a great tool
` bedtools genomecov -ibam <read>.sorted.bam -d | awk '$3>1' > <read>.sorted.bam.genomeCov.txt`

You will get a text file with the name of each alignment query and the bp position as well as how many times you had a bp align to it. You can now use this text file and put thorugh R or python and make your own coverage plots focusing on a gene of interest.
