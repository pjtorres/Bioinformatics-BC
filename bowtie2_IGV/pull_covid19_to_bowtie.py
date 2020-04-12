__author__= 'Pedro J. Torres'
import subprocess
import os,argparse

#--------Command line arguments------------------------------------------------------------------------------------
parser=argparse.ArgumentParser(description="Script allows you to generate insilico reads and add to wanted sample ")
parser.add_argument('-s','--sample', help='Input sampleectory, sample sampleectory',required=True)
parser.add_argument('-l','--loc', help='Location of sample in s3',required=True)
parser.add_argument('-i','--index', help='Location of sample in s3',required=True)


args = parser.parse_args()

sample=str(args.sample)
loc=str(args.loc)
index=str(args.index)

#aws s3 sync
aws = 'aws s3 sync s3://'+loc+'/'+sample +' '+sample+'/'
subprocess.call(aws, shell=True)

# find all reads that are alinnign to COVID19 and that are uniquel

collect = 'zcat '+sample+'/'+sample+'.centrifuge_taxa.classification.microbiome.tsv.gz | grep "NC_045512.2" | cut -f1 > ' +sample+'/'+sample+'.COVID19reads.txt'
print(collect)
subprocess.call(collect, shell=True)

# pull only reads that aling to COVID19 into their own fastq files
pull = 'filterbyname.sh in='+sample+'/'+sample+'.microbiome.pass.1.fq.gz in2='+sample+'/'+sample+'.microbiome.pass.2.fq.gz names='+sample+'/'+sample+'.COVID19reads.txt '+'include=t out='+sample+'/'+sample+'.COVID19.R1.fq.gz out2='+sample+'/'+sample+'.COVID19.R2.fq.gz'
print(pull)
subprocess.call(pull, shell=True)

#merge these reads
merge = 'bbmerge.sh in1='+sample+'/'+sample+'.COVID19.R1.fq.gz in2='+sample+'/'+sample+'.COVID19.R2.fq.gz out='+sample+'/'+sample+'.COVID19.PE.merged.fq.gz'
print(merge)
subprocess.call(merge, shell=True)

# not everything merges so I will add the forward reads sowill go back and pull foward reads that did not mapped
merged= "zcat "+sample+"/"+sample+".COVID19.PE.merged.fq.gz | grep '^@'| sed 's/@//g' > "+ sample+"/"+sample+".COVID19.PE.merged.reads.txt"
print(merged)
subprocess.call(merged, shell=True)

#now pull from foward reads
pull2 = 'filterbyname.sh in='+sample+'/'+sample+'.COVID19.R1.fq.gz names='+sample+'/'+sample+'COVID19.PE.merged.reads.txt include=f out='+sample+'/'+sample+'.COVID19.R1.notmerged.fq.gz'
print(pull2)
subprocess.call(pull2, shell=True)

#merge all 
zcat = 'zcat '+sample+"/"+sample+".COVID19.PE.merged.fq.gz "+ sample+"/"+sample+".COVID19.R1.notmerged.fq.gz > "+ sample+'/'+sample+".COVID19.FINAL.fq"
print(zcat)
subprocess.call(zcat, shell=True)

#bowtie2
bowtie = 'bowtie2 --end-to-end -x '+index + ' -q ' +sample+'/'+sample+".COVID19.FINAL.fq -S "+ sample+'/'+sample+".COVID19.FINAL.fq.sam"
print(bowtie)
subprocess.call(bowtie, shell=True)

#covnert to BAM 
bam = 'samtools view -S -b '+sample+'/'+sample+".COVID19.FINAL.fq.sam > "+sample+'/'+sample+".COVID19.FINAL.fq.sam.bam"
print(bam)
subprocess.call(bam, shell=True)

#SORT1
sort = 'samtools sort '+sample+'/'+sample+".COVID19.FINAL.fq.sam.bam -o "+sample+'/'+sample+".COVID19.FINAL.fq.sam.sorted.bam"
print(sort)
subprocess.call(sort, shell=True)

#index
bowtie_index = 'samtools index '+sample+'/'+sample+".COVID19.FINAL.fq.sam.sorted.bam"
print(bowtie_index)
subprocess.call(bowtie_index, shell=True)

print('DONE! :)'+'\n'+'\n')
print('Final file name is '+'\t'+sample+'/'+sample+".COVID19.FINAL.fq")
