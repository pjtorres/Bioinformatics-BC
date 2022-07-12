import os, argparse



# infile='1RAD_1027298.star.starLog.final.out'
# outfile='1RAD_1027298.star.starLog.final.out.percent.tsv'

#--------Command line arguments----------
parser=argparse.ArgumentParser(description="Script will allow you to percent mapped and unmaped reads in a STAR output log file.")
parser.add_argument('-i','--input', help='File should have suffuc star.starLog.final.out. ', required=True)



#Pass arguments
args = parser.parse_args()
infile=args.input
outfile = infile+'.percent.tsv'


fin = open(infile, 'r+')
fout = open(outfile,'w+')
for line in fin:
    line=line.strip().split('|')
    if line[0]=='Number of input reads ':
        total = line[1].strip()
    if line[0]=='Uniquely mapped reads % ':
        unique_mapped = line[1].strip()
    if line[0]=='% of reads mapped to multiple loci ':
        multi_mapped = line[1].strip()
    if line[0]=='% of reads mapped to too many loci ':
        toomany_mapped = line[1].strip()

    if line[0]=='% of reads unmapped: too short ':
        too_short = line[1].strip()

    if line[0]=='% of reads unmapped: other ':
        other_unmapped = line[1].strip()

total_unmapped = str(float(other_unmapped.strip('%')) +float(too_short.strip('%')))
total_mapped = str(float(unique_mapped.strip('%')) +float(multi_mapped.strip('%')) + float(toomany_mapped.strip('%')))


fout.write('percent_mapped'+'\t'+'percent_unmapped'+'\n')
fout.write(total_mapped+'\t'+total_unmapped)

fout.close()
fin.close()

        
