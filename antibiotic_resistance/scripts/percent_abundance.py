__author__= 'Pedro J. Torres'
import argparse
import os
import pandas,numpy

#-----------Command Line Arguments---------------------------------------
parser=argparse.ArgumentParser(description="Script allows you to convert raw counts from diamond and get relative abudnace.")
parser.add_argument('-i','--input', help=' Input txt file you want covnert to percentage',required=True)
parser.add_argument('-o','--out', help='Name of output file: jsut the mae e.g., "Percent_abundance"', required=True)#require later
args = parser.parse_args()
o_file=str(args.out)
txtfile=str(args.input) 

#-----------open text file in with pandas -------
print ("Retrieving percent abundance")
df=pandas.read_table(txtfile)
cols = [c for c in df.columns ] # get column names
df=df[cols]# convert columns to dataframe

#--------------- Convert raw counts to percentages -------
dft=(df.set_index('Atb_gene').T)
cols= df['Atb_gene'].tolist()
dft[cols] = dft[cols].div(dft[cols].sum(axis=1), axis=0).multiply(100)# remove this '.multiple(100)' if you only want relative abudnace
dft=dft.transpose()
dft.to_csv(o_file+'.csv')

print ("Done :)")
