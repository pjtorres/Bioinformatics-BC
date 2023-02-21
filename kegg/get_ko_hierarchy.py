#Created: 2023-02-14
#Author: Pedro J. Torres
#Last updated:
import os, argparse
import Bio.KEGG
from Bio.KEGG import REST
import argparse, re
import urllib.request
from urllib.error import HTTPError
import pandas as pd

#--------Command line arguments----------
parser=argparse.ArgumentParser(description="Get Ko symbol, function and pathway given a file where each row is a KO.")
parser.add_argument('-ko','--ko_file', help='file where each row is a KO.', required=True)
parser.add_argument('-o','--out', help='Output filename.', required=True)


#Pass arguments
args = parser.parse_args()
ko_file=args.ko_file
out=args.out


def get_ko_hierarchy(ko):
    ko_hierachy = []
    result = REST.kegg_get(ko)
    flag = False
    kodef = REST.kegg_list(ko).read()
    function_def  = kodef.split('\t')[1].split(';')[1].strip()
    symbol = kodef.split(';')[0].split('\t')[1].strip()
    for li in result:

        if flag and ('[BR:ko' in li or li[0]!=' '):

            break
        if flag:
            if ko in li:
                ko_hierachy.append(['ko:'+ko,symbol,function_def, l3, l2, l1])
            elif li[14] == ' ':
                l3 = li.strip('\n').strip('0123456789 ')
            elif li[13] == ' ':
                l2 = li.strip('\n').strip('0123456789 ')
            elif li[12] == ' ':
                l1 = li.strip('\n').strip('0123456789 ')

        else:
            if 'KEGG Orthology (KO)' in li:
                flag = True
    return (ko_hierachy)

print('loading.\n')
df = pd.read_csv(ko_file, header=None)
ko_list = df[0].values.tolist()
result = []
print('getting hierarchy. \n')
for ko in ko_list:
    try:
        result.extend(get_ko_hierarchy(ko))
    except:
        print(ko, 'not found')
pd.DataFrame(result, columns=['ko', 'symbol','function_def','l3', 'l2', 'l1'])

test_ko_df = pd.DataFrame(result, columns=['ko','symbol' ,'function_def','l3', 'l2', 'l1'])

test_ko_df.to_csv(out+'.tsv', sep='\t', index=False)

print('\nDone!\n')
