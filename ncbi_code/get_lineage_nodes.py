import pandas as pd
import argparse
import os



"""
Need to have two seq id , parental id, and names files already made. Below is an example of how you can make one
! echo -e "seq_id\tparent_id\tname" > species_seqid.txt
!grep species nodes.dmp | awk '{print $1, $2, $13}' FS='|' OFS='\t' | sed 's/\t\{2,\}/\t/g; s/\t$//' >> species_seqid.txt

! echo -e "seq_id\tparent_id\tname" > all_seqid.txt
!grep -v species nodes.dmp | awk '{print $1, $2, $13}' FS='|' OFS='\t' | sed 's/\t\{2,\}/\t/g; s/\t$//' >> all_seqid.txt
"""
# --------Command line arguments----------
parser = argparse.ArgumentParser(description="Make lineage tree from nodes file. ")
parser.add_argument('-s', '--species', help='Species seq and parental id name file.', required=True)
parser.add_argument('-a', '--all', help='All taxa except species seq and aprental id name file', required=True)
parser.add_argument('-o', '--output', help='The output file name.', required=True)

# Pass arguments
args = parser.parse_args()
species = args.species
all_l = args.all
outputf = args.output

species = pd.read_csv(species,sep='\t')
nodes = pd.read_csv(all_l, sep='\t')


seqid_name_dict = {}
df1 = species
df2 = nodes

# Iterate through each row in df1
for index, row in df1.iterrows():
#     print(row)
    seq_id = row["seq_id"]
    parent_id = row["parent_id"]
    name = row["name"]
    
    # Check if seq_id is not in the dictionary
    if seq_id not in seqid_name_dict:
        og_seqid= seq_id
        taxa_list = []
        taxa_list.append(name)
        
        # Iterate up the tree until seq_id is 1
        while seq_id != 1:
            name_row = df2[df2["seq_id"] == parent_id]
            if not name_row.empty:
                parent_id = name_row["parent_id"].values[0]
                name = name_row["name"].values[0]
                taxa_list.append(name)
                seq_id = name_row["seq_id"].values[0]
            else:
                break
        
        # Reverse the list to get the top-down order
        taxa_list.reverse()
        
        # Update the dictionary with seq_id and its corresponding name list
        seqid_name_dict[og_seqid] = taxa_list

# Print the seqid_name_dict
fout = open(outputf+'.txt','w+')
fout.write('species\tlineage\n')
for seq_id, taxa_list in seqid_name_dict.items():
    taxa_list = taxa_list[1:]
    sp = taxa_list[-1].split('s__')[1]
    print(sp + '\t'+';'.join(taxa_list)+'\n')
    fout.write(sp+'\t'+';'.join(taxa_list)+'\n')
    
fout.close()

print('\nDONE!\n')
