#!/usr/bin/python

#Created: 2021-04-15
#Author: Francine R. Camacho
#Last updated:  2021-04-15

import sys
import argparse
import os
import subprocess

"""
Database construction:
Derive the transcript fasta file from the Ensembl human genome fasta file using the Ensembl GTF file.
The RSEM package has a function called rsem-prepare-reference which will take a genome fasta file and a GTF file (--gtf argument therein)
and build an index for RSEM. The rsem-prepare-reference command will also produce a transcript fasta file obtained by applying the GTF annotations
to the genome.
Requirements: RSEM is downloaded and set to $PATH
"""
# Python usage: python generate_transcriptome.py -i 103 -o gxi_index_2.1.transcriptome

def prepare_reference(release_number, refname):

    current_dir = os.getcwd()
    data_dir = os.path.join(os.path.dirname(current_dir),'data')
    ensembl_dir = os.path.join(data_dir, 'Ensembl_' + release_number)
    ensemble_genome = os.path.join(ensembl_dir, 'Homo_sapiens.GRCh38.dna.primary_assembly.fa')
    ensemble_gtf = os.path.join(ensembl_dir, "Homo_sapiens.GRCh38." + release_number+".gtf")
    outdir = os.path.join(data_dir, 'Ensembl_' + release_number)

    if not os.path.exists(ensemble_genome) and not os.path.exists(ensemble_gtf):
        print('Ensembl database is missing ....')
        print('\t'+ ensemble_genome + '\n', '\t'+ensemble_gtf)
        print('Please download Ensembl database using download_ensembl.py')

    else:
        os.chdir(outdir)
        cmd = "rsem-prepare-reference --gtf " + ensemble_gtf + " "+ ensemble_genome + " " + refname + " > " + ensembl_dir + "/rsem-prepare.log"
        print("Running RSEM:", cmd)
        subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create transcriptome files from Ensembl database")
    parser.add_argument('-i', dest='release', help='Ensembl database release number', required=True)
    parser.add_argument('-o', dest='refname', help='outfile name for rsem-prepare-reference file', required=True)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
        prepare_reference(args.release, args.refname)
