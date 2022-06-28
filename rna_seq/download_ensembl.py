#!/usr/bin/python

#Created: 2021-04-15
#Author: Francine R. Camacho
#Last updated:  2021-04-15

from zipfile import ZipFile
import sys
import argparse
import wget
import os
import gzip
import shutil

# Python usage: python download_ensembl.py -i 101


"""
Unzipping Database: Function to unzip both Ensembl .fna and GTF Human Genome .gz files
"""
def gz_extract(directory):

    extension = ".gz"
    os.chdir(directory)

    for item in os.listdir(directory): # loop through items in dir
      if item.endswith(extension): # check for ".gz" extension
          gz_name = os.path.abspath(item) # get full path of files
          file_name = (os.path.basename(gz_name)).rsplit('.',1)[0] #get file name for file within
          print('\tUnzipping ', gz_name, '---> ', file_name,' ....')

          with gzip.open(gz_name,"rb") as f_in, open(file_name,"wb") as f_out:
              shutil.copyfileobj(f_in, f_out)
          os.remove(gz_name) # delete zipped file

"""
Database construction: Function downloads both Ensembl .fna and GTF Human Genome
"""
def download_ensembl(release_number, outdir):

    if outdir is None:
        outdir = os.getcwd()
        data_dir = os.path.join(os.path.dirname(outdir),'data')
    else:
        data_dir = outdir

    try:
        main_dir = os.path.join(data_dir,'Ensembl_'+ release_number)
        os.makedirs(main_dir, 0o777, True)
        print("Directory '% s' is created ...." % main_dir)

    except OSError as error:
        print("Directory '% s' already exist \n \tDownloading Ensembl there ...." % main_dir)
        pass

    print('Begin Ensembl download ....')

    ensembl_url = 'ftp://ftp.ensembl.org/pub/'
    human_genome = os.path.join(ensembl_url, "release-"+ release_number, "fasta/homo_sapiens/dna", "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
    human_gtf = os.path.join(ensembl_url, "release-"+ release_number, "gtf/homo_sapiens", "Homo_sapiens.GRCh38." + release_number+".gtf.gz")

    try:
        print("\tDownloading Ensembl human genome ....")
        # Download the latest human genome fasta file from Ensembl. It is recommended to use the "primary_assembly"
        wget.download(human_genome, out=main_dir)
    except:
        print (e)

    try:
        print("\n \tDownloading Ensembl annotation human genome ....")
        # Download the corresponding GTF annotation file for the human genome from Ensembl.
        wget.download(human_gtf, out=main_dir)
    except:
        print (e)

    print('Done with Ensembl database download!')
    gz_extract(main_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download the latest Ensembl database")
    parser.add_argument('-i', dest='release', help='Ensembl database release number', required=True)
    parser.add_argument('-o', dest='outdir', help='Directory for Ensembl database', required=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
        download_ensembl(args.release, args.outdir)
