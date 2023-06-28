import argparse
import os
import fnmatch
import pandas as pd
import boto3

# python download_aws_s3.py -b persephone -f data_warehouse/illumina_runs/wgs65/trimmed_default/ -p PBT-03090*gz -o wgs65
# --------Command line arguments----------
parser = argparse.ArgumentParser(description="Script to download R1 and R2 from AWS S3")
parser.add_argument('-b', '--bucket', help='The name of the S3 bucket', required=True)
parser.add_argument('-f', '--folder', help='The folder path in the S3 bucket', required=True)
parser.add_argument('-p', '--pattern', help='The partial file name pattern to match', required=True)
parser.add_argument('-o', '--output', help='The output directory path', required=True)

# Pass arguments
args = parser.parse_args()
bucket_name = args.bucket
s3_folder = args.folder
pattern = args.pattern
output_dir = args.output

def create_output_directory(output_dir):
    """
    Create the output directory if it does not exist
    Args:
        output_dir: the output directory path
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Output directory created: {output_dir}")
    else:
        print(f"Output directory already exists: {output_dir}")

def extract_sample_name(pattern):
    """
    Extracts the sample name from the pattern if it contains "*"
    Args:
        pattern: the partial file name pattern
    Returns:
        sample_name: the extracted sample name
    """
    if '*' in pattern:
        sample_name = pattern.split('*')[0]
        return sample_name
    else:
        return None

def download_s3_objects(bucket_name, s3_folder, pattern, output_dir):
    """
    Download S3 objects based on a partial file name pattern
    Args:
        bucket_name: the name of the S3 bucket
        s3_folder: the folder path in the S3 bucket
        pattern: the partial file name pattern to match
        output_dir: the output directory path
    Returns:
        r1_files: list of filenames with "R1" in their basename
        r2_files: list of filenames with "R2" in their basename
    """
    print('Running')
    s3 = boto3.client('s3')

    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=s3_folder)
    r1_files = []
    r2_files = []
    if 'Contents' in response:
        for obj in response['Contents']:
            key = obj['Key']
            if not obj['Key'].endswith('/'):
                if fnmatch.fnmatch(os.path.basename(key), pattern):
                    target = os.path.join(output_dir, os.path.basename(key))
                    s3.download_file(bucket_name, key, target)
                    print(f"Downloaded: {key} -> {target}")
                    if 'R1' in os.path.basename(key):
                        r1_files.append(os.path.join(output_dir, os.path.basename(key)))
                    if 'R2' in os.path.basename(key):
                        r2_files.append(os.path.join(output_dir, os.path.basename(key)))

    return r1_files, r2_files

def write_file(sample_name, r1_files, r2_files):
    """
    Write the filenames to a file
    Args:
        sample_name: the sample name
        r1_files: list of filenames with "R1" in their basename
        r2_files: list of filenames with "R2" in their basename
    """
    output_file = os.path.join(output_dir, f"{sample_name}.txt")
    with open(output_file, 'w') as file:
        for r1_file, r2_file in zip(r1_files, r2_files):
            file.write(f"{sample_name}\t{r1_file}\t{r2_file}\n")
    print(f"File written: {output_file}")

create_output_directory(output_dir)
sample_name = extract_sample_name(pattern)
r1_files, r2_files = download_s3_objects(bucket_name, s3_folder, pattern, output_dir)
write_file(sample_name, r1_files, r2_files)

print('Done!')
