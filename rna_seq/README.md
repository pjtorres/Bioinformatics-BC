# Comparing genome to transcripts in STAR


1. Download the latest human genome fasta file from Ensembl. It is recommended to use the "primary_assembly" version, e.g. (ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
2. Download the corresponding GTF annotation file for the human genome from Ensembl. e.g. ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

## Need to have
* [STAR installed](https://github.com/alexdobin/STAR)
* human genome indexed using STAR i.e. 

```bash
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles star_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa  --sjdbGTFfile star_index/Homo_sapiens.GRCh38.99.gtf  --sjdbOverhang 100 
```

and its associated gtf file

## Run STAR on fastq PE files. 
```bash
./run_STAR.sh -x /source/star_index/ -g /source/star_index/Homo_sapiens.GRCh38.99.gtf -1 <sample_R1.fastq> -2 <sample_R2.fastq> -o <sample_prefix> -p 8
```

## Get bam file of reads uniquely mapping to either known human transcripts or human genome only
```bash
./filter_genome_STAR.sh -t <sample_prefix>.starAligned.toTranscriptome.out.bam  -g <sample_prefix>.starAligned.sortedByCoord.out.bam  -o <sample_prefix>
```

**There will be two output directories**
1. <sample_prefix>.IN.transcript - bam and bai file with reads alinging to human trancriptome 
2. <sample_prefix>.NA.trancript.IN.genome - bam and bai file with reads aligning to human genome but not transcriptome

# If you want a quicker pseudo mapping alternative, that doesn’t align to the human genome but quickly to the human transcriptome. You can use Salmon. The results will account for GC bias, and will count the transcripts and genome for you using an expectation maximization algorithm to deal with multi-mapped transcripts. Below is an example script you can use

`bash
salmon quant			\
--libType	$libtypesalmon	\
-i		"$indexdir"	$infqarg	\
--seqBias			\
--gcBias			\
-p		"$numthreads"	\
-g		"$genemapfile"	\
--discardOrphansQuasi		\
--validateMappings		\
-o		"$outdir"	`
