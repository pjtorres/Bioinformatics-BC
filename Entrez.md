I use [NCBI Entrez Direct UNIX E-utilities regularly](https://www.ncbi.nlm.nih.gov/books/NBK179288/) for sequence and data retrieval from NCBI. These UNIX utils can be combined with any UNIX commands

# Looking for specific bugs:
```bash
esearch -db "nucleotide" -query "Faecalibacterium prausnitzii[ORGN]"|  efetch -format fasta
```

# rhino virus:

```bash
esearch -db "nucleotide" -query "Rhinovirus[ORGN]"|  efetch -format fasta | grep '>' | head
```
# specific genes:
```bash
esearch -db gene -query "Liver cancer AND Homo sapiens" |elink -target nuccore | efetch -format fasta
```
# filter all bacteria with a filter
```bash
esearch  -db "nucleotide" -query "Bacteria[Organism] AND Refseq[Filter]" | efetch -format fasta  
```
# download based on accession
```bash
esearch -db assembly -query GCF_000508965.1 | elink -target nucleotide -name assembly_nuccore_insdc | efetch -format fasta > GCF_000508965.1.fna
```

# get assembly status
```bash
search -db assembly -query "Veillonella sp. DNF00869" |  esummary |xtract -pattern DocumentSummary -element SpeciesName,assembly-status >> assembly_status.txt
```

[Great resource for downloading Entrez on your compute](https://dataguide.nlm.nih.gov/edirect/install.html)

[Tips on efetch values](https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly)
