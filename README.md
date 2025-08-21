## Overview
1. Download genome from NCBI
```
datasets download genome accession GCA_016772135.1 --include gff3,rna,cds,protein,genome,seq-report
unzip ncbi_dataset.zip
```
2. Build database
```
python gff2db.py \
    --gff ncbi_dataset/data/GCA_016772135.1/genomic.gff 
    --fasta ncbi_dataset/data/GCA_016772135.1/GCA_016772135.1_ASM1677213v1_genomic.fna 
    --gene FKS1
```
3. Add hotspot regions to database
- Open genes.json
- Update the `target` field with the text below
```
    "targets": [
      {
        "name": "hs1",
        "coords": [1902, 1928]
      },
      {
        "name": "hs2",
        "coords": [4047, 4071]
      },
      {
        "name": "hs3",
        "coords": [2070, 2073]
      }
    ]
```
4. Detect variants
Example reads can be downloaded with `fasterq-dump -S SRR14996559`
```
./fungalAMR \\
    --db genes.json \\
    --r1 SRR14996559_1.fastq.gz \\
    --r2 SRR14996559_2.fastq.gz \\
    --outdir results
```
5. Check results
Mutations can be found in `results/mutations.csv`.


## Using docker image
```
docker run --rm \
  -u $(id -u):$(id -g) \
  -e HOME=/tmp \
  -e XDG_CACHE_HOME=/tmp/.cache \
  -e MAMBA_NO_RC=true \
  -v "$PWD":/data \
  fungalamr:latest \
  --db /data/genes.json --r1 /data/SRR14996559_1.fastq.gz --r2 /data/SRR14996559_2.fastq.gz --outdir /data/results
```