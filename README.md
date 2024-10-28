# BIO312 Fall 2024 Project Repository

## Lab 3

1. uncompress the proteomes
```
gunzip proteomes/*.gz
```
2. put all the protein sequences into a single file
```
cat  proteomes/*.faa > allprotein.fas
```

3. build a blast database
```
makeblastdb -in allprotein.fas -dbtype prot
```
4. create a folder for your protein
```
mkdir ~/lab03-$MYGIT/smc3
```
5. Go to the folder
```
cd ~/lab04-$MYGIT/smc3
```
6. download query protein
```
ncbi-acc-download -F fasta -m protein "NP_005436.1"
```
7. perform blast search using the query protein
```
blastp -db ../allprotein.fas -query NP_005436.1.fa -outfmt 0 -max_hsps 1 -out smc3.blastp.typical.out
```
8. Perform blast search and request tabular output
```
blastp -db ../allprotein.fas -query NP_005436.1.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out smc3.blastp.detail.out
```
9. Filter output for high scoring putative homologs with e-value cutoff of < 1e-30
```
awk '{if ($6< 1e-30)print $1 }' smc3.blastp.detail.out > smc3.blastp.detail.filtered.out
```
10. Count total number of hits in the BLAST results after the filter
```
wc -l smc3.blastp.detail.filtered.out
```
11. find the number of paralogs in each species
```
grep -o -E "^[A-Z]\.[a-z]+" smc3.blastp.detail.filtered.out  | sort | uniq -c
```
