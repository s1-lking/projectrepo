# BIO312 Fall 2024 Project Repository

## Lab 3

1. uncompress the proteomes
```bash
gunzip proteomes/*.gz
```
2. put all the protein sequences into a single file
```bash
cat  proteomes/*.faa > allprotein.fas
```
3. build a blast database
```bash
makeblastdb -in allprotein.fas -dbtype prot
```
4. create a folder for your protein
```bash
mkdir ~/lab03-$MYGIT/smc3
```
5. Go to the folder
```bash
cd ~/lab03-$MYGIT/smc3
```
6. download query protein
```bash
ncbi-acc-download -F fasta -m protein "NP_005436.1"
```
7. perform blast search using the query protein
```bash
blastp -db ../allprotein.fas -query NP_005436.1.fa -outfmt 0 -max_hsps 1 -out smc3.blastp.typical.out
```
8. Perform blast search and request tabular output
```bash
blastp -db ../allprotein.fas -query NP_005436.1.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out smc3.blastp.detail.out
```
9. Filter output for high scoring putative homologs with e-value cutoff of < 1e-30
```bash
awk '{if ($6< 1e-30)print $1 }' smc3.blastp.detail.out > smc3.blastp.detail.filtered.out
```
10. Count total number of hits in the BLAST results after the filter
```bash
wc -l smc3.blastp.detail.filtered.out
```
11. find the number of paralogs in each species
```bash
grep -o -E "^[A-Z]\.[a-z]+" smc3.blastp.detail.filtered.out  | sort | uniq -c
```

## Lab 4

1. obtain the sequences that are in the BLAST output file
```bash
seqkit grep --pattern-file ~/lab03-$MYGIT/smc3/smc3.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/smc3/smc3.homologs.fas
```
2. make a multiple sequence alignment using muscle
```bash
muscle -align ~/lab04-$MYGIT/smc3/smc3.homologs.fas -output ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas
```
3. View the alignment and save it as a pdf
```bash
alv -kli  ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas | less -RS
```
```bash
alv -kli --majority ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas | less -RS
```
4. print your alignment to a large pdf file
```bash
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas
```
5. calculate the width (length) of the alignment
```bash
alignbuddy  -al  ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas
```
6. calculate the length of the alignment after removing any column with gaps
```bash
alignbuddy -trm all  ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas | alignbuddy  -al
```
7. calculate the length of the alignment after removing invariant (completely conserved) positions
```bash
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas | alignbuddy  -al
```
8. calculate the average percent identity using t_coffee
```bash
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas -output sim
```
9. calculating the average percent identity using alignbuddy
```bash
 alignbuddy -pi ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```

## Lab 5

1. remove any sequence that contains a duplicate label tag (I've indicated we did this by appending an f to homologsf in the output file name), and put a copy in your lab05 directory
```bash
sed 's/ /_/g'  ~/lab04-$MYGIT/smc3/smc3.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/smc3/smc3.homologsf.al.fas
```
2. find the maximum likehood tree estimate using IQ-TREE
```bash
iqtree -s ~/lab05-$MYGIT/smc3/smc3.homologsf.al.fas -bb 1000 -nt 2 
```
3. display the tree unrooted with a graphical display
```bash
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R  ~/lab05-$MYGIT/smc3/smc3.homologsf.al.fas.treefile ~/lab05-$MYGIT/smc3/smc3.homologsf.al.fas.treefile.pdf 0.4 15
```
4. midpoint rooting using gumtree
```bash
gotree reroot midpoint -i ~/lab05-$MYGIT/smc3/smc3.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/smc3/smc3.homologsf.al.mid.treefile
```
5. display the midpoint rooted tree as a svg
```bash
nw_order -c n ~/lab05-$MYGIT/smc3/smc3.homologsf.al.mid.treefile  | nw_display -w 1000 -b 'opacity:0' -s  >  ~/lab05-$MYGIT/smc3/smc3.homologsf.al.mid.treefile.svg -
```
6. convert the midpoint rooted tree svg into a pdf
```bash
convert ~/lab05-$MYGIT/smc3/smc3.homologsf.al.mid.treefile.svg  ~/lab05-$MYGIT/smc3/smc3.homologsf.al.mid.treefile.pdf
```
7. display the midpoint rooted tree as a cladogram
```bash
nw_order -c n ~/lab05-$MYGIT/smc3/smc3.homologsf.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > ~/lab05-$MYGIT/smc3/smc3.homologsf.al.midCl.treefile.svg -
convert ~/lab05-$MYGIT/smc3/smc3.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/smc3/smc3.homologsf.al.midCl.treefile.pdf
```

## Lab 6

1. Reconcile the gene and species tree using Notung
```bash
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/smc3/smc3.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/smc3/
```
2. generate a RecPhyloXML object and view the gene-within-species tree via thirdkind
```bash
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/smc3/smc3.homologsf.al.mid.treefile.rec.ntg --include.species
```
3. create a gene-reconciliation-within species tree reconciliation graphic using thirdkind
```bash
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/smc3/smc3.homologsf.al.mid.treefile.rec.ntg.xml -o  ~/lab06-$MYGIT/smc3/smc3.homologsf.al.mid.treefile.rec.svg
```
4. convert the reconciled tree svg into a pdf
```bash
convert  -density 150 ~/lab06-$MYGIT/smc3/smc3.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/smc3/smc3.homologsf.al.mid.treefile.rec.pdf
```

## Lab 8

1. make a copy of our raw unaligned sequence, removing the asterisk (stop codon) in the process
```bash
sed 's/*//' ~/lab04-$MYGIT/smc3/smc3.homologs.fas > ~/lab08-$MYGIT/smc3/smc3.homologs.fas
```
2. run RPS-Blast using the unaligned sequences as query and the given Pfam database
```bash
rpsblast -query ~/lab08-$MYGIT/smc3/smc3.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/smc3/smc3.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
3. Copy the final gene tree over from lab 5
```bash
cp ~/lab05-$MYGIT/smc3/smc3.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/smc3
```
4. plot the pfam domain predictions from rps-blast next to their cognate protein on the phylogeny
```bash
Rscript  --vanilla ~/lab08-$MYGIT/plotTreeAndDomains2.r ~/lab08-$MYGIT/smc3/smc3.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/smc3/smc3.rps-blast.out ~/lab08-$MYGIT/smc3/smc3.homologs.fas ~/lab08-$MYGIT/smc3/smc3.tree.rps.pdf
```
5. view the tab delimited annotations using the script below or viewing the pdf made in the previous step
```bash
mlr --inidx --ifs "\t" --opprint  cat ~/lab08-$MYGIT/smc3/smc3.rps-blast.out | tail -n +2 | less -S
```
6. examine the RPS-output file for proteins that have more than one annotation
```bash
cut -f 1 ~/lab08-$MYGIT/smc3/smc3.rps-blast.out | sort | uniq -c
```
7. examine the RPS-output file for the most commonly found Pfam domain annotation
```bash
cut -f 6 ~/lab08-$MYGIT/smc3/smc3.rps-blast.out | sort | uniq -c
```
8. examine the RPS-output file for the longest annotated protein domain
```bash
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/smc3/smc3.rps-blast.out |  sort  -k2nr
```
9. examine the RPS-output file for the protein that has a domain annotation with the best e-value
```bash
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/smc3/smc3.rps-blast.out
```
