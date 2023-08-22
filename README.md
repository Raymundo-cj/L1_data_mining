## L1_data_mining

In this note, i will talk about how to complete data mining.

### 1.data processing
in this part, I had download a big data from NCBI nr database,so the first step is to process it.
There is a reference code for how to download this dataset.（This is a traditional method,so it will take a long time）
```
#!/bin/bash
<<COM
for i in {00..67};
do wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/nr.$i.tar.gz;
   wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/nr.$i.tar.gz.md5;
   md5sum -c nr.$i.tar.gz.md5;
   tar -zxvf nr.$i.tar.gz.md5 -C;
   echo "nr.$i has done";
done
COM
wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5
```
**1.data splite**

let data splite 9 parts.
```
#!/bin/bash
perl fasta-splitter.pl --n-parts 9 ~/caojian/blast/nr.fasta
echo 'finished'
```
**2.len>=200**

Amino acid length less than 200 has no effect on the analysis and needs to be removed.
```
#!/bin/bash
python find_protein_len200.py nr.part-1.fasta nr-1_len200.fasta &
··· &
wait
echo "finished"
```
**3.remove repeats**

blast will report errors when encountering sequences with large segment repeats, and large segment repeats need to be filtered out.
```
#!/bin/bash
python delete_repeat_fasta3.py nr-1_len200.fasta nr_1_report.fasta &
... &
wait
echo "finished"
```
### 2.data align

**1.psi-blast**
```
makeblastdb -in /data/caojian/L1/seed.fasta -dbtype prot -title seed.pep -parse_seqids -out seed.pep
# 构建比对的种子序列数据库
```
```
#!/bin/bash
psiblast -query /data/caojian/L1/ncbi/fast_file/nr_1_report.fasta -evalue 1e-5 -inclusion_ethresh .002 -db /data/caojian/L1/blast/seed.pep -num_iterations 5 -seg yes -outfmt '7 std qseq sseq stitle' -out nr_bl_1.output -max_target_seqs 500 
```
这一步执行结束后会得到这样的一些结果：

![image](https://github.com/Raymundo-cj/L1_data_mining/assets/64938817/fd77b331-9639-4651-af35-89cec8187813)

根据下面的代码筛选出比对结果中的蛋白序列号及蛋白序列文件：

```
sed '/#/d' nr_bl_1.output|awk '{print $1}'|sort|uniq > nr_bl_1.list

seqtk subseq /data/caojian/L1/ncbi/fast_file/nr_1_report.fasta nr_bl_1.list > nr_bl_1.fasta

cat nr_bl_1.fasta ... >> nr_bl.fasta
```

**2.mmseq**

```
mmseqs createdb /public1/home/scb8190/yumeixia/Ago/PIWI/database/pAgo_PIWI.fasta pAgo_DB
mmseqs createdb ~/shared_public_data/metagenome_data//public1/home/scb8190/shared_public_data/metagenome_data/split_12_repeat_filter/metagenome.pep_9_filter.fasta metagenome.pep_1_DB
mmseqs search metagenome.pep_1_DB pAgo_DB result_1_DB 1.tmp --start-sens 4 --sens-steps 4 -s 7 -e 1e-6
mmseqs convertalis metagenome.pep_1_DB pAgo_DB result_1_DB  mmseqs_1_results.txt
```
### 3.后续分析

**1.pfam注释**

```
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
hmmscan --domtblout nr_blast.dbl --tblout nr_blast.tbl --noali -E 1e-2 --cpu 64 /public1/home/scb8190/Test/CRISPR/Pfam/Pfam-A.full.uniprot.hmm nr_bl.fasta
echo "finished"
```
根据pfam结果筛选含有“Endonuclease-reverse transcriptase”domain的序列

```
grep 'Endonuclease-reverse transcriptase' nr_blast.dbl|awk '{print $4}'|sort|uniq >domain_bl.list
seqtk subseq nr_bl.fasta domain_bl.list > nr_bl_pfam.fasta
```
**2.比对物种全基因组**
```
makeblastdb -in Arabidopsis_thaliana.fasta -dbtype nucl -parse_seqids -input_type fasta -out Arabidopsis_thaliana

tblastn -db ~/caojian/non-ltr/orf2_tblastn/gene/Arabidopsis_thaliana/Arabidopsis_thaliana -query ~/caojian/L1/seed/nr_bl_pfam.fasta -outfmt 6 -out tblastn_Arabidopsis_thaliana.result
```
根据上面得到的tblastn结果，将比对的序列根据比对长度筛选大于800，获取其在基因组上的位置，在对比位置处上下游各延长2K，得到后续的文件

```
mkdir result
python from_outresult_getfa.py tblastn_Arabidopsis_thaliana.result ~/caojian/non-ltr/orf2_tblastn/gene/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
# python file save at Data_processing
```





