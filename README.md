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

**1.**




