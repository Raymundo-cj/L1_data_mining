from Bio import SeqIO

import sys
input_file =sys.argv[1]
output_file =sys.argv[2]
# 读取蛋白序列文件
sequences = SeqIO.parse(input_file, "fasta")

# 定义最小和最大序列长度
min_length = 200

# 筛选符合长度要求的序列
filtered_sequences = [seq_record for seq_record in sequences if min_length <= len(seq_record.seq)]

SeqIO.write(filtered_sequences, output_file, "fasta")