from Bio import SeqIO
import os
import glob

# 定义fasta文件列表
fasta_folder="/public1/home/scb8190/caojian/L1/seed/Arabidopsis_thaliana/Arabidopsis_thaliana_800/orf_result"
fasta_files = glob.glob(os.path.join(fasta_folder,"*.fasta"))

# 定义输出结果的TXT文件路径
output_file = "output.txt"

# 打开输出文件
with open(output_file, "w") as f:
    # 遍历每个fasta文件
    for fasta_file in fasta_files:
        # 记录当前fasta文件的序列数和长度
        seq_count = 0
        seq_length = 0

        # 打开fasta文件
        with open(fasta_file) as handle:
            # 使用SeqIO解析fasta文件
            for record in SeqIO.parse(handle, "fasta"):
                # 增加序列计数
                seq_count += 1
                # 增加序列长度
                seq_length += len(record.seq)
                # 写入序列编号和长度到输出文件
                f.write(f"{record.id}\t{len(record.seq)}\n")

        # 输出fasta文件的统计信息
        f.write(f"文件: {fasta_file}\n")
        f.write(f"总序列数: {seq_count}\n")
        f.write("\n")
