from Bio import SeqIO
import sys
# 定义一个函数，用于检测是否存在重复片段
def check_repeat(seq, n):
    """
    检查序列中是否存在长度为n的重复片段
    """
    for i in range(len(seq) - n):
        if seq[i:i+n] in seq[i+n:]:
            return True
    return False

# 读取FASTA文件
input_file =sys.argv[1]
output_file =sys.argv[2]

with open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        # 检查序列是否包含重复片段
        if not check_repeat(str(record.seq), 10): #在此处，长度为20的片段被认为是“重复片段”。
            # 将没有重复片段的蛋白质序列写回到输出文件中
            SeqIO.write(record, output_handle, "fasta")
