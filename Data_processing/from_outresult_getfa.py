import csv,sys,pysam,os

input_file=sys.argv[1]
fasta=sys.argv[2]
output_folder = "result"

def get_seq(chr,start,end):
	fasta_open = pysam.Fastafile(fasta)
	seq = fasta_open.fetch(chr,start,end)
	fasta_open.close()
	return seq

with open(input_file) as f:
   reader = csv.reader(f, delimiter='\t')
   for row in reader:
     if int(row[3]) > 800:
        if int(row[8]) >int(row[9]):
           start=int(row[9]) -2000
           end=int(row[8]) +2000
        else:
           start=int(row[8]) -2000
           end=int(row[9]) +2000
        seq_id=row[1]
        tem = get_seq(seq_id,start,end)
        file_path = os.path.join(output_folder, seq_id + "_" + str(start) + "_" + str(end) + ".fasta")
        file = open(file_path,"w")
        file.write(">" + seq_id + str(start) + "-" + str(end) + "\n")
        file.write(tem)
        file.close()
