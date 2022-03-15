#import snakemake
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

x=[]

for i in snakemake.input:
    for record in SeqIO.parse(snakemake.input[0], "fasta"):
        def replace_stop_codons(record, codon_stop_array = ["TAG", "TGA", "TAA"]):
            tempRecordSeq = list(record.seq)
            for index in range(0, len(record.seq), 3):
                    codon = record.seq[index:index+3]
                    if codon in codon_stop_array:
                        tempRecordSeq[index:index+3] = '?','?','?'
            record.seq = Seq("".join(tempRecordSeq))
            x.append(record)
        replace_stop_codons(record)
    with open(snakemake.output[0], "w") as output_handle:
        SeqIO.write(x, output_handle, "fasta")
