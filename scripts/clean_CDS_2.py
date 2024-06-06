# this writes only files which do not contain "N" nor "X" nor "?" in any record of the fasta sequence imported with seqio

import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools
import os

os.makedirs(snakemake.output[0])

input=snakemake.input

for i in input:
    file=os.path.basename(i)
    cumlen=0
    x=[]
    for record in SeqIO.parse(i, "fasta"):
        sequence=record.seq
#        print(sequence)
        # pattern matching "N"
        lenN = len(re.findall("N", str(sequence)))
        lenX = len(re.findall("X", str(sequence)))
        lenQM = len(re.findall("\?", str(sequence)))
        cumlen = cumlen + lenN + lenX + lenQM
#        print(lenN)
        x.append(record)
    print("cumlen is", cumlen)
#select only input files which do not contain "N" in any of their record's sequence
    if cumlen == 0:
#write it to output file only if it does not contain "N"
        with open(str(os.path.join(snakemake.output[0], file)), "w") as output_handle:
            SeqIO.write(x, output_handle, "fasta")
