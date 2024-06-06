# this writes only files which do not contain "N" nor "X" nor "?" in any record of the fasta sequence imported with seqio and which do not contain any
# record of the fasta sequence imported with seqio whichs length is not a multiple of 3

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
    seqcumlen=0
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
#    print(cumlen)
        sequencelength = len(str(sequence))
#        print(sequencelength)
        k = 3
        if (sequencelength % k == 0): #"is an exact multiple of k, no remainder from the division"
            seqcumlen = seqcumlen
        else:
            seqcumlen = seqcumlen + 1
    print("seqcumlen is", seqcumlen)
#select only input files which do not contain "N" nor "X" nor "?" or there length is NOT a multiple of 3 in any of their record's sequence
    if ((cumlen == 0) and (seqcumlen == 0)):
#write it to output file only if it does not contain "N" nor "X" nor "?" nor NOT a multiple of 3
        with open(str(os.path.join(snakemake.output[0], file)), "w") as output_handle:
            SeqIO.write(x, output_handle, "fasta")
