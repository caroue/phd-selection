### use a python JSON parser such that you can call the key name directly:
import json
import re

pos=[]

for i in snakemake.input:
#    f=open(snakemake.input[i], "r")
    print(i)
    sa=re.sub('absrel/', '', i)
    sample=re.sub('_absrel.json', '', sa)
    print(sample)
    f=open(i, "r")
    data=json.load(f)
    dict=data['branch attributes']['0']
    for key in dict:
        pos.append(sample + '\t' + str(key) + '\t' + str(dict[key]['Uncorrected P-value']) + '\t' + str(dict[key]['Corrected P-value']))

with open(snakemake.output[0], mode='wt', encoding='utf-8') as myfile:
    myfile.write('\n'.join(pos))
