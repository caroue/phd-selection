# phd-selection
works on snakemake version 7.0.0

# data/: 
softlink here the "Conserved_dataset_matrix.tsv", which contains the OG in the first column and the corresponding protein names in the second column 
softlink here the fasta file which contains the CDS for all proteins of all species. NOTE: the fasta headers need to contain the same protein ID as the "Conserved_dataset_matrix.tsv". Must not contain ":"
 softlink here a list of all OGs you want to be analyzed. e.g. "OGs_in_conserved_dataset.txt"

# envs/:
specifies the envs for each programme which is installed by conda.

# software NOT installed with conda:
pre- and post-msa.bf scripts (codon-msa) from hyphy developer "spond" in: https://github.com/veg/hyphy-analyses.git
/scratch1/cropbio/uebermuth/software/hyphy-analyses/codon-msa/pre-msa.bf and /scratch1/cropbio/uebermuth/software/hyphy-analyses/codon-msa/post-msa.bf

# start the workflow:
write the name of the list containing the OGs you want to be analyzed into the first line of the Snakefile
start the workflow with as many cores as you want and use the --use-conda flag for the necessary envs to be installed on the run
snakemake --cores 90 --use-conda
