# insert the list name (OGs_x.txt) of the OG set to analyze
with open("data/OGs_full_list_removed_1_HOG.txt", "r") as file:
    SAMPLES = [sample for sample in file.read().split("\n") if len(sample) > 0]

rule all:
    input:
        "OGs_number_of_branches_under_pos_selection.tsv",
        "OGs_pos_selection_pvalues.tsv"


rule extract_protein_names:
    input:
        "data/Conserved_dataset_matrix.tsv"
    output:
        "proteins_per_OG/Proteins_{sample}.txt"
    shell:
        """
        grep "{wildcards.sample}" {input} | awk -F '\\t' '{{print $2}}' > {output}
        """

#rule protein_to_transcript_ID:
#    input:
#        "proteins_per_OG/Proteins_{sample}.txt"
#    output:
#        "transcripts_per_OG/Transcripts_{sample}.txt"
#    shell:
#        """
#        sed "s/BRADI.*\./Bd_transcript_/g" {input} | sed 's/HORVU/Hv_transcript_HORVU/g' | sed 's/SORBI.*\./Sb_transcript_/g' | sed 's/SETIT.*\./Si_transcript_/g' | sed '/Zm/ s/P/T/g'| sed 's/Zm/transcript_Zm/g' | sed 's/Traes/transcript_Traes/g' > {output}
#        """
# use double quotes for the sed command

rule extract_CDS:
    input:
        cds="data/All_species_CDS_longestiso_gene_name.fa",
        proteins="proteins_per_OG/Proteins_{sample}.txt"
    output:
        "CDS/CDS_{sample}.fasta"
    conda:
        "envs/fasomerecords.yaml"
    shell:
        """
        faSomeRecords {input.cds} {input.proteins} {output}
        """

### checkpoint: check CDSs for "N" in the sequence, remove files which contain "N" in the sequence
### and only move to rule pre_msa the retained files
#We pretend that the number of clusters is unknown beforehand. Hence, the checkpoint only defines an output directory
checkpoint clean_CDS:
    input:
        expand("CDS/CDS_{sample}.fasta", sample=SAMPLES)
    output:
        directory("clean_CDS")
    conda:
        "envs/Biopython.yaml"
    script:
        "scripts/clean_CDS_2.py"
#        script which opens each CDS_{sample}.fasta, reads the CDS, outputs only if no "N" is found

rule pre_msa:
    input:
        "clean_CDS/CDS_{simple}.fasta"
    output:
        AA="pre-msa/{simple}_AA.fasta",
        NT="pre-msa/{simple}_NT.fasta",
        fil="pre-msa/{simple}_filtered.json",
        cop="pre-msa/{simple}_copies.json"
    conda:
        "envs/hyphy.yaml"
    shell:
#        "hyphy /scratch1/cropbio/uebermuth/software/hyphy-analyses/codon-msa/pre-msa.bf --input {input} --protein {output.AA} --rna {output.NT} --filter {output.fil} --copies {output.cop}"
### try with --E 0.05 option to allow alignment of sequences with low homology 
        "hyphy /scratch1/cropbio/uebermuth/software/hyphy-analyses/codon-msa/pre-msa.bf --E 0.05 --input {input} --protein {output.AA} --rna {output.NT} --filter {output.fil} --copies {output.cop}"

rule msa:
    input:
        "pre-msa/{simple}_AA.fasta"
    output:
        "msa/{simple}_AA_al.fasta"
    conda:
        "envs/muscle.yaml"
    shell:
        "muscle -align {input} -output {output}"

rule post_msa:
    input:
        protal="msa/{simple}_AA_al.fasta",
        NT="pre-msa/{simple}_NT.fasta"
    output:
        "post-msa/{simple}_NT_al.fasta"
    conda:
        "envs/hyphy.yaml"
    shell:
        "hyphy /scratch1/cropbio/uebermuth/software/hyphy-analyses/codon-msa/post-msa.bf --protein-msa {input.protal} --nucleotide-sequences {input.NT} --output {output} --compress No"

def aggregate_input(wildcards):
    ck_output = checkpoints.clean_CDS.get(**wildcards).output[0]
    file_names= expand("post-msa/{simple}_NT_al.fasta",
                       simple = glob_wildcards(os.path.join(ck_output, "CDS_{simple}.fasta")).simple)
    return file_names

### checkpoint: check nucleotide msa files for "?" in the sequence, remove files which contain "?" in the sequence
### and only move to rule tree the retained files
#We pretend that the number of clusters is unknown beforehand. Hence, the checkpoint only defines an output directory
checkpoint clean_msa:
    input:
        aggregate_input
    output:
        directory("clean_msa")
    conda:
        "envs/Biopython.yaml"
    script:
        "scripts/clean_CDS_3.py"
#        script which opens each {sample}_NT_al.fasta, reads the nucleotide alignment, outputs only if no "N" nor "?" is found

rule tree:
    input:
        "clean_msa/{sumple}_NT_al.fasta"
    output:
        "tree/{sumple}_NT_al.nwk"
    conda:
        "envs/fasttree.yaml"
    shell:
        "fasttree -nt {input} > {output}"

rule absrel:
    input:
        al="post-msa/{sumple}_NT_al.fasta",
        tree="tree/{sumple}_NT_al.nwk"
    output:
        json="absrel/{sumple}_absrel.json",
        absrelout="absrel/{sumple}_absrel.out"
    conda:
        "envs/hyphy.yaml"
    threads:
        8
    shell:
        "hyphy absrel --alignment {input.al} --tree {input.tree} --output {output.json}  > {output.absrelout}"

### the following rules require input with wildcards generates by checkpoints and output only textfiles (which are finally queried by rule_all. 
### Hence, we need an input function here which defines the wildcards created by checkpoint clean_CDS

def aggregate_input_2(wildcards):
    ck_output = checkpoints.clean_msa.get(**wildcards).output[0]
    file_names= expand("absrel/{sumple}_absrel.json",
                       sumple = glob_wildcards(os.path.join(ck_output, "{sumple}_NT_al.fasta")).sumple)
    return file_names



rule create_table_numbers:
    input:
        aggregate_input_2
    output:
        "OGs_number_of_branches_under_pos_selection.tsv"
    script:
        "scripts/Create_OGs_number_of_pos_selection.py"

rule create_table_pvalues:
    input:
        aggregate_input_2
    output:
        "OGs_pos_selection_pvalues.tsv"
    script:
        "scripts/Create_OGs_pos_selection_pvalue.py"
