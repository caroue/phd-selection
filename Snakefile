with open("data/OGs.txt", "r") as file:
    SAMPLES = [sample for sample in file.read().split("\n") if len(sample) >0]

rule all:
    input:
        "OGs_number_of_branches_under_pos_selection.tsv"


rule extract_protein_names:
    input:
        "data/Conserved_dataset_matrix.tsv"
    output:
        "proteins_per_OG/Proteins_{sample}.txt"
    shell:
        """
        grep "{wildcards.sample}" {input} | awk -F '\\t' '{{print $2}}' > {output}
        """

rule protein_to_transcript_ID:
    input:
        "proteins_per_OG/Proteins_{sample}.txt"
#        expand("proteins_per_OG/Proteins_{sample}.txt", sample=SAMPLES)
    output:
        "transcripts_per_OG/Transcripts_{sample}.txt"
    shell:
        """
        sed "s/BRADI.*\./Bd_transcript_/g" {input} | sed 's/HORVU/Hv_transcript_HORVU/g' | sed 's/SORBI.*\./Sb_transcript_/g' | sed 's/SETIT.*\./Si_transcript_/g' | sed '/Zm/ s/P/T/g'| sed 's/Zm/transcript_Zm/g' | sed 's/Traes/transcript_Traes/g' > {output}
#        sed "s/BRADI.*\./Bd_transcript:/g" {input} | sed 's/HORVU/Hv_transcript:HORVU/g' | sed 's/SORBI.*\./Sb_transcript:/g' | sed 's/SETIT.*\./Si_transcript:/g' | sed '/Zm/ s/P/T/g'| sed 's/Zm/transcript:Zm/g' | sed 's/Traes/transcript:Traes/g' > {output}
        """
# use double quotes for the sed command
#        "scripts/protein_to_transcript_IDs.sh < {input} > {output}"
    
rule extract_CDS:
    input:
        cds="data/All_species_CDS_edited_2.fa",
        transcripts="transcripts_per_OG/Transcripts_{sample}.txt"
    output:
        "CDS/CDS_{sample}.fasta"
    conda:
        "envs/fasomerecords.yaml"
    shell:
        """
        faSomeRecords {input.cds} {input.transcripts} {output}
        """

#alignment funktionierte nicht, weil die transkripte im fasta header ein ":" hatten. Habe ich in der rule protein_to_transcript_ID und im All_species_CDS_edited_2.fa auf "_" geändert
#alignment fasta hat long lines, die nicht gebrochen sind!??
rule alignment:
    input:
        "CDS/CDS_{sample}.fasta"
    output:
        nt="alignment/{sample}_NT_al.fasta",
        aa="alignment/{sample}_AA_al.fasta"
    shell:
        "java -jar -Xmx600m /scratch1/cropbio/uebermuth/software/macse_v2.06.jar -prog alignSequences "
        "-seq {input} -out_NT {output.nt} -out_AA {output.aa}"

rule remove_stop_codon_from_NT_AL:
    input:
        "alignment/{sample}_NT_al.fasta"
    output:
        "alignment/nostop/{sample}_NT_al_nostop.fasta"
    conda:
#        "envs/biopython_and_snakemake.yaml"
        "envs/Biopython.yaml"
    script:
        "scripts/replace_stop_codons.py"

rule tree:
    input:
        "alignment/nostop/{sample}_NT_al_nostop.fasta"
    output:
        "tree/{sample}_NT_al_nostop.nwk"
    conda:
        "envs/fasttree.yaml"
    shell:
        "fasttree -nt {input} > {output}"

rule absrel:
    input:
        al="alignment/nostop/{sample}_NT_al_nostop.fasta",
        tree="tree/{sample}_NT_al_nostop.nwk"
    output:
        json="absrel/{sample}_absrel.json",
        absrelout="absrel/{sample}_absrel.out"
    conda:
        "envs/hyphy.yaml"
    shell:
        "hyphy absrel --alignment {input.al} --tree {input.tree} --output {output.json}  > {output.absrelout}"

rule create_table:
    input:
        expand("absrel/{bla}_absrel.json", bla=SAMPLES)
    output:
        "OGs_number_of_branches_under_pos_selection.tsv"
    script:
        "scripts/Create_OGs_number_of_pos_selection.py"
