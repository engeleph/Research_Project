# splits big fasta file into small chunks to test the other pipeline

##################### Parameters #####################

# load config
configfile: "snakemake_config.yaml" # command line can overwrite these parameters

big_fasta_file = config["ref_data"]
splitted_sequence = config["splitted_data"]
parts=list(range(1,3751))


##################### Rules #####################

# standalone rule to put fastas into fasta directory
rule split ref:
    input: big_fasta_file
    output:
        splitted_sequence
    message: "Splits fasta file {input} into reads into file '{output}'"
    shell: "seqkit sliding {input} -s 990 -W 1000 -o {output}"
