# splits big fasta file into small chunks to test the other pipeline

##################### Parameters #####################

# load config
configfile: "snakemake_config.yaml" # command line can overwrite these parameters

big_fasta_file = config["splitted_data"]
fasta_dir = config["fasta_dir"]


 ##################### Rules #####################

# standalone rule to put fastas into fasta directory
rule populate_fasta_dir:
    input: big_fasta_file
    output: fasta_dir
    message: "Splitting fasta file {big_fasta_file} into parts into directory '{fasta_dir}'"
    shell: 
        '''
        seqkit split2 {big_fasta_file:q} -s 3750 -o split -O {fasta_dir:q} -f
        gzip {output}/*
        '''

