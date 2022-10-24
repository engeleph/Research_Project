
# Initializes the merged graph creation pipeline

from snakemake_helpers import *


##################### Parameters #####################

configfile: "snakemake_config.yaml" # todo3: remove

# passed via command line
existing_bf_graph = config["existing_bf_graph"]

K = config["K"]
bloom_num_hash_fcns = config["bloom_num_hash_fcns"]
bloom_filter_size_exp = config["bloom_filter_size_exp"]
metagraph_exec = config.get("metagraph_exec", "metagraph")


##################### Define vars #####################

existing_processed_files_filename = get_processed_files_filename(existing_bf_graph)


##################### Rules #####################

# standalone rule to build an empty bloom filter graph, if desired
rule build_empty_bf:
    params:
        existing_bf_graph_no_ext = remove_graph_extension(existing_bf_graph)
    message: "Putting empty graph to location '{existing_bf_graph}'"
    shell: 
        """
        empty_file=empty_file.fasta.gz # metagraph requires file extension, so process substitution '<(true)' does not work
        trap 'rm -f -- "$empty_file"' EXIT
        touch "$empty_file"
        {metagraph_exec} build -v -p {threads} -k {K} --count-kmers -o {params.existing_bf_graph_no_ext} "$empty_file"
        true > {existing_processed_files_filename:q}
        """
