# Merges new fasta data into existing bloom graph

from snakemake_helpers import *


##################### Parameters #####################

# load config
configfile: "snakemake_config.yaml" # command line can overwrite these parameters

# passed via command line
unique_id = config["unique_id"] # using 'unique_id = time.time()' does not work because code referring to it in the "run:" or "shell:" will re-evaluate the python code, so also the unique_id -> must be passed in config!!
existing_bf_graph = config["existing_bf_graph"]

K = config["K"] # for succinct graph building
fasta_dir = config["fasta_dir"]
base_output_dir = config.get("base_output_dir", "outputs")
metagraph_exec = config.get("metagraph_exec", "metagraph")


##################### Define vars #####################

# output directories

output_dir = get_output_path(base_output_dir, unique_id)
output_dir.mkdir(parents=True, exist_ok=True)
log_dir = output_dir / "logs"
benchmark_dir = output_dir / "benchmarks"

# tracking processed files and files in process
files_processing_filename = get_files_processing_filename(output_dir, unique_id)

# intermediate and final files
# put into separate directories for easy cleaning (because each graph file may have an associated weights file)
succinct_graph = output_dir / "uncleaned_succinct" / f"graph.dbg"
cleaned_fasta = output_dir / "cleaned_fasta" / "cleaned.fasta.gz"
# existing graph
assert(Path(existing_bf_graph).exists())
existing_processed_files_filename = get_processed_files_filename(existing_bf_graph)
# new merged graph
merged_graph = get_merged_graph_path(output_dir, unique_id)
merged_processed_files_filename = get_processed_files_filename(merged_graph)


##################### Rules #####################

# first rule = default rule
# triggers pipeline
rule all:
    input: merged_graph


##################### Standalone Rules #####################

# standalone rule to clean everything output by this pipeline (only makes sense if a unique_id is provided)
rule clean:
    # do not delete existing since the graph may not have been created by this pipeline
    message: "Cleaning output dir '{output_dir:q}'"
    shell: "rm -rf -- {output_dir:q} || true"

rule clean_intermediate:
    message: "Cleaning intermediate files: '{succinct_graph:q}' succinct dir, '{cleaned_fasta:q}' cleaned fasta dir, '{benchmark_dir:q}' benchmark dir, '{log_dir:q}' log dir"
    shell: """rm -rf -- "$(dirname {succinct_graph:q})" "$(dirname {cleaned_fasta:q})" {benchmark_dir:q} {log_dir:q} || true"""

##################### Rules to merge in new data into an existing graph #####################

# save files to process to a file, so that the pipeline can be resumed if it fails at an intermediate stage
if files_processing_filename.exists():
    # pipeline was previously run, so read files that were processed in that pipeline
    new_files = read_list(files_processing_filename)
    print("Read new files {}".format(', '.join(new_files)))
else:
    new_files = get_new_fa_gz_files(fasta_dir, existing_processed_files_filename)
    print("Computed new files {}".format(', '.join(new_files)))
if len(new_files) == 0:
    print("No new files to process")
    sys.exit(1)

rule write_files_to_process:
    input: existing_processed_files_filename
    output: files_processing_filename
    message: "Getting new files to process, reading from '{existing_processed_files_filename}', writing difference to '{files_processing_filename}'"
    run: write_list(files_processing_filename, new_files)

rule build_succinct_dbg:
    input: files_processing_filename
    output: succinct_graph
    params: 
        succinct_graph_no_ext = remove_graph_extension(succinct_graph), 
        fasta_files = lambda wildcards, output, input: read_list(str(input))
    message: "Building succinct graph '{succinct_graph}' from files: {input:q}"
    benchmark: benchmark_dir / "build_succinct_dbg.tsv"
    shell: "{metagraph_exec} build -v -p {threads} -k {K} -o {params.succinct_graph_no_ext} --count-kmers {params.fasta_files:q}"

# cleans and writes to fasta
rule succinct_to_fasta:
    input: succinct_graph
    output: cleaned_fasta
    params:
        cleaned_fasta_without_ext = remove_fasta_gz_ext(cleaned_fasta)
    message: "Converting succinct graph {succinct_graph} to cleaned fasta file {cleaned_fasta} with weights"
    benchmark: benchmark_dir / "succinct_to_fasta.tsv"
    shell: "{metagraph_exec} clean -v -p {threads} --to-fasta --prune-tips $((2*{K})) --prune-unitigs 0 --fallback 2 -o {params.cleaned_fasta_without_ext} {succinct_graph}"

#todo3
# # helper rule to add empty processed files if inexistent
# rule add_processed_files_if_inexistent:
#     input: existing_bf_graph # graph should exist
#     message: "Add processed_files for existing graph at location '{existing_processed_files_filename}'"
#     output: existing_processed_files_filename
#     shell: "true > {existing_processed_files_filename:q}"

rule extend_bloom_with_fasta:
    input: existing_bf_graph, existing_processed_files_filename, cleaned_fasta
    output: merged_graph, merged_processed_files_filename
    params:
        merged_graph_no_ext = remove_graph_extension(merged_graph),
    message: "Extending existing graph '{existing_bf_graph}' with cleaned fasta file '{cleaned_fasta}' to produce '{merged_graph}'"
    benchmark: benchmark_dir / "extend_bloom_with_fasta.tsv"
    shell:
        '''
        {metagraph_exec} extend -v -p {threads} -i {existing_bf_graph:q} -o {params.merged_graph_no_ext} {cleaned_fasta:q}
        cat {existing_processed_files_filename:q} {files_processing_filename:q} > {merged_processed_files_filename:q}
        '''

