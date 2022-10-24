# Given a directory of merged graphs, it computes the completeness of the graph with respect to a reference graph
# The reference graph is built as a succinct DBG, if it does not exist.

import glob
import os

from snakemake_helpers import *
import compute_graph_stats

##################### Parameters #####################

# load config
configfile: "snakemake_config.yaml" # command line can overwrite these parameters

base_output_dir = config["base_output_dir"] # base directory where merged graphs are in
ref_graph_path = config["ref_graph_path"] # where the reference graph is, or where it should be built to
ref_data = config["splitted_data"]
output_dir = config["output_dir"]
K = config["K"] # for succinct graph building, >= K used for merged graphs (ideally equal)

metagraph_exec = config.get("metagraph_exec", "metagraph")
count_distribution_exec = config.get("count_distribution_exec", "count_distribution")


##################### Define vars #####################

# output directories for ref graph construction
log_dir = "{}/logs".format(base_output_dir)
benchmark_dir = "{}/benchmarks".format(base_output_dir)

if ref_data is None:
    assert ref_graph_path.exists()
#else:
#    output_dir.mkdir(exist_ok=True)
    

##################### Rules #####################

# targets cannot include wildcards, so we glob first
merged_graph_dirs = glob.glob(str(get_output_path(base_output_dir, "*")))
output_files = expand(
    [str(Path(p) / "{output_filename}") for p in merged_graph_dirs],
    output_filename = ["counts.txt", "completeness.txt", "counts_hist.pdf"]
)

# default rule
rule all:
    input:
        completeness = os.path.join(base_output_dir,"{unique_id}","completeness.txt"),
        counts_hist = os.path.join(base_output_dir, "{unique_id}","counts_hist.pdf"),

rule clean:
    message: "Cleaning statistics in subdirectories of {merged_graph_dirs:q}"
    shell: "set -x; rm -rf -- {output_files:q}"

# build ref graph
rule build_ref_graph:
    output: ref_graph_path
    params:
        ref_graph_no_ext = remove_graph_extension(ref_graph_path),
    message: "Building succinct ref graph '{ref_graph_path}' from ref data '{ref_data}'"
    benchmark: os.path.join(benchmark_dir, "build_ref_graph.tsv")
    shell: "{metagraph_exec} build -v -p {threads} -k {K} --count-kmers -o {params.ref_graph_no_ext:q} {ref_data:q}"

# compute counts file
rule compute_counts:
    input: 
        ref_graph = ref_graph_path, 
        #merged_graph = get_output_path(base_output_dir, "{unique_id}") / "merged_{unique_id}.bfdbg",
        merged_graph = os.path.join(base_output_dir,"{unique_id}","merged_{unique_id}.bfdbg")
    output: 
        counts = os.path.join(base_output_dir, "{unique_id}","counts.txt"),
    message: "Computing counts for graph '{input.merged_graph}'"
    benchmark: os.path.join(base_output_dir, "{unique_id","compute_counts.tsv")
    shell: "{count_distribution_exec} {input.ref_graph:q} {input.merged_graph:q} > {output.counts:q}"

# compute completion rate and histogram over counts vs ref counts
rule compute_stats:
    input: 
        ref_graph = ref_graph_path, 
        counts = os.path.join(base_output_dir,"{unique_id}","counts.txt"),
    output: 
        completeness = os.path.join(base_output_dir,"{unique_id}","completeness.txt"),
        counts_hist = os.path.join(base_output_dir, "{unique_id}","counts_hist.pdf"),
    params:
        output_dir = lambda wildcards: os.path.join(base_output_dir, f"{wildcards.unique_id}")
    message: "Computing stats for directory '{params.output_dir}'"
   # benchmark: os.path.join(base_output_dir, "{unique_id}", "compute_stats.tsv")
    run: compute_graph_stats.main(input.counts, output.completeness, output.counts_hist)
