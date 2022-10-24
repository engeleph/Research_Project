import toml
from pathlib import Path
import time
import snakemake
import sys
import shutil

from snakemake_helpers import *

"""
This script keeps the bloom filter graph (BFG) updated by looking for new data (.fa.gz) and building a new BFG.

It obtains the current BFG from the config file, builds a new BFG, and updates the config file accordingly.
This means that this script can be resumed from an existing metagraph.

If the current BFG does not exist, it creates an empty one.

It simulates arriving fastas if the `fasta_dir` is initially empty.

Run with:
        rm -rf fasta_data outputs
        python keep_graph_updated.py
        # you can Ctrl-C interrupt it and then restart it, so it resumes
        python keep_graph_updated.py
        # then compute statistics, see readme

        # merge into one file, Mac OS X: pdfunite
        pdfunite outputs/pipeline_*/counts_hist.pdf outputs/merged_counts_hist.pdf
        open outputs/merged_counts_hist.pdf
"""

############ Python helpers ############

# read bf graph location
def parse_graph_path(config_path):
    # print(f"Parsing metagraph reference from file '{config_path}'")
    data = toml.load(config_path)
    p = Path(data["conditions"]["reference"])
    assert(p.name.endswith(".bfdbg"))
    return p

# write bf graph location
def write_graph_path(config_path, new_graph_filename):
    data = toml.load(config_path)
    data["conditions"]["reference"] = str(Path(new_graph_filename).resolve()) # otherwise PosixPath(filename) will be written to config file
    with open(config_path, "w") as f:
        toml.dump(data, f)

# create live config by copying configName.suffix to new file named configName_live.suffix
def create_live_path(config_path):
    config_path = Path(config_path)
    new_config_path = config_path.parent / (config_path.stem + "_live" + config_path.suffix)
    shutil.copy(config_path, new_config_path)
    return new_config_path

def is_empty_dir(direc):
    it = Path(direc).iterdir()
    for _ in it:
        return False
    return True

# first splits big fasta file into small files, then moves them into directory at fixed time intervals, optionally with initial wait
def simulate_arriving_fastas(big_fasta_file, fasta_dir, n_parts, time_interval, initial_wait=1):
    import os
    import tempfile
    import threading

    def f():
        time.sleep(initial_wait) # not required
                                    
        temp_fasta_dir = Path(tempfile.mkdtemp())
        exit_code = os.system(f"""seqkit split2 {big_fasta_file} -p {n_parts} -e .gz -o split -O {str(temp_fasta_dir)} -f""")
        assert exit_code == 0, "seqkit command unsuccessful"

        for file in sorted(temp_fasta_dir.iterdir()):
            source_file = temp_fasta_dir / file.name
            target_file = fasta_dir / file.name
            if target_file.exists():
                # useful when the program is resumed
                print(f"Skipping move since target file '{target_file}' already exists", flush=True)
                source_file.unlink() # so that directory can be deleted at the end
                continue

            print(f"Moving file {source_file} -> {target_file}")
            shutil.move(source_file, target_file)
            time.sleep(time_interval)
        temp_fasta_dir.rmdir() # should be empty

    # f()
    t = threading.Thread(target=f, daemon=True) # set daemon so that main program does not wait for this one
    t.start()


############ params ############
# readfish_config_path = create_live_path("readfish_config.toml")
readfish_config_path = "readfish_config_live.toml"
snake_config_path = "snakemake_config.yaml"
create_graph_if_inexistent = True
delete_intermediate_outputs = True # whether to delete intermediate outputs after building the merged graph
delete_previous_graph = False # whether to delete the graph once a merged graph is obtained
simulate_fastas = True # whether to simulate fasta files arriving in fasta dir
simulation_config = {
    "time_interval": 10, # time interval at which fasta will be deposited
    "n_parts": 6, # number of fasta files to split into
}
new_fasta_check_interval = 3 # interval in seconds at which to check for new fasta files
threads = None # None: get max available


############ vars ############
snake_config = snakemake.load_configfile(snake_config_path)
fasta_dir = Path(snake_config["fasta_dir"])

existing_graph_path = parse_graph_path(readfish_config_path)
existing_processed_fastas_path = Path(get_processed_files_filename(existing_graph_path))

snake_config['max_threads'] = threads if threads else snakemake.available_cpu_count()

############ init pipeline, create empty graph if necessary; if it exists, check that processed_fastas files also exists ############
if create_graph_if_inexistent and not existing_graph_path.exists():
    # snakemake --cores 4 --snakefile pipelines/pipeline_init.smk --config existing_bf_graph="outputs/pipeline_zero/merged_zero.bfdbg"
    print(f"Graph '{existing_graph_path}' does not exist, so creating it")
    existing_graph_path.parent.mkdir(parents=True, exist_ok=True)
    assert snakemake.snakemake("init_pipeline.smk", config={**snake_config, "existing_bf_graph": existing_graph_path}), "Could not initialize pipeline"

if not existing_graph_path.exists() or not existing_processed_fastas_path.exists():
    print(f"Graph '{existing_graph_path}' or processed_files '{existing_processed_fastas_path}' does not exist.")
    sys.exit(1)

processed_fastas = set(str(Path(p).resolve()) for p in read_list(existing_processed_fastas_path))
print("Already processed fastas: {}.".format(", ".join(processed_fastas)))
existing_processed_fastas_path = None # will not be updated in loop, so set to None

fasta_dir.mkdir(parents=True, exist_ok=True)
if simulate_fastas:
    if not is_empty_dir(fasta_dir):
        print("Warning: Simulation of fastas into non-empty directory, probably because simulation is resumed.")
    simulate_arriving_fastas(snake_config["ref_data"], fasta_dir=fasta_dir, **simulation_config)


#todo3: remove todos
#todo3: remove config.yaml
# todo3: remove
# Nanopore device setup
# setup device: https://community.nanoporetech.com/to/minion-start
# https://community.nanoporetech.com/support/articles/why-do-i-need-a-ssd-can-i
# MinKNOW needs to be installed on a SSD with at least 1TB of space
# https://community.nanoporetech.com/support/articles/is-there-a-way-to-view-my-mk1c-sequencing-run-on-a-bigger-screen
# download MinION Software for Mk1b
# MinKNOW GUI cannot run Mk1b devices
# open MinKNOW, not MinKNOW-UI (only has UI, no basecalling, which is done on device)
# https://community.nanoporetech.com/posts/software-patch-release-21-9359: You can start an experiment with a CTC or used flow cell for the purposes of this test.
# todo3: configure GPU: https://community.nanoporetech.com/docs/prepare/library_prep_protocols/experiment-companion-minknow/v/mke_1013_v1_revcl_11apr2016/installing-gpu-version-of-guppy-with-minknow-for-minion
# todo3: basecalling several times (twice): fast basecalling or not?, no fast5 output (just fastq)
# try setup on cluster

# pass and fail folders
# todo3: allow fasta and fastq and gz etc

# todo3: try simulating a run

while True:
    new_fasta_files = get_fasta_files(fasta_dir) - processed_fastas
            
    if len(new_fasta_files) == 0:
        time.sleep(new_fasta_check_interval)
        print("Did not discover new fasta (.fa.gz) files")
        # print(f"Current content: {get_fasta_files(fasta_dir)}")
        continue
                                                    
    print("Discovered new fasta files: {}".format(", ".join(new_fasta_files)))
                                                            
    unique_id = time.time()
                                                                    
    # new paths
    output_dir = get_output_path(snake_config["base_output_dir"], unique_id)
    graph_path = get_merged_graph_path(output_dir, unique_id)
                                                                                    
    print(f"Running pipeline with id {unique_id} and output_dir '{output_dir}'")
    # write files to process to output_dir, so that the pipeline picks them up
    files_processing_filename = get_files_processing_filename(output_dir, unique_id)
    files_processing_filename.parent.mkdir(parents=True, exist_ok=True)
    write_list(files_processing_filename, new_fasta_files)
                    
    # snakemake --cores 4 --snakefile pipelines/merge_pipeline.smk --config existing_bf_graph="outputs/pipeline_zero/merged_zero.bfdbg" unique_id="first"
    assert snakemake.snakemake("pipelines/merge_pipeline.smk", config={**snake_config, "existing_bf_graph": existing_graph_path, "unique_id": unique_id}), "Could not merge graph"

    # update toml (before deleting previous graph; otherwise, if the program crashes right after deleting old graph, but before modifying the toml, the toml is no longer valid)
    print(f"Updating toml {readfish_config_path}")
    write_graph_path(readfish_config_path, graph_path)

    if delete_intermediate_outputs or delete_previous_graph:
        # snakemake --cores 4 --snakefile pipelines/merge_pipeline.smk --config existing_bf_graph="outputs/pipeline_zero/merged_zero.bfdbg" unique_id="first" -- clean_intermediate
        print("Deleting intermediate graph files")
        assert snakemake.snakemake("pipelines/merge_pipeline.smk", targets=["clean_intermediate"], config={**snake_config, "existing_bf_graph": existing_graph_path, "unique_id": unique_id}), "Could not clean intermediate graph dirs"
    if delete_previous_graph:
        existing_graph_path.unlink()
        Path(get_processed_files_filename(existing_graph_path)).unlink()
        try:
            existing_graph_path.parent.rmdir()
            print(f"Deleted previous graph's dir '{existing_graph_path.parent}'")
        except:
            # there were more files (happens when the initial graph lives in a directory not created by this program)
            print(f"Directory '{existing_graph_path.parent}' not empty, so not deleting it")

    # update graph path and processed fastas
    existing_graph_path = graph_path
    processed_fastas = processed_fastas.union(new_fasta_files)


