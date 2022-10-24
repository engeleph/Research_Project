##################### Python Helpers #####################

from pathlib import Path 

# extract unique_id from metagraph file
# filename of the form 'graphname_<unique_id>.graphext'
def extract_unique_id(filename):
        return Path(filename).stem.split("_")[-1]

# remove last extension, e.g. '.dbg'
def remove_graph_extension(filename):
    p = Path(filename)
    return p.parent / p.stem
                
def remove_fasta_gz_ext(filename):
    filename = str(filename)
    ext = ".fasta.gz"
    assert filename.endswith(ext)
    return filename[:-len(ext)]

def read_list(filename):
    with open(filename, "r") as f:
        return list(filter(lambda x: x != "", map(lambda x: x.strip(), f.readlines())))

def write_list(filename, lst):
    with open(filename, "w") as f:
        print("\n".join(lst), file=f)

# returns set of absolute paths
def get_fasta_files(fasta_dir):
    return set(str(Path(p).resolve()) for p in Path(fasta_dir).iterdir() if p.name.endswith(".fa.gz"))

# get new .fa.gz files in directory that were not already processed (written in file)
def get_new_fa_gz_files(fasta_dir, processed_files_filename):
    return get_fasta_files(fasta_dir) - set(read_list(processed_files_filename))

processed_files_format = "processed_files_{}.txt" # format string for filename of already processed files
# the processed_files lives in the same directory as the graph
def get_processed_files_filename(existing_graph_filename):
    return Path(existing_graph_filename).parent / processed_files_format.format(extract_unique_id(existing_graph_filename))
# files that are currently being processed
get_files_processing_filename = lambda output_dir, unique_id: output_dir / f"files_processing_{unique_id}.txt"

# directory paths
get_output_path = lambda base_output_dir, unique_id: Path(base_output_dir) / f"pipeline_{unique_id}"
get_merged_graph_path = lambda output_dir, unique_id: output_dir / f"merged_{unique_id}.bfdbg"
