# Compute stats comparing ref graph with other graph, requires kmer_counts file

import matplotlib
matplotlib.use('Agg') # so we don't display figures (segfaults in snakemake)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#compute_completeness
# completeness = fraction of k-mers in reference graph that are also in other graph
# also report absolute numbers
def compute_completeness(df):
    nb_kmers = df.shape[0] # in reference
    nb_kmers_other = (df["bf_counts"] > 0).sum()
    fraction = nb_kmers_other/nb_kmers

        
    print(f"Fraction: {fraction}")
    print(f"nb_kmers_other: {nb_kmers_other}")
    print(f"nb_kmers in ref: {nb_kmers}")

# histogram over counts in reference graph and non-ref graph
def create_count_histogram(df, save_filename):
    fig, ax = plt.subplots()
    fig.tight_layout()
    ax.set_title("Histogram over k-mer counts")
    sns.histplot(df.melt(var_name = 'Graph'), x='value', hue='Graph', ax=ax, multiple='dodge', shrink=.75, bins=np.arange(0, 50)-0.5)
    ax.set_yscale('log')
    ax.set_xlabel('k-mer count')
    ax.set_ylabel('# k-mers with count')
    fig.savefig(save_filename, bbox_inches='tight')

def main(kmer_counts_filename, completeness_out_file, out_image_file):
    from contextlib import redirect_stdout
    import pandas as pd

    df = pd.read_csv(kmer_counts_filename, sep = " ", names = ['ref_counts', 'bf_counts'])

    with open(completeness_out_file, "w") as f:
        with redirect_stdout(f):
            compute_completeness(df)

    create_count_histogram(df, out_image_file)

if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])

