#imports

import gseapy
import pandas as pd
import json

#defining the function that executes the enrichment analysis

def run_enrichr(
    gene_list, gene_sets, comparison, n_genes, go_lookup, desc_col, direction=None
):
    enr2 = gseapy.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        background=n_genes,
        outdir="enrichr_test/",
        cutoff=0.05,
        no_plot=True,
        verbose=False,
    )
    results = (
        enr2.results.rename(columns={"Adjusted P-value": "padj"})
        .set_index("Term")
        .sort_values("P-value", ascending=True)
    )
    results["Gene_set"] = go_lookup.loc[results.index, desc_col]
    if direction is not None:
        results["Direction"] = direction
    significant = results.query("padj < 0.05")
    return (results, significant)

#Calling the function on specific DE genes differently for up and down regulated genes

if __name__ == "__main__":
    lookup = pd.read_csv('go_lookup_table.csv', index_col=0)
    de_results = pd.read_csv('new_analysis/deseq_results/results_12_18_A_C.csv', index_col='uniprothit')
    comparison = '12 AZ v/s 18 Control'
    up_genes = de_results.query("padj < 0.05 and log2FoldChange > 1")
    down_genes = de_results.query("padj < 0.05 and log2FoldChange < -1")
    
    with open('go_gene_list.json', "r") as f:
        gene_sets = json.load(f)
    
    
    up_genes_list = [x for x in up_genes.index.values if not pd.isna(x)]
    down_genes_list = [x for x in down_genes.index.values if not pd.isna(x)]
    

    up_results = run_enrichr(
        up_genes_list,
        gene_sets,
        comparison,
        de_results.shape[0],
        lookup,
        'desc',
        direction="+",
    )
    down_results = run_enrichr(
        down_genes_list,
        gene_sets,
        comparison,
        de_results.shape[0],
        lookup,
        'desc',
        direction="-",
    )
    for i, out_fn in enumerate(
        ["new_analysis/enrichr_results/results_12_18.csv", "new_analysis/enrichr_results/significant_12_18.csv"]
    ):
        pd.concat([up_results[i], down_results[i]]).to_csv(out_fn)

