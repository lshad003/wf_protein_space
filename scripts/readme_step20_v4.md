### Step 20. Taxonomic context of the unknown fraction

Each gene inherited the taxonomy of its contig, from existing UniRef50-based
contig assignments, to ask whether unannotated genes sit on less classifiable DNA.

**Result.** Contig taxonomy is available for 48 of 88 metagenomes, covering
888,183 of 2,067,011 genes (42.97%). Contigs are unclassified for 74.6% of K and
KWP genes, 74.4% of GU and 91.6% of EU genes. EU genes are 4.5% eukaryotic and
3.9% bacterial, against about 23% bacterial and 2% eukaryotic in the other
classes. Genes with no database match sit on DNA that is itself largely
unclassified, consistent with lineages lacking sequenced relatives.

| File | Purpose |
|---|---|
| `scripts/step20c_gene_tax_n88.py` | Contig taxonomy by gene and class |

Output: `results/gene_taxonomy_n88.tsv`

