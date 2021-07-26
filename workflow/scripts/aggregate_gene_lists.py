samples = snakemake.params.samples
with open(snakemake.output.transactions, "wt") as out:
    for l in snakemake.input.lists:
        genes = [l.strip() for l in open(l, "rt").readlines()]
        print(",".join(gene for gene in genes if gene), file=out, end="\n")
with open(snakemake.output.samples, "wt") as out:
    print("\n".join(samples), file=out)