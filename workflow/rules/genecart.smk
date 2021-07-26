rule genecart:
    input:
        transactions="results/transactions/genes_per_sample.list",
        samples="results/transactions/samples.list",
    output:
        frequent_sets="results/genecart/frequent_sets.tsv",
        association_rules="results/genecart/association_rules.tsv",
    log:
        "logs/genecart.log",
    params:
        min_support=get_min_support,
        max_len=get_max_len,
        algorithm=get_frequent_set_algorithm,
        association_rule_metric=get_association_rule_metric,
        association_rule_threshold=get_association_rule_threshold,
    conda:
        "../envs/genecart.yaml"
    script:
        "../scripts/genecart.py"


rule aggregate_gene_lists:
    input:
        lists=expand("results/genes/{sample}.list", sample=SAMPLES),
    output:
        transactions="results/transactions/genes_per_sample.list",
        samples="results/transactions/samples.list",
    params:
        samples=SAMPLES,
    conda:
        "../envs/aggregate_gene_lists.yaml"
    log:
        "logs/aggregate_gene_lists.log",
    script:
        "../scripts/aggregate_gene_lists.py"


rule extract_genes:
    input:
        bcf="results/variants/{sample}.bcf",
    output:
        list="results/genes/{sample}.list",
    params:
        filter=get_filter_expression,
        symbol=get_symbol_expression,
    log:
        "logs/extract_genes/{sample}.log",
    conda:
        "../envs/vembrane.yaml"
    shell:
        """
        vembrane filter '{params.filter}' {input.bcf} | vembrane table '{params.symbol}' | tail -n +2 | sort | uniq > {output.list} 2> {log}
        """
