import pandas as pd
from typing import List, Tuple
from functools import reduce
from mlxtend.frequent_patterns import (
    fpgrowth,
    fpmax,
    apriori,
    association_rules,
)
from mlxtend.preprocessing import TransactionEncoder
from contextlib import redirect_stdout, redirect_stderr


def read_samples(path: str) -> List[str]:
    with open(path, "rt") as s:
        samples = [sample.strip() for sample in s.readlines()]
    return samples


def read_transactions(path: str) -> List[List[str]]:
    with open(path, "rt") as t:
        transactions = [l.strip().split(",") for l in t.readlines()]
    return transactions


def genecart(
    samples: List[str], transactions: List[List[str]], **params
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    assert len(samples) > 1
    assert len(samples) == len(transactions)
    observed_genes = reduce(set.union, map(set, transactions))
    print(f"{len(samples)} samples, {len(observed_genes)} genes")

    encoder = TransactionEncoder()
    matrix = encoder.fit_transform(transactions, sparse=True)
    print(f"encoded transactions")

    df = pd.DataFrame.sparse.from_spmatrix(
        matrix, index=samples, columns=encoder.columns_
    )
    print("built pandas sparse dataframe")
    print(df)

    algo = params["algorithm"].lower()
    algorithm = (
        fpgrowth if algo == "fpgrowth" else (fpmax if algo == "fpmax" else apriori)
    )

    print(f"determining frequent sets with {str(algorithm)}")
    frequent_sets = algorithm(
        df,
        min_support=params["min_support"],
        use_colnames=True,
        max_len=params["max_len"],
    )
    print("finished finding frequent sets")

    print("generating association rules")
    rules = association_rules(
        frequent_sets,
        metric=params["association_rule_metric"],
        min_threshold=params["association_rule_threshold"],
    )
    print("finished generating association rules")

    return frequent_sets, rules


def main(snakemake):
    samples = read_samples(snakemake.input.samples)
    transactions = read_transactions(snakemake.input.transactions)
    frequent_sets, rules = genecart(samples, transactions, **snakemake.params)
    frequent_sets.to_csv(snakemake.output.frequent_sets, sep="\t")
    rules["antecedents"] = rules["antecedents"].apply(lambda s: ','.join(map(str, s)))
    rules["consequents"] = rules["consequents"].apply(lambda s: ','.join(map(str, s)))
    rules.to_csv(snakemake.output.association_rules, sep="\t")



if __name__ == "__main__":
    import sys
    with open(snakemake.log[0], "wt") as f:
        with redirect_stdout(f), redirect_stderr(f):
            main(snakemake)
