import pandas as pd

SAMPLES = (
    pd.read_csv(
        config["samples"], sep="\t", comment="#", dtype={"sample": str, "bcf": str}
    )
    .set_index("sample", drop=False)
    .sort_index()
)
validate(SAMPLES, schema="../schemas/samples.schema.yaml")


def get_sample_bcf(wildcards):
    return SAMPLES.loc[wildcards.sample]["bcf"]


def get_frequent_set_algorithm(wildcards):
    return config["genecart"].get("frequent_sets", dict()).get("algorithm", "fpgrowth")


def get_min_support(wildcards):
    return float(
        config["genecart"].get("frequent_sets", dict()).get("min_support", 0.2)
    )


def get_max_len(wildcards):
    max_len = config["genecart"].get("frequent_sets", dict()).get("max_len", None)
    if max_len:
        return int(max_len)
    else:
        return None


def get_association_rule_metric(wildcards):
    return (
        config["genecart"].get("association_rules", dict()).get("metric", "confidence")
    )


def get_association_rule_threshold(wildcards):
    return float(
        config["genecart"].get("association_rules", dict()).get("threshold", 0.8)
    )


def get_filter_expression(wildcards):
    return config["genecart"].get(
        "filter_expression", 'ANN["IMPACT"] in ("MODERATE", "HIGH")'
    )


def get_symbol_expression(wildcards):
    return config["genecart"].get("symbol_expression", "symbol")


def get_frequent_sets_write_samples(wildcards):
    return config["genecart"].get("frequent_sets", dict()).get("write_samples", True)


def get_association_rules_write_samples(wildcards):
    return (
        config["genecart"].get("association_rules", dict()).get("write_samples", True)
    )
