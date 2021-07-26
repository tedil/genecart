(SAMPLES,) = glob_wildcards("results/variants/{sample}.bcf")
assert len(SAMPLES) > 1


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
