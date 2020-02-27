rule run_coloc:
    input:
        "output/simulation/{simulation}/{settings}/gene_{gene}/data"
    output:
        "output/simulation/{simulation}/{settings}/gene_{gene}/coloc"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/run_coloc.R"

