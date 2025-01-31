rule save_config:
    input:
        config      = "../config/config.yaml",
        metadata    = config["metadata"]
    output:
        config      = os.path.join(RESULT_DIR, "config_used_in_this_experiment.yaml"),
        metadata    = os.path.join(RESULT_DIR, "metadata_used_in_this_experiment.csv")
    shell:
        """
        cp {input.config} {output.config}
        cp {input.metadata} {output.metadata}
        """