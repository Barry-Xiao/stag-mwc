# vim: syntax=python expandtab
# MultiQC

spa_config = config["strainphlan"]

if (config["multiqc_report"] and
    not spa_config["clade_of_interest"]):
    citations.add(publications["MultiQC"])

    mqc_config = config["multiqc"]
    rule multiqc:
        input:
            all_outputs
        output:
            report=report(f"{OUTDIR}/multiqc/multiqc_report.html",
                category="Sequencing data quality",
                caption="../../report/multiqc.rst"),
        log:
           f"{LOGDIR}/multiqc/multiqc.log"
        shadow:
            "shallow"
        conda:
            config["conda"] if config["conda"] else "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads: 1
        params:
            extra=mqc_config["extra"],
        shell:
            """
            multiqc {OUTDIR} \
                --filename {output.report} \
                --force \
                2> {log}
            """

    # Appended after the rule definition to avoid circular dependency
    all_outputs.append(f"{OUTDIR}/multiqc/multiqc_report.html")
