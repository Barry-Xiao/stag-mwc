# vim: syntax=python expandtab
# Metagenomic strain-level population genomics with StrainPhlAn.
 
from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:

mpa_config = config["metaphlan"]
spa_config = config["strainphlan"]
if config["strain_level_profiling"]["strainphlan"]:
    if not mpa_config["bt2_db_dir"] or not Path(mpa_config["bt2_db_dir"]).exists():
        err_message = "No MetaPhlAn database dir found at: '{}'!\n".format(mpa_config["bt2_db_dir"])
        err_message += "bt2_db_dir and bt2_index are required to run StrainPhlAn as it uses the output from MetaPhlAn as input.\n"
        err_message += "Specify relevant paths in the metaphlan section of config.yaml.\n"
        err_message += "If you do not want to run MetaPhlAn or StrainPhlAn, set 'metaphlan: False' and 'strainphlan: false' in config.yaml"
        raise WorkflowError(err_message)
    if not spa_config["clade_of_interest"]:
        available_clades=f"{OUTDIR}/strainphlan/print_clades_only.tsv",
        all_outputs.append(available_clades)
        user_messages.warn("Clade of interest not specified in strainphlan section of config.yaml.")
        user_messages.warn("Based on your samples strainphlan will create a list of available clades in output/strainphlan/print_clades_only.tsv.")
        user_messages.warn("If you still want to run strainphlan, please update config.yaml e.g. \"clade_of_interest: s__Bifidobacterium_longum\".")
    if spa_config["clade_of_interest"]:
        n_sample = len(SAMPLES)
        min_prop = float(spa_config["clade_min_prop"])
        min_num = round(n_sample*min_prop)
        available_clades = []
        with open(f"{OUTDIR}/strainphlan/print_clades_only.tsv") as handler:
            for line in handler.readlines():
                if re.match("^Clade.*",line):
                    continue
                temp = line.rstrip().split('\t')
                if int(temp[1]) >= min_num:
                    available_clades.append(temp[0])

        if len(available_clades) == 0:
            err_message = f"No clade met the criteria of being found in {min_prop*100}% of total samples"
            raise WorkflowError(err_message)

        ref_markers = expand(f"{OUTDIR}/strainphlan/ref_markers/{{clade}}.fna", clade=available_clades),
        spa_alignment = expand(f"{OUTDIR}/strainphlan/{{clade}}/{{clade}}.StrainPhlAn4_concatenated.aln", clade=available_clades),
        spa_tree=expand(f"{OUTDIR}/strainphlan/{{clade}}/RAxML_bestTree.{{clade}}.StrainPhlAn4.tre", clade=available_clades),
        all_outputs.append(spa_alignment)
        all_outputs.append(spa_tree)
        user_messages.info("If strainphlan failed, ensure your clade_of_interest is present in output/strainphlan/print_clades_only.tsv.")
        citations.add(publications["MetaPhlAn"])
        citations.add(publications["StrainPhlAn"])


rule consensus_markers:
    """Generate consensus markers"""
    input:
        sam=f"{OUTDIR}/metaphlan/{{sample}}.sam.bz2", 
    output:
        consensus_markers=f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}/{{sample}}.pkl",
    log:
        stdout=f"{LOGDIR}/strainphlan/sample2markers.{{sample}}.strainphlan.stdout.log",
        stderr=f"{LOGDIR}/strainphlan/sample2markers.{{sample}}.strainphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        config["conda"] if config["conda"] else "../../envs/metaphlan.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:biobakery"+singularity_branch_tag
    threads: 8
    params:
        output_dir=f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}",
        database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_index']}.pkl",
    shell:
        """
        sample2markers.py \
             -d {params.database} \
             -i {input.sam} \
             -o {params.output_dir} \
             -n 8 \
             > {log.stdout} \
             2> {log.stderr}
        """

rule print_clades:
    """Print available clades"""
    input:
        consensus_markers=expand(f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}/{{sample}}.pkl", sample=SAMPLES),
    output:
       available_clades=f"{OUTDIR}/strainphlan/print_clades_only.tsv",
    log:
        stdout=f"{LOGDIR}/strainphlan/available_clades.stdout",
        stderr=f"{LOGDIR}/strainphlan/available_clades.stderr",
    shadow:
        "shallow"
    conda:
        config["conda"] if config["conda"] else "../../envs/metaphlan.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:biobakery"+singularity_branch_tag
    threads: 8
    params:
        out_dir=f"{OUTDIR}/strainphlan",
        database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_index']}.pkl",
    shell:
        """
        strainphlan \
             -s {input.consensus_markers} \
             --print_clades_only \
             -d {params.database} \
             -o {params.out_dir} \
             -n {threads} \
             > {log.stdout} \
             2> {log.stderr}

        """

rule extract_markers:
    """Extract marker sequences for clade of interest"""
    input:
        available_clades=f"{OUTDIR}/strainphlan/print_clades_only.tsv",
        consensus_markers=expand(f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}/{{sample}}.pkl", sample=SAMPLES),
    output:
        reference_markers=f"{OUTDIR}/strainphlan/ref_markers/{{clade}}.fna",
    log:
        stdout=f"{LOGDIR}/strainphlan/extract_markers.{{clade}}.strainphlan.stdout.log",
        stderr=f"{LOGDIR}/strainphlan/extract_markers.{{clade}}.strainphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        config["conda"] if config["conda"] else "../../envs/metaphlan.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:biobakery"+singularity_branch_tag
    threads: 8
    params:
        clade=f"{{clade}}",
        out_dir=f"{OUTDIR}/strainphlan/ref_markers",
        database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_index']}.pkl",
    shell:
        """
        extract_markers.py \
             -c {params.clade} \
             -o {params.out_dir} \
             -d {params.database} \
             > {log.stdout} \
             2> {log.stderr}
        """

rule strainphlan:
    """Generate tree and alignment"""
    input:
        consensus_markers=expand(f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}/{{sample}}.pkl", sample=SAMPLES),
        reference_markers=f"{OUTDIR}/strainphlan/ref_markers/{{clade}}.fna",
    output:
        alignment=f"{OUTDIR}/strainphlan/{{clade}}/{{clade}}.StrainPhlAn4_concatenated.aln",
        tree=f"{OUTDIR}/strainphlan/{{clade}}/RAxML_bestTree.{{clade}}.StrainPhlAn4.tre",
    log:
        stdout=f"{LOGDIR}/strainphlan/alignment.{{clade}}.strainphlan.stdout.log",
        stderr=f"{LOGDIR}/strainphlan/alignment.{{clade}}.strainphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        config["conda"] if config["conda"] else "../../envs/metaphlan.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:biobakery"+singularity_branch_tag
    threads: 8
    params:
        clade=f"{{clade}}",
        out_dir=f"{OUTDIR}/strainphlan/{{clade}}",
        database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_index']}.pkl",
        mode=spa_config["mode"],
        extra=spa_config["extra"],  # This is extremely useful if you want to include a reference genome
    shell:
        """
        echo "please compare your clade_of_interest to list of available clades in print_clades_only.tsv" > {log.stderr}

        strainphlan \
             -s {input.consensus_markers} \
             -m {input.reference_markers} \
             {params.extra} \
             -d {params.database} \
             -o {params.out_dir} \
             -n {threads} \
             -c {params.clade} \
             --phylophlan_mode accurate \
             --mutation_rates \
             > {log.stdout} \
             2>> {log.stderr}
        """

