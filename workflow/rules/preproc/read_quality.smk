# vim: syntax=python expandtab
# Read pre-processing 
# TODO: Remove superfluous str conversions of paths in expand and log statements
#      when Snakemake is pathlib compatible.


run_fp_and_kd = config["qc_reads"]["fastp"] and config["qc_reads"]["kneaddata"]

if run_fp_and_kd:
    err_message = "Running both fastp and kneadata for quality control is not supported"
    raise WorkflowError(err_message)

#########################################
#              fastp
#########################################


fastp_config = config["quality_control"]["fastp"]
if config["qc_reads"]["fastp"]:
    # Add final output files from this module to 'all_outputs' from
    # the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
    trimmed_qc = expand(str(OUTDIR/"fastp/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1, 2])

    if fastp_config["keep_output"]:
        all_outputs.extend(trimmed_qc)

    citations.add(publications["fastp"])

    rule fastp:
        input:
            read1=INPUT_read1,
            read2=INPUT_read2,
        output:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz" if fastp_config["keep_output"] else temp(OUTDIR/"fastp/{sample}_1.fq.gz"),
            read2=OUTDIR/"fastp/{sample}_2.fq.gz" if fastp_config["keep_output"] else temp(OUTDIR/"fastp/{sample}_2.fq.gz"),
            json=LOGDIR/"fastp/{sample}.fastp.json",
            html=LOGDIR/"fastp/{sample}.fastp.html",
        log:
            stdout=str(LOGDIR/"fastp/{sample}.stdout.log"),
            stderr=str(LOGDIR/"fastp/{sample}.stderr.log"),
        shadow:
            "shallow"
        conda:
            config["conda"] if config["conda"] else "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads: 4
        params:
            extra=fastp_config["extra"],
        shell:
            """
            fastp \
                --in1 {input.read1} \
                --in2 {input.read2} \
                --out1 {output.read1} \
                --out2 {output.read2} \
                --json {output.json} \
                --html {output.html} \
                --thread {threads} \
                {params.extra} \
                > {log.stdout} \
                2> {log.stderr}
            """

#########################################
#              kneaddata
#########################################

kd_config = config["quality_control"]["kneaddata"]
if config["qc_reads"]["kneaddata"]:
    # Add final output files from this module to 'all_outputs' from
    # the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
    processed_qc = expand(str(OUTDIR/"kneaddata/{sample}_{readpair}.fastq"),
            sample=SAMPLES,
            readpair=[1, 2])

    if fastp_config["keep_output"]:
        all_outputs.extend(processed_qc)

    citations.add(publications["kneaddata"])
    citations.add(publications["Trimmomatic"])
    citations.add(publications["TRF"])
    citations.add(publications["FastQC"])
    citations.add(publications["Bowtie2"])

    rule kneaddata:
        input:
            read1=INPUT_read1,
            read2=INPUT_read2,
        output:
            paired1=OUTDIR/"kneaddata/{sample}_1.fastq" if kd_config["keep_fastq"] else temp(OUTDIR/"kneaddata/{sample}_1.fastq"),
            paired2=OUTDIR/"kneaddata/{sample}_2.fastq" if kd_config["keep_fastq"] else temp(OUTDIR/"kneaddata/{sample}_2.fastq"),
            unmatched1=OUTDIR/"kneaddata/{sample}_unmatched_1.fastq" if kd_config["keep_unmatched"] else temp(OUTDIR/"kneaddata/{sample}_unmatched_1.fastq"),
            unmatched2=OUTDIR/"kneaddata/{sample}_unmatched_2.fastq" if kd_config["keep_unmatched"] else temp(OUTDIR/"kneaddata/{sample}_unmatched_2.fastq"),
            runlog=LOGDIR/"kneaddata/{sample}.log",
            fastqcout=directory(OUTDIR/"fastqc/{sample}"),
        log:
            stdout=str(LOGDIR/"kneaddata/{sample}.kneaddata.log"),
        shadow:
            "shallow"
        conda:
            config["conda"] if config["conda"] else "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads: 8
        params:
            db=kd_config["db_path"],
            outdir=lambda w: f"{OUTDIR}/kneaddata/{w.sample}_out",
            output1=lambda w: f"{OUTDIR}/kneaddata/{w.sample}_out/{w.sample}_paired_1.fastq",
            output2=lambda w: f"{OUTDIR}/kneaddata/{w.sample}_out/{w.sample}_paired_2.fastq",
            output3=lambda w: f"{OUTDIR}/kneaddata/{w.sample}_out/{w.sample}_unmatched_1.fastq",
            output4=lambda w: f"{OUTDIR}/kneaddata/{w.sample}_out/{w.sample}_unmatched_2.fastq",
            outputlog=lambda w: f"{OUTDIR}/kneaddata/{w.sample}_out/{w.sample}.log",
            outprefix=lambda w: f"{w.sample}",
            outputfastqc=lambda w: f"{OUTDIR}/kneaddata/{w.sample}_out/fastqc/*",
            trimmomatic=kd_config["trim_jar"],
            extra=kd_config["extra"],
        shell:
            """
            kneaddata \
                --input1 {input.read1} \
                --input2 {input.read2} \
                -db {params.db} \
                --output {params.outdir} \
                --run-trim-repetitive \
                --run-fastqc-start \
                --threads {threads} \
                --output-prefix {params.outprefix} \
                --trimmomatic {params.trimmomatic}\
                {params.extra} \
            > {log.stdout}

            mv {params.output1} {output.paired1}
            mv {params.output2} {output.paired2}
            mv {params.output3} {output.unmatched1}
            mv {params.output4} {output.unmatched2}
            mv {params.outputlog} {output.runlog}
            mkdir {output.fastqcout}
            mv {params.outputfastqc} {output.fastqcout}

            rm -r {params.outdir}

            """


###no filtering process

if not (config["qc_reads"]["kneaddata"] or config["qc_reads"]["fastp"]):
    trimmed_qc = expand(str(OUTDIR/"fastp/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1, 2])

    if fastp_config["keep_output"]:
        all_outputs.extend(trimmed_qc)

    localrules:
        skip_fastp,

    rule skip_fastp:
        input:
            read1=INPUT_read1,
            read2=INPUT_read2,
        output:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        log:
            stderr=str(LOGDIR/"fastp/{sample}.stderr.log"),
        shell:
            """
            ln -sv $(readlink -f {input.read1}) {output.read1} >> {log.stderr}
            ln -sv $(readlink -f {input.read2}) {output.read2} >> {log.stderr}
            """
