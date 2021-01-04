RUNS, = glob_wildcards("fastq/{run}.fastq.gz")

SALMON = "/proj/milovelab/bin/salmon-1.4.0_linux_x86_64/bin/salmon"

rule all:
  input: expand("quants/{run}/quant.sf", run=RUNS)

rule salmon_index:
    input: "gencode.v36.transcripts.fa.gz"
    output: directory("gencode.v36-salmon_1.4.0")
    shell: "{SALMON} index --gencode -p 8 -t {input} -i {output}"

rule salmon_quant:
    input:
        read = "fastq/{sample}.fastq.gz",
        index = "gencode.v36-salmon_1.4.0"
    output:
        "quants/{sample}/quant.sf"
    params:
        dir = "quants/{sample}"
    shell:
        "{SALMON} quant -i {input.index} -l A -p 8 --gcBias "
        "-o {params.dir} -r {input.read}"
