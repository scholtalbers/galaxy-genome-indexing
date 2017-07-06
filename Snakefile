# Copy this file and config.json to a new directory e.g. ``cp Snakefile config.json dm6/``
# Then edit the INPUTS and OUTPUTS section of the config.json

# TODO: add liftOver, etc

import os

# best practice will be to keep the snakefile and config.json in one separate directory to know what has been produced
# with which files. i.e. no need to update this snakefile, just update the config.json
configfile: "config.json"

# INPUTS - Please change in the config.json
LONG_NAME = config["INPUTS"]["LONG_NAME"]
DBKEY = config["INPUTS"]["DBKEY"]
FULL_NAME = "{FULL_NAME} ({DBKEY})".format(**config["INPUTS"])
DOWNLOAD_URL = config["INPUTS"]["DOWNLOAD_URL"]

# OUTPUTS - Might want to change the galaxy_root
GALAXY_ROOT = config["OUTPUTS"]["GALAXY_ROOT"]
INDICES_OUT = config["OUTPUTS"]["INDICES_OUT"]
GENOME_OUT = config["OUTPUTS"]["GENOME_OUT"]

FASTA_OUT = "{0}/{1}/{2}/fasta/".format(GENOME_OUT,LONG_NAME,DBKEY)
FASTA_GZ = "{0}{1}.fa.gz".format(FASTA_OUT, DBKEY)
FASTA = "{0}{1}.fa".format(FASTA_OUT, DBKEY)
FASTA_BASE = "{0}.fa".format(DBKEY)
FASTA_DICT = "{0}{1}.dict".format(FASTA_OUT, DBKEY)
FASTA_INDEX = "{0}{1}.fa.fai".format(FASTA_OUT, DBKEY)
FASTA_LEN = "{0}{1}.len".format(FASTA_OUT, DBKEY)
FASTA_2BIT = "{0}{1}.2bit".format(FASTA_OUT, DBKEY)

BOWTIE_OUT = "{0}/bowtie2_v2_1/{1}/{2}/".format(INDICES_OUT,LONG_NAME, DBKEY)
BOWTIE_FINISHED = "{0}finished.txt".format(BOWTIE_OUT)

GSNAP_BASE_OUT = "{0}/gsnap_v2012_03_23/{1}".format(INDICES_OUT,LONG_NAME)
GSNAP_OUT = "{0}/gsnap_v2012_03_23/{1}/{2}".format(INDICES_OUT,LONG_NAME, DBKEY)
GSNAP_FINISHED = "{0}/finished.txt".format(GSNAP_OUT)

BWA_OUT = "{0}/galaxy_bwa/{1}/{2}/".format(INDICES_OUT, LONG_NAME, DBKEY)
BWA_FINISHED = "{0}finished.txt".format(BWA_OUT)

print(FASTA_DICT, FASTA_LEN, FASTA, BOWTIE_FINISHED, BWA_FINISHED, GSNAP_FINISHED)

rule all:
    input:
        FASTA_DICT,
        FASTA_LEN,
        FASTA,
        BOWTIE_FINISHED,
        GSNAP_FINISHED,
        BWA_FINISHED

rule report:
    input:
        FASTA,
        FASTA_LEN,
        FASTA_DICT,
	FASTA_2BIT,
        BOWTIE_FINISHED,
        GSNAP_FINISHED,
        BWA_FINISHED
    output:
        html = "{0}_report.html".format(DBKEY),
        txt = "add_{0}_to_galaxy.sh".format(DBKEY)
    run:
        from snakemake.utils import report

        txt = open(output.txt, 'w')
        report_text = """
# Generated FASTA, BOWTIE, BWA and GNSAP indices. You still need to add the build key to Galaxy by going through
# these commands::

# BUILD KEY:
echo -e "{DBKEY}\\t{FULL_NAME}\\t{LONG_NAME}\\n" >> {GALAXY_ROOT}/tool-data/shared/ucsc/builds.txt
# FASTA:
echo -e "{DBKEY}\\t{DBKEY}\\t{FULL_NAME}\\t{FASTA}\\n" >> {GALAXY_ROOT}/tool-data/all_fasta.loc
echo -e "{DBKEY}\\t{DBKEY}\\t{FULL_NAME}\\t{FASTA}\\n" >> {GALAXY_ROOT}/tool-data/fasta_indexes.loc
echo -e "index\\t{DBKEY}\\t{FASTA}\\n" >> {GALAXY_ROOT}/tool-data/sam_fa_indices.loc
echo -e "{DBKEY}\\t{DBKEY}\\t{FULL_NAME}\\t{FASTA}\\n" >> {GALAXY_ROOT}/tool-data/picard_index.loc
echo -e "{DBKEY}\\t{FASTA_2BIT}\\n" >> {GALAXY_ROOT}/tool-data/twobit.loc
echo -e "seq\\t{DBKEY}\\t{FASTA_2BIT}\\n" >> {GALAXY_ROOT}/tool-data/alignseq.loc

# GMAP Indices:
echo -e "{DBKEY}\\t{DBKEY}\\t{FULL_NAME}\\t12,13,14,15\\tsplicesites,introns,snps\\tsnps,dbsnp\\t{GSNAP_OUT}\\n" >> {GALAXY_ROOT}/tool-data/gmap_indices.loc
# BWA
echo -e "{DBKEY}\\t{DBKEY}\\t{FULL_NAME}\\t{BWA_OUT}{FASTA_BASE}\\n" >> {GALAXY_ROOT}/tool-data/bwa_mem_index.loc
# BOWTIE
echo -e "{DBKEY}\\t{DBKEY}\\t{FULL_NAME}\\t{BOWTIE_OUT}{DBKEY}\\n" >> {GALAXY_ROOT}/tool-data/bowtie2_indices.loc

# For trackster - you can also replace rsync by cp if needed...
rsync -aq {FASTA_LEN} {GALAXY_ROOT}/tool-data/shared/ucsc/chrom/

""".format(**globals())
        txt.write(report_text)
        txt.close()
        report = report("""
===============================================================
Index Generation for {FULL_NAME}
===============================================================

{report_text}

""", output.html, **input)

rule download_reference:
    output:
        fasta = FASTA_GZ,
        fasta_out = FASTA_OUT
    shell:
        "mkdir -p {output.fasta_out} && wget -O {output.fasta} '{DOWNLOAD_URL}'"

rule expand_fasta:
    input:
        FASTA_GZ
    output:
        fasta_out=FASTA
    shell:
        "gunzip -c {input} > {output.fasta_out}"

rule fasta_index:
    input:
        FASTA
    output:
        FASTA_INDEX
    shell:
        "samtools faidx {input}"

rule sequence_dictionary:
    input:
        FASTA_GZ
    output:
        FASTA_DICT
    shell:
        "java -jar {GALAXY_ROOT}/tool-data/shared/jars/picard/CreateSequenceDictionary.jar \
        REFERENCE={input} OUTPUT={output}"

rule sequence_len:
    input:
        FASTA_INDEX
    output:
        FASTA_LEN
    shell:
        "cut -f1,2 {input} > {output}"

rule sequence_2bit:
    input:
        FASTA
    output:
        FASTA_2BIT
    shell:
        "faToTwoBit-3.4 {input} {output}"

rule bowtie2_index:
    input:
        FASTA
    params: abs_fasta = srcdir(FASTA)
    output:
        finished = BOWTIE_FINISHED,
    shell:
        "mkdir -p {BOWTIE_OUT} && \
        cd {BOWTIE_OUT} && \
        bowtie2-build-2.1.0 {params.abs_fasta} {DBKEY} && \
        cd - && touch {output.finished}"

rule gsnap_index:
    input:
        FASTA
    params: abs_fasta = srcdir(FASTA),
    output:
        finished = GSNAP_FINISHED
    shell:
        "gmap_build-2012_03_23 -d {DBKEY} -D {GSNAP_BASE_OUT} {params.abs_fasta} && \
        ln -s {GSNAP_OUT}/{DBKEY}.ref153positions {GSNAP_OUT}/{DBKEY}.ref12153positions && \
        touch {output.finished}"

# We will actually use the bwa version that galaxy uses for building indexes
rule bwa_index:
    input:
        FASTA
    params:
        abs_fasta = srcdir(FASTA),
        index_type = 'bwtsw'
    output:
        finished = BWA_FINISHED
    shell:
        "mkdir -p {BWA_OUT} && \
        cd {BWA_OUT} && \
        ln -s {params.abs_fasta} && \
        {GALAXY_ROOT}/dependencies/bwa/0.7.10.039ea20639/devteam/package_bwa_0_7_10_039ea20639/5b9aca1e1c07/bwa index -a {params.index_type} {FASTA_BASE} && \
        cd - && \
        touch {output.finished}"
