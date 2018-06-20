# Vivi : Validation of in vivo Loci through targeted sequencing of genomic DNA.
#
# Author : Christopher Nobles, Ph.D.

import os
import sys
import re
import yaml
import configparser
from tools.pytools.defs import *

if not config:
    raise SystemExit(
        "No config file specified. Feel free to use the default config as a"
        "template to generate a config file, and specify with --configfile")

# Import sampleInfo
if ".csv" in config["Sample_Info"]:
    delim = ","
elif ".tsv" in config["Sample_Info"]:
    delim = "\t"
else:
    raise SystemExit("Sample Info file needs to contain extention '.csv' or '.tsv'.")

# Sample information
sampleInfo = import_sample_info(
    config["Sample_Info"], config["Sample_Name_Column"], delim)

SAMPLES=sampleInfo[config["Sample_Name_Column"]]
TYPES=config["Read_Types"]
READS=config["Genomic_Reads"]

wildcard_constraints:
    samples=SAMPLES

# Working paths
RUN = config["Run_Name"]
ROOT_DIR = config["Install_Directory"]
RUN_DIR = config["Install_Directory"] + "/analysis/" + RUN

# Check for directory paths.
if not os.path.isdir(ROOT_DIR):
    raise SystemExit("Path to vivi is not found. Check configuration file.")

# Target Rules
rule all:
    input: expand(RUN_DIR + "/processData/{sample}.bai", sample = SAMPLES)

# Architecture Rules
include: "rules/arch.rules"

# Processing Rules
include: "rules/demulti.rules"
include: "rules/trim.rules"
include: "rules/filt.rules"
include: "rules/assembly.rules"
include: "rules/consol.rules"
include: "rules/align.rules"
