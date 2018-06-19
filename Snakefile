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

# Working paths
RUN = config["Run_Name"]
ROOT_DIR = config["Install_Directory"]
RUN_DIR = config["Install_Directory"] + "/analysis/" + RUN

# Check for directory paths.
if not os.path.isdir(ROOT_DIR):
    raise SystemExit("Path to iGUIDE is not found. Check configuration file.")

# Target Rules
rule all:
    input: RUN_DIR + "/reports/vivi.report." + RUN + ".html"

# Architecture Rules
include: "rules/arch.rules"

# Processing Rules
include: "rules/demulti.rules"
include: "rules/trim.rules"
include: "rules/filt.rules"
include: "rules/consol.rules"
if (config["Aligner"] == "BLAT" or config["Aligner"] == "blat"):
    include: "rules/align.blat.rules"
elif (config["Aligner"] == "BWA" or config["Aligner"] == "bwa"):
    raise SystemExit("BWA aligner not supported yet.")
else:
    "Aligner: " + config["Aligner"] + " not supported."
    "Please choose a supported option: BLAT or BWA."
include: "rules/process.rules"

