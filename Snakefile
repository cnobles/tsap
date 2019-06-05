# TsAP : Targeted sequencing Analysis Pipeline.
#
# Author : Christopher Nobles, Ph.D.

import os
import sys
import re
import yaml
import configparser
from pathlib import Path
from tsaplib import import_sample_info, choose_sequence_data


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
    config["Sample_Info"], config["Sample_Name_Column"], delim
)

SAMPLES=sampleInfo[config["Sample_Name_Column"]]
TYPES=config["Read_Types"]
READS=config["Genomic_Reads"]
PANEL=str(".").join(str(Path(config["Panel_Path"]).name).split(".")[:-1])

wildcard_constraints:
    sample="[^/^\\.]+"

# Working paths
RUN = config["Run_Name"]
ROOT_DIR = ""
try:
    ROOT_DIR = os.environ["TSAP_DIR"]
except KeyError:
    raise SystemExit(
        "TSAP_DIR environment variable not defined. Are you sure you "
        "activated the tsap conda environment?")
RUN_DIR = ROOT_DIR + "/analysis/" + RUN

# Check for sequencing platform.
if not "Platform" in config:
    config["Platform"] = "illumina"

# Check for directory paths.
if not os.path.isdir(ROOT_DIR):
    raise SystemExit("Path to tsap is not found. Check configuration file.")

# Set parameters
BC5p=choose_sequence_data(config["BC5p"], sampleInfo)
BC3p=choose_sequence_data(config["BC3p"], sampleInfo)

BC5list=[]
for bc in BC5p:
    BC5list.append(len(BC5p[bc]))
BC5len=int(set(BC5list).pop())

BC3list=[]
for bc in BC3p:
    BC3list.append(len(BC3p[bc]))
BC3len=int(set(BC3list).pop())


# Summary Report File
report_output = RUN_DIR + "/reports/report." + RUN
if (config["format"] == "pdf"):
    report_output = report_output + ".pdf"
elif (config["format"] == "html"):
    report_output = report_output + ".html"

# Target Rules
rule all:
    input: report_output

# Processing Rules
if (config["Platform"] == "illumina"):
    include: "rules/demulti.rules"
    include: "rules/trim.rules"
    include: "rules/filt.rules"
    include: "rules/assembly.rules"
    include: "rules/consol.rules"
    include: "rules/align.rules"
    include: "rules/process.rules"
elif (config["Platform"] == "nanopore"):
    include: "rules/trim.np.rules"
    include: "rules/filt.np.rules"
    include: "rules/consol.np.rules"
    include: "rules/align.np.rules"
    include: "rules/process.np.rules"
else:
    raise SystemExit("Sequencing [Platform](s) supported: 'illumina' and 'nanopore'.")
