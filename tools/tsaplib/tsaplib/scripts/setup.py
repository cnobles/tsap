import os
import sys
import argparse
import subprocess

from ruamel.yaml import YAML
from pathlib import Path

def main( argv = sys.argv ):
    """Create a new project directory with necessary subdirectories."""

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your TsAP "
            "environment and try this command again.")

    usage_str = "\n  tsap %(prog)s <path/to/config.file> <options>"
    
    description_str = (
        "Setup a new TsAP project given a project configuration file."
    )
    
    parser = argparse.ArgumentParser(
        prog = "setup", 
        usage = usage_str,
        description = description_str
    )

    parser.add_argument(
        "config", 
        help = ("name of config file (%(default)s)"),
        metavar = "CONFIG_FILE"
    )

    parser.add_argument(
        "-i", "--tsap_dir", 
        default = os.getenv("TSAP_DIR", os.getcwd()),
        help = "Path to TsAP installation")

    # The remaining args will not be used
    args, remaining = parser.parse_known_args(argv)

    # TsAP directory
    tsap_directory = Path(args.tsap_dir)
    
    if not tsap_directory.exists():
        sys.stderr.write(
            "Error: could not find TsAP directory '{}'.\n".format(
                args.tsap_dir))
        sys.exit(1)
    
    # Load config yaml file
    yaml = YAML(typ = 'safe')
    config = yaml.load(open(args.config, "r"))
    
    analysis_directory = tsap_directory / "analysis" / config['Run_Name']
    
    # Check for existing project directory
    if analysis_directory.exists():
        sys.stderr.write(
            "Error: Project directory currently exists: '{}'.\n".format(
                str(analysis_directory)))
        sys.exit(1)
    else:
        os.makedirs(str(analysis_directory))
    
    # Construct directory tree
    sub_directories = [
        "input_data", "logs", "processData", "output", "reports"
    ]
    
    for sub_dir in sub_directories:
        os.makedirs(str(analysis_directory / sub_dir))
        
    # Check for input files
    read_types = config["Read_Types"] 
    
    if config["Platform"] == "illumina": 
        for type in read_types: 
            check_existing_fastq(Path(config["Seq_Path"]) / config[type])
    elif config["Platform"] == "nanopore":
        for type in read_types: 
            check_existing_fastq(Path(config["Seq_Path"]))
    else:
        sys.stderr.write(
            "Error: Only 'illumina' and 'nanopore' sequencing platforms supported, not: '{}'.\n".format(
                str(config["Platform"])))
        sys.exit(1)

    # Create symbolic link to config
    config_path = Path(args.config).absolute()
    
    if config_path.exists():
        os.symlink(str(config_path), str(analysis_directory / "config.yml"))
    else:
        sys.stderr.write(
            "Error: could not locate aboslute path to config file: '{}'.\n".format(
                str(config_path)))
        sys.exit(1)
    
    if analysis_directory.exists():
        print("  '{}' setup has completed.".format(config["Run_Name"]))
    else:
        sys.stderr.write(
            "Error: could not setup project: '{}'.\n".format(
                str(config["Run_Name"])))
        sys.exit(1)
        

def check_existing_fastq(path, force=False):
    if path.is_file() and not force:
        print("Sample file '{}' found.".format(path))
    else:
        print("Warning: specified sample file '{}' does not exist. "
              "Make sure it exists before running tsap run.".format(path))
