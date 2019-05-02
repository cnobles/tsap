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
            "Could not determine Conda prefix. Activate your VivI "
            "environment and try this command again.")

    usage_str = "\n  vivi %(prog)s <path/to/config.file> <options>"
    
    description_str = (
        "Setup a new VivI project given a project configuration file."
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
        "-i", "--vivi_dir", 
        default = os.getenv("VIVI_DIR", os.getcwd()),
        help = "Path to VivI installation")

    # The remaining args will not be used
    args, remaining = parser.parse_known_args(argv)

    # VivI directory
    vivi_directory = Path(args.vivi_dir)
    
    if not vivi_directory.exists():
        sys.stderr.write(
            "Error: could not find VivI directory '{}'.\n".format(
                args.vivi_dir))
        sys.exit(1)
    
    # Load config yaml file
    yaml = YAML(typ = 'safe')
    config = yaml.load(open(args.config, "r"))
    
    analysis_directory = vivi_directory / "analysis" / config['Run_Name']
    
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
    
    output_sub_directories = ["unique_aligns", "chimera_aligns"]
    
    for sub_dir in sub_directories:
        os.makedirs(str(analysis_directory / sub_dir))
        
    for out_sub in output_sub_directories:
        os.makedirs(str(analysis_directory / "output" / out_sub))

    # Check for input files
    read_types = config["Read_Types"] 
     
    for type in read_types: 
        check_existing_fastq(Path(config["Seq_Path"]) / config[type])

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
              "Make sure it exists before running vivi run.".format(path))
