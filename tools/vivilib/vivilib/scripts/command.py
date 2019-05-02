import sys
import argparse
import subprocess

from vivilib import __version__
from vivilib.scripts.run import main as Run
from vivilib.scripts.setup import main as Setup
#from vivilib.scripts.config import main as Config
from vivilib.scripts.clean import main as Clean

def main():

    usage_str = "\n  %(prog)s [-h/--help,-v/--version] <subcommand> <path/to/config.file> <options> -- <snakemake.options>"
    description_str = (
        "list of subcommands:\n\n"
        "  primary:\n"
        "    setup        \tCreate a new project directory using config file.\n"
        "    run          \tExecute the VivI pipeline.\n\n"
        "  accessory:\n"
        "    clean        \tCleanup project directory to reduce size. Keeps terminal files.\n\n"
    ).format(version=__version__)

    parser = argparse.ArgumentParser(
        prog = "vivi",
        usage = usage_str,
        description = description_str,
        epilog = "For more help, see the docs at http://vivi.readthedocs.io.",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        add_help = False
    )

    parser.add_argument(
        "command", default = "None", help = argparse.SUPPRESS, nargs = "?"
    )
    
    parser.add_argument(
        "-v", "--version", action = "version",
        version = "%(prog)s {}".format(__version__)
    )

    args, remaining = parser.parse_known_args()
    
    sub_cmds = ["setup", "run", "clean"]
    
    if not args.command in sub_cmds:
        parser.print_help()
        if not args.command in ['None']:
            sys.stderr.write("  Unrecognized subcommand, '{}'.\n".format(
                args.command
            ))
        sys.exit(1)

    if args.command == "setup":
        Setup(remaining)
    elif args.command == "run":
        Run(remaining)
    elif args.command == "clean":
        Clean(remaining)
    else:
        parser.print_help()
