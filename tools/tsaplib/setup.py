from os import getenv, getcwd
from subprocess import run, PIPE
from setuptools import setup, find_packages

def get_tsap_version(with_hash = False):
    tsap_version_path = getenv("TSAP_DIR", getcwd()) + "/.version"
    tsap_version = open(tsap_version_path, "r").readlines()[0].rstrip()
    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
      )
    commit_str = commit_hash.stdout.decode('utf-8').rstrip()
    if with_hash:
        return tsap_version + "+" + commit_str
    else:
        return tsap_version

setup(
    name = "tsap",
    version = get_tsap_version(),
    packages = find_packages(),
    entry_points = { 'console_scripts': [
        'tsap = tsaplib.scripts.command:main'
    ] }
)
