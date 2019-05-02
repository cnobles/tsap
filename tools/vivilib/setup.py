from os import getenv, getcwd
from subprocess import run, PIPE
from setuptools import setup, find_packages

def get_vivi_version(with_hash = False):
    vivi_version_path = getenv("VIVI_DIR", getcwd()) + "/.version"
    vivi_version = open(vivi_version_path, "r").readlines()[0].rstrip()
    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
      )
    commit_str = commit_hash.stdout.decode('utf-8').rstrip()
    if with_hash:
        return vivi_version + "+" + commit_str
    else:
        return vivi_version

setup(
    name = "vivi",
    version = get_vivi_version(),
    packages = find_packages(),
    entry_points = { 'console_scripts': [
        'vivi = vivilib.scripts.command:main'
    ] }
)
