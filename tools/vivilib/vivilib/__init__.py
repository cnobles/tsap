from os import getenv, getcwd, path, chdir
from subprocess import run, PIPE

def import_sample_info(filePath, sampleColName, delim):
    sampleInfo = {}
    with open(filePath, 'r') as info:
        data = info.readlines()
    listData = [row.rstrip().split(delim) for row in data]
    mcols = listData[0]
    samCol = mcols.index(sampleColName)
    samNames = [row.rstrip().split(delim)[samCol] for row in data[1:]]
    for m in mcols:
        ind = mcols.index(m)
        vals = []
        for row in listData[1:]:
            vals.append(row[ind])
            colData = dict(zip(samNames, vals))
            sampleInfo[m] = colData
    return sampleInfo

def choose_sequence_data(config_input, sampleInfo):
    if "sampleInfo" in config_input:
        colnam = config_input.split(":")[1]
        if not colnam in sampleInfo:
            raise SystemExit(print("Cannot find ", colnam, "in sampleInfo."))
        seq = sampleInfo[colnam]
    else:
        initial_col = list(sampleInfo)[0]
        samples = list(sampleInfo[initial_col])
        seq = dict(zip(samples, [config_input] * len(samples)))
    return seq

def get_vivi_version(with_hash = False):
    vivi_path = getenv("VIVI_DIR", None)

    if vivi_path is None:
        raise SystemExit(
          print("  VIVI_DIR cannot be found as an environmental variable.\n"
                "  Check to make sure your VivI environment is active,   \n"
                "  you may need to restart your environment, update, or    \n"
                "  reinstall VivI with the install.sh script.")
        )
    else:
        vivi_version_path = vivi_path + "/.version"

    if not path.exists(vivi_version_path):
        raise SystemExit(
          print("  VivI version cannot be located. Check environmental\n"
                "  variables, such as VIVI_DIR, otherwise you may want\n"
                "  to restart your environment, update, or reinstall    \n"
                "  VivI using the install.sh script.")
        )

    vivi_version = open(vivi_version_path, "r").readlines()[0].rstrip()

    wd = getcwd()
    chdir(str(vivi_path))

    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
    )

    chdir(wd)

    commit_str = commit_hash.stdout.decode('utf-8').rstrip()

    if with_hash:
        return vivi_version + "+" + commit_str
    else:
        return vivi_version

__version__ = get_vivi_version( with_hash = True )
