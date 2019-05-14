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

def get_tsap_version(with_hash = False):
    tsap_path = getenv("TSAP_DIR", None)

    if tsap_path is None:
        raise SystemExit(
          print("  TSAP_DIR cannot be found as an environmental variable.\n"
                "  Check to make sure your TsAP environment is active,   \n"
                "  you may need to restart your environment, update, or    \n"
                "  reinstall TsAP with the install.sh script.")
        )
    else:
        tsap_version_path = tsap_path + "/.version"

    if not path.exists(tsap_version_path):
        raise SystemExit(
          print("  TsAP version cannot be located. Check environmental\n"
                "  variables, such as TSAP_DIR, otherwise you may want\n"
                "  to restart your environment, update, or reinstall    \n"
                "  TsAP using the install.sh script.")
        )

    tsap_version = open(tsap_version_path, "r").readlines()[0].rstrip()

    wd = getcwd()
    chdir(str(tsap_path))

    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
    )

    chdir(wd)

    commit_str = commit_hash.stdout.decode('utf-8').rstrip()

    if with_hash:
        return tsap_version + "+" + commit_str
    else:
        return tsap_version

__version__ = get_tsap_version( with_hash = True )
