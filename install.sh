#!/usr/bin/env bash

__conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

read -r -d '' __usage <<-'EOF'
  -e --environment [arg]  Environment to install to. Default: "tsap"
  -s --tsap_dir [arg]     Location of tsap source code. Default: this directory
  -c --conda [arg]        Location of Conda installation. Default: ${PREFIX}
  -u --update [arg]       Update tsap [lib]rary, conda [env], or [all].
  -r --requirements       Install from requirements rather than build (slow).
  -t --test               After installation, run test to check functionality.
  -v --verbose            Show subcommand output
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

read -r -d '' __helptext <<-'EOF'
 This script installs or upgrades TsAP, including Conda (if not installed).
 To upgrade, pass the '--upgrade all' option, then be sure to update your config
 files using 'tsap config update'.
EOF


# Load BASH3Boilerplate for command-line parsing and logging
source etc/b3bp.sh

function __err_report() {
    local error_code
    error_code=${?}
    error "Error in ${__file} in function ${1} on line ${2}"
    exit ${error_code}
}
trap '__err_report "${FUNCNAME:-.}" ${LINENO}' ERR  

# help mode
if [[ "${arg_h:?}" = "1" ]]; then
    # Help exists with code 1
    help "Help using ${0}:"
fi

# verbose mode
if [[ "${arg_v:?}" = "1" ]]; then
    LOG_LEVEL="7"
fi

# debug mode
if [[ "${arg_d:?}" = "1" ]]; then
    set -o xtrace
    LOG_LEVEL="7"
fi

function debug_capture () {
    debug "$(echo -e "$(${@})")"
}

function installation_error () {
    error "${1} failed!"
    if [[ "${arg_v:?}" != 1 && "${arg_d:?}" != 1 ]]; then
        error "Try re-running with -v or -d, or file an issue on Github."
    fi
    exit 1
}

# Set variables
__conda_path="${arg_c:-${HOME}/miniconda3}"
__tsap_dir="${arg_s:-$(readlink -f ${__dir})}"
__tsap_env="${arg_e:-tsap}"
__run_tsap_tests=false
__reqs_install=false
__update_lib=false
__update_env=false
__req_r_version="3.4.1"
__old_path=$PATH
__output=${2-/dev/stdout}

PATH=$PATH:${__conda_path}/bin

if [[ "${arg_t:?}" = "1" ]]; then
    __run_tsap_tests=true
fi

if [[ "${arg_r:?}" = "1" ]]; then
    __reqs_install=true
fi

if [[ "${arg_u}" = "all" || "${arg_u}" = "env" ]]; then
    __update_lib=true
    __update_env=true
elif [[ "${arg_u}" = "lib" ]]; then
    __update_lib=true
fi

function __test_conda() {
    command -v conda &> /dev/null && echo true || echo false
}

function __detect_conda_install() {
    local discovered=$(__test_conda)

    if [[ $discovered = true ]]; then
      	local conda_path="$(which conda)"
        echo ${conda_path%'/condabin/conda'}
    fi
}

function __test_env() {
    if [[ $(__test_conda) = true ]]; then
        $(conda env list \
        | cut -f1 -d' ' \
        | grep -Fxq $__tsap_env > /dev/null) && \
        echo true || echo false
    else
      	echo false
    fi
}

function __test_r_version () {
    activate_tsap

    local sem_version=$(R --version | grep 'R version' | cut -d ' ' -f 3)

    if (( ${sem_version//./} >= ${__req_r_version//./} )); then
        echo true
    else
        echo false
    fi

    deactivate_tsap
}

function __test_r_packages () {
    activate_tsap

    $(Rscript ${__tsap_dir}/tools/rscripts/check_for_required_packages.R \
        > /dev/null) && echo true || echo false

    deactivate_tsap
}

function __test_tsaplib() {
    if [[ $(__test_env) = true ]]; then
      	activate_tsap
      	command -v tsap &> /dev/null && echo true || echo false
      	deactivate_tsap
    else
      	echo false
    fi
}

function __test_tsap() {
    if [[ $(__test_env) = true ]]; then
      	$(bash ${__tsap_dir}/etc/tests/test.sh ${__tsap_env} &> /dev/null) && \
      	    echo true || echo false
    else
      	echo "fail"
    fi
}

function __test_bashrc () {
  grep conda.sh ~/.bashrc > /dev/null && echo true || echo false
}

function activate_tsap () {
    set +o nounset
    conda activate $__tsap_env
    set -o nounset
}

function deactivate_tsap () {
    set +o nounset
    conda deactivate
    set -o nounset
}

function install_conda () {
    local tmpdir=$(mktemp -d)

    debug "Downloading miniconda..."
    debug_capture wget -q ${__conda_url} -O ${tmpdir}/miniconda.sh 2>&1
    debug "Installing miniconda..."
    debug_capture bash ${tmpdir}/miniconda.sh -b -p ${__conda_path} 2>&1

    # Source conda into path
    source ${__conda_path}/etc/profile.d/conda.sh

    if [[ $(__test_conda) != true ]]; then
      	installation_error "Environment creation"
    fi

    rm ${tmpdir}/miniconda.sh
}

function install_environment () {
    if [[ $__reqs_install == "true" ]]; then
        local install_options="--quiet --file etc/requirements.yml"
        debug_capture conda env update --name=$__tsap_env ${install_options} 2>&1
    else
        local install_options="--quiet --yes --file etc/build.b0.3.0.txt"
        debug_capture conda create --name=$__tsap_env ${install_options} 2>&1
    fi

    if [[ $(__test_env) != true ]]; then
      	installation_error "Environment creation"
    fi

    if [[ $(__test_r_version) != true ]]; then
        installation_error "Insufficient R-program version"
    fi

    if [[ $(__test_r_packages) != true ]]; then
        installation_error "R-package installation"
    fi
}

function install_env_vars () {
    activate_tsap

    echo -ne "#/bin/sh\nexport TSAP_DIR=${__tsap_dir}" > \
	      ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

    mkdir -p ${CONDA_PREFIX}/etc/conda/deactivate.d/

    echo -ne "#/bin/sh\nunset TSAP_DIR" > \
	      ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh

	  deactivate_tsap
}

function install_tsaplib () {
    activate_tsap

    debug_capture pip install --upgrade ${__tsap_dir}/tools/tsaplib/ 2>&1

    if [[ $(__test_tsaplib) != true ]]; then
      	installation_error "Library installation"
    fi

    deactivate_tsap
}

info "Starting TsAP installation..."
info "    Conda path:  ${__conda_path}"
info "    TsAP src:  ${__tsap_dir}"
info "    TsAP env:  '${__tsap_env}'"

debug "Components detected:"
__conda_installed=$(__test_conda)
debug "    Conda:         ${__conda_installed}"
__env_exists=$(__test_env)
debug "    Environment:   ${__env_exists}"
__tsaplib_installed=$(__test_tsaplib)
debug "    Library:       ${__tsaplib_installed}"

__env_changed=false


# Install Conda if necessary
if [[ $__conda_installed = true ]]; then
    if [[ $(__detect_conda_install) != $__conda_path ]]; then
        warning "Found pre-existing Conda installation in $(__detect_conda_install)".
        warning "Ignoring specified Conda path in favor of existing Conda install."
        __conda_path=$(__detect_conda_install)
    fi
    info "Conda already installed."
else
    info "Installing Conda..."
    install_conda
    __env_changed=true
fi

# Source conda into shell
source ${__conda_path}/etc/profile.d/conda.sh

# Create Conda environment for TsAP
if [[ $__env_exists = true && $__update_env = false ]]; then
    info "Specified environment already exists (use '--update env' to update)"
else
    if [[ $__reqs_install = "true" ]]; then
        __build_source="etc/requirements.yml"
    else
        __build_source="etc/build.b0.3.0.txt"
    fi

    info "Creating TsAP environment..."
    info "    Building from: $__build_source"
    install_environment
    __env_changed=true
    info "$__tsap_env environment created."
fi


# Always update the env_vars.sh in the TsAP environment
debug "Updating \$TSAP_DIR variable to point to ${__tsap_dir}"
info "Setting environmental variables..."
install_env_vars


# Install tsaplib into environment if changed or requested
if [[ $__env_changed = true ]]; then
    info "Environment installed/updated; (re)installing TsAP library..."
    install_tsaplib
elif [[ $__tsaplib_installed = false ]]; then
    info "Installing TsAP library..."
    install_tsaplib
elif [[ $__update_lib = true ]]; then
    info "Updating TsAP library..."
    install_tsaplib
else
    info "TsAP library already installed (use '--update lib' to update)"
fi


# Run tests if requested
if [[ $__run_tsap_tests = true ]]; then
    info "Running TsAP tests...(this may take a few mins)"
    __tsap_tested=$(__test_tsap)
    
    if [[ $__tsap_tested = true ]]; then
        info "    TsAP Tests:  passed"
    else
        warning "    TsAP Tests:  FAILED"
        warning "    Try running the test outside of the install to confirm."
        warning "    Just run 'bash etc/tests/test.sh $__tsap_env'."
    fi
fi


# Check if sourcing conda in .bashrc

if [[  $(__test_bashrc) = false ]]; then
    warning "** Conda was not detected on your PATH. **"
    warning "This is normal if you have not installed Conda before."
    warning "To add it to your current .bashrc to be sourced during shell"
    warning "startup, run:"
    warning "   'echo \"# Added to activate conda within shell\" >> ~/.bashrc"
    warning "   'echo \"source ${__conda_path}/etc/profile.d/conda.sh\" >> ~/.bashrc"
    warning "and close and re-open your terminal session to apply."
    warning "When finished, run 'conda activate ${__tsap_env}' to begin."
else
    info "Done. Run 'conda activate ${__tsap_env}' to begin."
fi
