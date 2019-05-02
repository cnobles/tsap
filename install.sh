#!/usr/bin/env bash

__conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

read -r -d '' __usage <<-'EOF'
  -e --environment  [arg] Environment to install to. Default: "vivi"
  -s --vivi_dir   [arg] Location of vivi source code. Default: this directory
  -c --conda  [arg]       Location of Conda installation. Default: ${PREFIX}
  -u --update [arg]       Update vivi [lib]rary, conda [env], or [all].
  -r --requirements       Install from requirements rather than build (slow).
  -t --test               After installation, run test to check functionality.
  -v --verbose            Show subcommand output
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

read -r -d '' __helptext <<-'EOF'
 This script installs or upgrades VivI, including Conda (if not installed).
 To upgrade, pass the '--upgrade all' option, then be sure to update your config
 files using 'vivi config update'.
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
__vivi_dir="${arg_s:-$(readlink -f ${__dir})}"
__vivi_env="${arg_e:-vivi}"
__run_vivi_tests=false
__reqs_install=false
__update_lib=false
__update_env=false
__req_r_version="3.4.1"
__old_path=$PATH
__output=${2-/dev/stdout}

PATH=$PATH:${__conda_path}/bin

if [[ "${arg_t:?}" = "1" ]]; then
    __run_vivi_tests=true
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
        | grep -Fxq $__vivi_env > /dev/null) && \
        echo true || echo false
    else
      	echo false
    fi
}

function __test_r_version () {
    activate_vivi

    local sem_version=$(R --version | grep 'R version' | cut -d ' ' -f 3)

    if (( ${sem_version//./} >= ${__req_r_version//./} )); then
        echo true
    else
        echo false
    fi

    deactivate_vivi
}

function __test_r_packages () {
    activate_vivi

    $(Rscript ${__vivi_dir}/tools/rscripts/check_for_required_packages.R \
        > /dev/null) && echo true || echo false

    deactivate_vivi
}

function __test_vivilib() {
    if [[ $(__test_env) = true ]]; then
      	activate_vivi
      	command -v vivi &> /dev/null && echo true || echo false
      	deactivate_vivi
    else
      	echo false
    fi
}

function __test_vivi() {
    if [[ $(__test_env) = true ]]; then
      	$(bash ${__vivi_dir}/etc/tests/test.sh ${__vivi_env} &> /dev/null) && \
      	    echo true || echo false
    else
      	echo "fail"
    fi
}

function __test_bashrc () {
  grep conda.sh ~/.bashrc > /dev/null && echo true || echo false
}

function activate_vivi () {
    set +o nounset
    conda activate $__vivi_env
    set -o nounset
}

function deactivate_vivi () {
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
        debug_capture conda env update --name=$__vivi_env ${install_options} 2>&1
    else
        local install_options="--quiet --yes --file etc/build.b0.2.0.txt"
        debug_capture conda create --name=$__vivi_env ${install_options} 2>&1
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
    activate_vivi

    echo -ne "#/bin/sh\nexport VIVI_DIR=${__vivi_dir}" > \
	      ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

    mkdir -p ${CONDA_PREFIX}/etc/conda/deactivate.d/

    echo -ne "#/bin/sh\nunset VIVI_DIR" > \
	      ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh

	  deactivate_vivi
}

function install_vivilib () {
    activate_vivi

    debug_capture pip install --upgrade ${__vivi_dir}/tools/vivilib/ 2>&1

    if [[ $(__test_vivilib) != true ]]; then
      	installation_error "Library installation"
    fi

    deactivate_vivi
}

info "Starting VivI installation..."
info "    Conda path:  ${__conda_path}"
info "    VivI src:  ${__vivi_dir}"
info "    VivI env:  '${__vivi_env}'"

debug "Components detected:"
__conda_installed=$(__test_conda)
debug "    Conda:         ${__conda_installed}"
__env_exists=$(__test_env)
debug "    Environment:   ${__env_exists}"
__vivilib_installed=$(__test_vivilib)
debug "    Library:       ${__vivilib_installed}"

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

# Create Conda environment for VivI
if [[ $__env_exists = true && $__update_env = false ]]; then
    info "Specified environment already exists (use '--update env' to update)"
else
    if [[ $__reqs_install = "true" ]]; then
        __build_source="etc/requirements.yml"
    else
        __build_source="etc/build.b0.2.0.txt"
    fi

    info "Creating VivI environment..."
    info "    Building from: $__build_source"
    install_environment
    __env_changed=true
    info "$__vivi_env environment created."
fi


# Always update the env_vars.sh in the VivI environment
debug "Updating \$VIVI_DIR variable to point to ${__vivi_dir}"
info "Setting environmental variables..."
install_env_vars


# Install vivilib into environment if changed or requested
if [[ $__env_changed = true ]]; then
    info "Environment installed/updated; (re)installing VivI library..."
    install_vivilib
elif [[ $__vivilib_installed = false ]]; then
    info "Installing VivI library..."
    install_vivilib
elif [[ $__update_lib = true ]]; then
    info "Updating VivI library..."
    install_vivilib
else
    info "VivI library already installed (use '--update lib' to update)"
fi


# Run tests if requested
if [[ $__run_vivi_tests = true ]]; then
    info "Running VivI tests...(this may take a few mins)"
    __vivi_tested=$(__test_vivi)
    
    if [[ $__vivi_tested = true ]]; then
        info "    VivI Tests:  passed"
    else
        warning "    VivI Tests:  FAILED"
        warning "    Try running the test outside of the install to confirm."
        warning "    Just run 'bash etc/tests/test.sh $__vivi_env'."
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
    warning "When finished, run 'conda activate ${__vivi_env}' to begin."
else
    info "Done. Run 'conda activate ${__vivi_env}' to begin."
fi
