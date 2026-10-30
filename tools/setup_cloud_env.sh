#!/bin/bash
# One-time setup of a fresh Claude cloud container for this project:
#   - system packages for building OpenFOAM-12 (Foundation)
#   - OpenFOAM-12 + ThirdParty-12 from GitHub, built in /opt (~2-3 h on 4 cores)
#   - Python packages: numpy, matplotlib (plots), cadquery (geometry -> STEP/STL)
#   - /opt/of12env.sh: source it to get the OpenFOAM environment
#
# Idempotent: every step is skipped when it is already done, so it can be
# re-run after an interruption. Usage (as root, from the repository root):
#
#   bash tools/setup_cloud_env.sh            # everything
#   bash tools/setup_cloud_env.sh --no-foam  # only the Python packages
#
# Run it in the background and wait (it takes hours):
#   nohup bash tools/setup_cloud_env.sh > /tmp/setup_cloud_env.log 2>&1 &
#
# If dl.openfoam.org is allowed in the environment's network settings, the
# binary package is much faster and is tried first (apt package openfoam12,
# installed to /opt/openfoam12).
set -u

log() { echo "[setup $(date +%H:%M:%S)] $*"; }

NO_FOAM=0
[ "${1:-}" = "--no-foam" ] && NO_FOAM=1

# ------------------------------------------------------------------ python
log "python packages"
for pkg in numpy matplotlib cadquery; do
    if python3 -c "import $pkg" > /dev/null 2>&1; then
        log "  $pkg: already installed"
    else
        log "  $pkg: installing"
        pip install -q "$pkg" 2>&1 | tail -2 || log "  $pkg: FAILED (pip install $pkg)"
    fi
done

[ "$NO_FOAM" = 1 ] && { log "done (--no-foam)"; exit 0; }

# ------------------------------------------------------------------ OpenFOAM
FOAM_SRC=/opt/OpenFOAM-12
FOAM_BIN=$FOAM_SRC/platforms/linux64GccDPInt32Opt/bin/foamRun
ENV_FILE=/opt/of12env.sh

if [ -x /opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/foamRun ] || [ -x "$FOAM_BIN" ]; then
    log "OpenFOAM-12 already present"
else
    # 1) binary package, if the network policy allows dl.openfoam.org
    if curl -fsS -m 20 -o /tmp/openfoam.gpg https://dl.openfoam.org/gpg.key 2>/dev/null; then
        log "dl.openfoam.org reachable: installing the binary package openfoam12"
        cp /tmp/openfoam.gpg /etc/apt/trusted.gpg.d/openfoam.asc
        . /etc/os-release
        echo "deb https://dl.openfoam.org/ubuntu ${VERSION_CODENAME} main" > /etc/apt/sources.list.d/openfoam.list
        apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq openfoam12 \
            && log "openfoam12 package installed in /opt/openfoam12"
    fi

    # 2) otherwise build from source
    if [ ! -x /opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/foamRun ]; then
        log "building OpenFOAM-12 from source (2-3 hours)"
        # a fresh container has empty package lists
        apt-get update -qq > /tmp/apt_of12.log 2>&1
        DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
            build-essential flex bison cmake zlib1g-dev libopenmpi-dev openmpi-bin \
            libscotch-dev libptscotch-dev libreadline-dev libncurses-dev libxt-dev \
            libboost-system-dev libboost-thread-dev libgmp-dev libmpfr-dev >> /tmp/apt_of12.log 2>&1 \
            || { log "apt-get failed, see /tmp/apt_of12.log"; exit 1; }
        cd /opt || exit 1
        [ -d OpenFOAM-12 ] || git clone -q --depth 1 https://github.com/OpenFOAM/OpenFOAM-12.git
        [ -d ThirdParty-12 ] || git clone -q --depth 1 https://github.com/OpenFOAM/ThirdParty-12.git
        # system scotch, no ParaView (nothing to display in the container)
        sed -i 's/^export SCOTCH_TYPE=ThirdParty/export SCOTCH_TYPE=system/; s/^export ParaView_TYPE=system/export ParaView_TYPE=none/' \
            OpenFOAM-12/etc/bashrc
        cd OpenFOAM-12 || exit 1
        bash -c "source /opt/OpenFOAM-12/etc/bashrc && ./Allwmake -j $(nproc) -q" > /opt/of12build.log 2>&1
        if [ -x "$FOAM_BIN" ]; then
            log "OpenFOAM-12 built"
        else
            log "BUILD FAILED, see /opt/of12build.log (tail):"; tail -20 /opt/of12build.log; exit 1
        fi
    fi
fi

# environment file used by the project notes: source /opt/of12env.sh
if [ -x /opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/foamRun ]; then
    echo 'source /opt/openfoam12/etc/bashrc; export FOAM_SIGFPE=1; export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1' > $ENV_FILE
else
    echo 'source /opt/OpenFOAM-12/etc/bashrc; export FOAM_SIGFPE=1; export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1' > $ENV_FILE
fi
log "environment: source $ENV_FILE"

# ------------------------------------------------------------------ check
bash -c "source $ENV_FILE && for t in foamRun blockMesh snappyHexMesh surfaceFeatures checkMesh decomposePar reconstructPar foamToVTK foamPostProcess foamListTimes foamDictionary; do command -v \$t > /dev/null || echo MISSING \$t; done; ls \$FOAM_LIBBIN | grep -q lagrangianParcel || echo MISSING liblagrangianParcel" \
    && log "check finished (no MISSING lines above = ready)"
