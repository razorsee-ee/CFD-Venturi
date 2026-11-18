#!/bin/bash
# OpenFOAM-12 (Foundation) in a fresh Claude cloud container (Ubuntu, root).
# Project-independent: copy this file into any repository and run
#
#   nohup bash setup_openfoam12.sh > /tmp/setup_openfoam12.log 2>&1 &
#
# then wait for the line "check finished" and source /opt/of12env.sh.
#
# Order of attempts:
#   1) binary package openfoam12 (/opt/openfoam12, a few minutes). Needs
#      dl.openfoam.org in the environment's allowed hosts. dl.openfoam.org
#      redirects the .deb to http://downloads.sourceforge.net, which apt refuses
#      (https -> http), so the .deb is fetched from SourceForge over https and
#      checked against the MD5 from the signed apt index before installing.
#   2) build from GitHub sources (/opt/OpenFOAM-12, ~1 h on 4 cores).
#
# Idempotent: re-running skips what is already done.
set -u

log() { echo "[of12 $(date +%H:%M:%S)] $*"; }

PKG_BIN=/opt/openfoam12/platforms/linux64GccDPInt32Opt/bin/foamRun
SRC_BIN=/opt/OpenFOAM-12/platforms/linux64GccDPInt32Opt/bin/foamRun
ENV_FILE=/opt/of12env.sh
export DEBIAN_FRONTEND=noninteractive

if [ -x "$PKG_BIN" ] || [ -x "$SRC_BIN" ]; then
    log "OpenFOAM-12 already present"
else
    # ---------------------------------------------------------- 1) package
    if curl -fsS -m 20 -o /tmp/openfoam.gpg https://dl.openfoam.org/gpg.key 2>/dev/null; then
        log "dl.openfoam.org reachable: binary package"
        cp /tmp/openfoam.gpg /etc/apt/trusted.gpg.d/openfoam.asc
        . /etc/os-release
        echo "deb https://dl.openfoam.org/ubuntu ${VERSION_CODENAME} main" > /etc/apt/sources.list.d/openfoam.list
        apt-get update -qq > /tmp/apt_of12.log 2>&1
        line=$(apt-get install --print-uris -qq openfoam12 2>/dev/null | grep dl.openfoam.org | head -1)
        url=$(echo "$line" | awk '{print $1}' | tr -d "'")
        md5=$(echo "$line" | awk '{print $4}' | sed 's/^MD5Sum://')
        sf=${url/https:\/\/dl.openfoam.org\/ubuntu\//https://downloads.sourceforge.net/project/foam/foam/ubuntu/}
        deb=/tmp/$(basename "$url")
        if [ -n "$url" ] && curl -fsSL -m 1200 -o "$deb" "$sf" \
                && echo "$md5  $deb" | md5sum -c --quiet - ; then
            log "  $(basename "$deb"): MD5 OK, installing"
            apt-get install -y -qq "$deb" >> /tmp/apt_of12.log 2>&1 \
                && log "  installed in /opt/openfoam12" \
                || log "  apt-get install failed, see /tmp/apt_of12.log"
        else
            log "  could not fetch/verify the .deb, falling back to the source build"
        fi
    fi

    # ---------------------------------------------------------- 2) source
    if [ ! -x "$PKG_BIN" ]; then
        log "building OpenFOAM-12 from source (~1 h on 4 cores)"
        # a fresh container has empty package lists
        apt-get update -qq >> /tmp/apt_of12.log 2>&1
        apt-get install -y -qq \
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
        if [ -x "$SRC_BIN" ]; then
            log "OpenFOAM-12 built"
        else
            log "BUILD FAILED, see /opt/of12build.log (tail):"; tail -20 /opt/of12build.log; exit 1
        fi
    fi
fi

# ------------------------------------------------------------ environment
if [ -x "$PKG_BIN" ]; then ROOT=/opt/openfoam12; else ROOT=/opt/OpenFOAM-12; fi
echo "source $ROOT/etc/bashrc; export FOAM_SIGFPE=1; export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1" > $ENV_FILE
log "environment: source $ENV_FILE"

# ------------------------------------------------------------ check
bash -c "source $ENV_FILE
for t in foamRun blockMesh snappyHexMesh surfaceFeatures checkMesh setFields decomposePar reconstructPar \
         foamToVTK foamPostProcess foamListTimes foamDictionary mpirun; do
    command -v \$t > /dev/null || echo MISSING \$t
done
for l in lagrangianParcel multiphaseEuler incompressibleDenseParticleFluid multicomponentFluid; do
    ls \$FOAM_LIBBIN | grep -q \$l || echo MISSING lib\$l
done" && log "check finished (no MISSING lines above = ready)"
