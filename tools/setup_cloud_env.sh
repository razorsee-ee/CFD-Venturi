#!/bin/bash
# One-time setup of a fresh Claude cloud container for this project:
#   - system packages for building OpenFOAM-12 (Foundation)
#   - OpenFOAM-12 via tools/setup_openfoam12.sh: binary package (minutes) or
#     source build in /opt (~1 h on 4 cores)
#   - Python packages: numpy, matplotlib (plots), cadquery (geometry -> STEP/STL)
#   - /opt/of12env.sh: source it to get the OpenFOAM environment
#
# Idempotent: every step is skipped when it is already done, so it can be
# re-run after an interruption. Usage (as root, from the repository root):
#
#   bash tools/setup_cloud_env.sh            # everything
#   bash tools/setup_cloud_env.sh --no-foam  # only the Python packages
#
# Run it in the background and wait (minutes to ~1 h):
#   nohup bash tools/setup_cloud_env.sh > /tmp/setup_cloud_env.log 2>&1 &
#
# If dl.openfoam.org is allowed in the environment's network settings, the
# binary package is much faster and is tried first (see setup_openfoam12.sh).
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
# package or source build, /opt/of12env.sh and a basic check
bash "$(dirname "$0")/setup_openfoam12.sh" || exit 1
ENV_FILE=/opt/of12env.sh

# ------------------------------------------------------------------ check
bash -c "source $ENV_FILE && for t in foamRun blockMesh snappyHexMesh surfaceFeatures checkMesh decomposePar reconstructPar foamToVTK foamPostProcess foamListTimes foamDictionary; do command -v \$t > /dev/null || echo MISSING \$t; done; ls \$FOAM_LIBBIN | grep -q lagrangianParcel || echo MISSING liblagrangianParcel" \
    && log "check finished (no MISSING lines above = ready)"
