#!/bin/bash
set -e

# get the number of *physical* cores to use as the number of threads
CORES=$(getconf _NPROCESSORS_ONLN)
echo "Running with $CORES threads."

# get the path of this script, since we want to resolve the sysimage and run_cminject.jl
# relative to it
SCRIPTPATH="${BASH_SOURCE%/*}"
SYSIMAGE="$SCRIPTPATH/build/CMInject.so"

if [ ! -f "$SYSIMAGE" ]; then
    echo "The sysimage could not be found at $SYSIMAGE -- please run build/create_sysimage.sh."
    exit 1
fi

# then start an interactive Julia session with the sysimage, using $CORES threads and declaring
# this script's directory as the project
julia -t "$CORES"\
    -J"$SYSIMAGE" --sysimage-native-code yes\
    --project="$SCRIPTPATH" "$SCRIPTPATH/cminject.jl"\
    $@
