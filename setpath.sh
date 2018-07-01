# Bash command to set ALTA environment for ease of use.
# Invoke this script one GNU/Linux or OSX using the source
# command:
#
#    source setpath.sh
#
# ALTA command line programs and plugins will be available
# to the shell.
#
sources=`pwd`

if [ "$BASH_VERSION" ]; then
        sources=`dirname $BASH_SOURCE`
        sources=`builtin cd ${sources}; builtin pwd`
elif [ "$ZSH_VERSION" ]; then
        sources=$(dirname "$0:A")
fi

builddir="$sources/build"

external=${sources}/external/build/lib
scripts=${sources}/sources/scripts/

export ALTA=$sources
export ALTA_DIR=$sources/sources
export ALTA_PLUGIN_PATH="$builddir/plugins"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$builddir/plugins:$external"
export PATH="$PATH:$builddir/softs:$scripts"
export PYTHONPATH="$PYTHONPATH:$builddir/python"
export JUPYTER_PATH="$JUPYTER_PATH:$builddir/python"
