#!/bin/sh

acqu_dir=$(dirname $(readlink -f $BASH_SOURCE)) #${BASH_SOURCE%/*} to get relative path

acqu_sys=$acqu_dir/acqu_core
acqu=$acqu_dir/acqu_user
CALIB=$acqu_dir/CaLib
OSCAR=$acqu_dir/OSCAR
Worker=$acqu_dir/Worker/AR
# add some paths
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CALIB/lib:$OSCAR/lib
PATH=$PATH:$acqu/bin:$Worker/bin

export acqu_sys acqu CALIB OSCAR Worker LD_LIBRARY_PATH PATH

