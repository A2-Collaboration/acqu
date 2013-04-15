#!/bin/sh
export acqu_dir=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
# echo $acqu_dir
export acqu_sys=$acqu_dir/acqu_core
export acqu=$acqu_dir/acqu_user
export CALIB=$acqu_dir/CaLib
export OSCAR=$acqu_dir/OSCAR
export Worker=$acqu_dir/Worker/AR
# use a widespread compiler
# only make in acqu_core and acqu_user need this
export CCCOMP=g++
# acqu_dir is not needed
#acqu_dir=

# add some library path
LD_LIBRARY_PATH=$CALIB/lib:$OSCAR/lib:$LD_LIBRARY_PATH
# add AcquRoot and Worker Dictory to $PATH
PATH=$acqu/bin:$Worker/bin:$PATH
