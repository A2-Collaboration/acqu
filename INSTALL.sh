#!/bin/sh

error_exit()
{
	echo "$1" 1>&2
	exit 1
}

#check for ncurses
if [ ! -f /usr/include/ncurses.h ]; then
	error_exit "ncurses.h not found! Please install libncurses dev package."
fi
# ensure that you execute this script in a ROOT enabled shell
# something like this could do the job
# cd /opt/root && . ./bin/thisroot.sh && cd - > /dev/null
if [ "x$ROOTSYS" = "x" ]; then
    error_exit "ERROR: ROOTSYS not found, shell not ROOT enabled?"
fi
#check if root is build with mysql support
root-config --features | grep mysql > /dev/null
if [ $? -ne 0 ]; then
	error_exit "ROOT built without MySQL support! Please install mysql(client) dev package and make ROOT again."
fi

#get current acqu location (relative path)
d=${0%/*} #$(dirname $0)

# source the environment "usual" variables
. $d/setup.sh

# 1a) Build AcquRoot in acqu_core
cd $acqu_sys
make AcquRoot || error_exit "Cannot compile acqu_core, check previous output" 

# 1b) Build OSCAR
cd $OSCAR
make || error_exit "Cannot compile OSCAR, check previous output" 

# 1c) Build CaLib
cd $CALIB
make || error_exit "Cannot compile CaLib, check previous output" 

# 2) Build minimal acqu_user code
cd $acqu
make AcquRoot || error_exit "Cannot compile acqu_user, check previous output" 

echo
echo
echo "#########################################################################"
echo
echo "SUCCESSFULLY COMPILED"
echo "You might want to add 'source $acqu_dir/setup.sh' to your ~/.bashrc"
echo "This automatically sets the appropriate environment variables to run acqu"
echo

