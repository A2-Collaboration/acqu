export CCCOMP = g++
export CCOMP  = gcc
SHELL = /bin/sh
#ARCH := $(shell uname -m)
#OS := $(subst -,,$(shell uname))

#define some environment variables used for compiling
acqu_dir := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
export acqu_sys := $(acqu_dir)/acqu_core
export acqu := $(acqu_dir)/acqu_user
export CALIB := $(acqu_dir)/CaLib
export OSCAR := $(acqu_dir)/OSCAR
export Worker := $(acqu_dir)/Worker/AR
export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):$(CALIB)/lib:$(OSCAR)/lib

all: check acqu_core calib worker acqu_user
	@echo
	@echo
	@echo "#########################################################################"
	@echo
	@echo "SUCCESSFULLY COMPILED"
	@echo
	@echo "If you're using bash, you may want to add:"
	@echo "'source $(acqu_dir)/setup.sh' to your ~/.bashrc file"
	@echo
	@echo "Or, if you're using tcsh, you may want to add:" 
	@echo "'source $(acqu_dir)/setup.csh' to your ~/.tcshrc file"
	@echo
	@echo "This automatically sets the appropriate environment variables to run acqu"
	@echo
	@echo "#########################################################################"
	@echo
	@echo 


check:
#	check for ncurses
	@if [ ! -f /usr/include/ncurses.h ]; then echo "ncurses.h not found! Please install libncurses dev package."; exit 1; fi
#	ensure that you execute this script in a ROOT enabled shell
	@if [ "x$(ROOTSYS)" = "x" ]; then echo "ERROR: \$$ROOTSYS not found, shell not ROOT enabled?"; exit 1; fi
#	check if root is build with mysql support
	@if [ -z "$(filter mysql,$(shell root-config --features))" ]; then echo "ROOT built without MySQL support! Please install mysql(client) dev package and make ROOT again."; exit 1; fi #filter ignores matching strings in the middle of others in contrast to findstring which searches for a matching string everywhere
#	Check if the GNU gold linker is installed on Ubuntu systems newer than 11.04. It is used for a successful compilation of the Calibration code due to some lib problems with the GNU linker in the newer versions
	@if [ $(shell echo $(shell (sed '/DISTRIB_RELEASE/!d' /etc/lsb-release 2>/dev/null || echo 0) | sed 's/DISTRIB_RELEASE=//') '>' 11.04 | bc -l) -eq 1 ]; then if [ -z "$(filter gold,$(shell ld -v))" ]; then echo "Please install the binutils-gold package to get the GNU gold linker, otherwise the compilation fails on Ubuntu > 11.04."; exit 1; fi; fi

acqu_DAQ: check
	$(MAKE) -C $(acqu_sys) AcquDAQ

acqu_core: check
	$(MAKE) -C $(acqu_sys) AcquRoot
	$(MAKE) -C $(acqu_sys) AcquMC

acqu_user: check
	$(MAKE) -C $(acqu) AcquRoot

calib: check
	$(MAKE) -C $(OSCAR)
	$(MAKE) -C $(CALIB)

worker:
	$(MAKE) -C $(Worker)

clean:
	$(MAKE) -C $(acqu_sys) clean
	$(MAKE) -C $(OSCAR) clean
	$(MAKE) -C $(CALIB) clean
	$(MAKE) -C $(acqu) clean
	$(MAKE) -C $(Worker) clean

ultraclean: clean
	$(MAKE) -C $(acqu_sys) ultraclean
	$(MAKE) -C $(acqu) ultraclean

