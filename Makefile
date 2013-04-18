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

all:
#	check for ncurses
	@if [ ! -f /usr/include/ncurses.h ]; then echo "ncurses.h not found! Please install libncurses dev package."; exit 1; fi
#	ensure that you execute this script in a ROOT enabled shell
	@if [ "x$(ROOTSYS)" = "x" ]; then echo "ERROR: \$$ROOTSYS not found, shell not ROOT enabled?"; exit 1; fi
#	check if root is build with mysql support
	@if [ -z "$(filter mysql,$(shell root-config --features))" ]; then echo "ROOT built without MySQL support! Please install mysql(client) dev package and make ROOT again."; exit 1; fi #filter ignores matching strings in the middle of others in contrast to findstring which searches for a matching string everywhere
	@echo Start compiling...
	@echo
	$(MAKE) -C $(acqu_sys) AcquRoot
	$(MAKE) -C $(OSCAR)
	$(MAKE) -C $(CALIB)
	$(MAKE) -C $(acqu) AcquRoot
	$(MAKE) -C $(Worker)
	@echo
	@echo
	@echo "#########################################################################"
	@echo
	@echo "SUCCESSFULLY COMPILED"
	@echo "You might want to add 'source $(acqu_dir)/setup.sh' to your ~/.bashrc"
	@echo "This automatically sets the appropriate environment variables to run acqu"
	@echo

#CaLib is needed to run properly
acqu:
	@if [ "x$(ROOTSYS)" = "x" ]; then echo "ERROR: \$$ROOTSYS not found, shell not ROOT enabled?"; exit 1; fi
	$(MAKE) -C $(acqu_sys) AcquRoot
	$(MAKE) -C $(acqu) AcquRoot

worker:
	$(MAKE) -C $(Worker)

user:
	$(MAKE) -C $(acqu) AcquRoot

clean:
	$(MAKE) -C $(acqu_sys) clean
	$(MAKE) -C $(OSCAR) clean
	$(MAKE) -C $(CALIB) clean
	$(MAKE) -C $(acqu) clean
	$(MAKE) -C $(Worker) clean

