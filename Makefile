#CCCOMP = g++
#CCOMP  = gcc
#OS := $(subst -,,$(shell uname))
#ARCH = $(shell uname -m)

all:
	$(MAKE) -C $(acqu_sys) AcquRoot
	$(MAKE) -C $(OSCAR)
	$(MAKE) -C $(CALIB)
	$(MAKE) -C $(acqu) AcquRoot

clean:
	$(MAKE) -C $(acqu_sys) clean
	$(MAKE) -C $(OSCAR) clean
	$(MAKE) -C $(CALIB) clean
	$(MAKE) -C $(acqu) clean

