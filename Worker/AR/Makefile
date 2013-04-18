#Architecture
ARCH = x86_64
OPTI = -O3 #-march=opteron

#ROOT includes
IR = $(ROOTSYS)/include

#ROOT libraries
ROOTLIBS   := $(shell root-config --glibs) -lEG -lstdc++
ROOTCFLAGS := $(shell root-config --auxcflags)

#Directories
WORKER := $(shell pwd)
O  = $(WORKER)/obj
B  = $(WORKER)/bin
L  = $(WORKER)/lib
MC = $(WORKER)/src
D  = $(WORKER)/deps

#Compilers
CCCOMP = g++
CCOMP  = gcc

#Compilation flags
CXXFLAGS  = -fPIC $(ROOTCFLAGS) -w $(OPTI) -c
CFLAGS    = -fPIC $(ROOTCFLAGS) -w $(OPTI)
LDXXFLAGS = -shared -o $L/$@
INCFLAGS  = -I$(IR)

CLEANMOD = $B/Worker

#Build Phobos event generator
all: Worker

#Clean up and remove object files
clean:
	@rm -f $(D)/*.d
	@rm -f $(O)/*.o
	@rm -f $(CLEANMOD)
	@rm -f $(L)/lib*.so
	@rm -f $(MC)/WorkerDict.cc
	@echo Binaries and dependencies deleted
	@echo "--------------------------------------------------------------------------------"

#Clean up and remove additional object files
veryclean:
	@rm -f $(D)/*.d
	@rm -f $(O)/*.o
	@rm -f $(O)/$(ARCH)/*.o
	@rm -f $(CLEANMOD)
	@rm -f $(L)/lib*.so
	@rm -f $(MC)/WorkerDict.cc
	@echo Binaries, libraries and dependencies deleted
	@echo "--------------------------------------------------------------------------------"

vpath %.cc $(MC)
vpath %.c $(MC)
vpath %.h $(MC)
vpath %.o $(O)
vpath %.d $(D)

#Phobos Classes
MCCl   = $(wildcard $(MC)/T*.cc)
MCSrc  = $(MCCl) $(MC)/WorkerDict.cc
MCInc  = $(patsubst $(MC)/%.cc, $(MC)/%.h, $(MCCl))
MCObjD = $(patsubst $(MC)/%.cc, $(O)/%.o, $(MCSrc))
MCObj  = $(notdir $(MCObjD))
MCDep  = $(patsubst $(MC)/%.cc, $(D)/%.d, $(MCCl)) $(D)/Worker.d
INCFLAGS += -I$(MC)

#Main library for Phobos
libWorker.so: $(MCObj)
	@rm -f $(MC)/WorkerDict.cc
	$(CCCOMP) $(LDXXFLAGS) $(MCObjD)
	@echo "$@ done"
	@echo "--------------------------------------------------------------------------------"

#Phobos executable
Worker: Worker.o libWorker.so
	$(CCCOMP) $(O)/Worker.o $(L)/libWorker.so \
	$(ROOTLIBS) -lThread \
	-o $(B)/Worker
	@echo "$@ done"
	@echo "--------------------------------------------------------------------------------"


#Phobos ROOT dictionary
$(MC)/WorkerDict.cc: $(MCInc) WorkerLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint  -f $@ -c $(INCFLAGS) $^
	@echo "--------------------------------------------------------------------------------"

%.o: %.cc
	$(CCCOMP) $(CXXFLAGS) $(INCFLAGS) -o $(O)/$@ -c $<

$(O)/$(ARCH)/%.o: $(MC)/%.c
	$(CCOMP) $(OPTI) -c $< -o $@

$(D)/%.d: %.cc
	@echo Making dependencies for $<
	@$(SHELL) -ec 'gcc -MM $(INCFLAGS) -c $< \
		| sed '\''s%^.*\.ox%$*\.o%g'\'' \
		| sed '\''s%\($*\)\.o[ :]*%\1.o $@ : %g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

-include $(MCDep)
