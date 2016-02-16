os = $(shell uname -s)

INCFLAGS      = -I$(ROOTSYS)/include -I$(STARPICOPATH)
INCFLAGS      += -I./src
INCFLAGS      += -I./include
INCFLAGS      += -I${PPJET}/include
INCFLAGS      += -I${PPJET}/src


ifeq ($(os),Linux)
#ly CXXFLAGS      = 
CXXFLAGS      = -fPIC #ly	get complaint from complier and add it 
else
#CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
## for debugging:
 CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       = -g
LDFLAGSS      = -g --shared 
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = clang
endif


ROOTLIBS      = $(shell root-config --libs)

LIBPATH       = $(ROOTLIBS)  -L$(STARPICOPATH)
LIBS          = -lTStarJetPico


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################

###############################################################################
# standard rules
%.o : %.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

%  : %.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $(LIBS) $^ -o $@

###############################################################################
###############################################################################
############################# Main Targets ####################################
###############################################################################
all    : Ntrk  


###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo 
	@echo CLEANING
	rm -vf Ntrk.o
	rm -vf Ntrk 


