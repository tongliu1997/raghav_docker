os = $(shell uname -s)

INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(STARPICOPATH) -I${JETREADER}/include
#INCFLAGS      += -I$(JAPATH)/src
INCFLAGS      += -I./src


ifeq ($(os),Linux)
CXXFLAGS      = -O2 -fPIC -pipe -Wall -std=c++1z
CXXFLAGS     += -Wno-unused-variable
CXXFLAGS     += -Wno-unused-but-set-variable
CXXFLAGS     += -Wno-sign-compare
CXXFLAGS     += -D_VANILLA_ROOT_ # for star library 
# # for gprof -- cannot combine with -g!
# CXXFLAGS     += -pg
# # for valgrind, gdb
# CXXFLAGS     += -g
else
CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
CXXFLAGS     += -Wno-return-type-c-linkage 
CXXFLAGS     += -D_VANILLA_ROOT_ # for star library 

## for debugging:
# CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       =
# # for gprof -- cannot combine with -g!
# LDFLAGS      += -pg
# # for valgrind, gdb
# LDFLAGS      += -g

LDFLAGSS      = --shared 
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

LDFLAGS += -lEG

# # uncomment for debug info in the library
# CXXFLAGS     += -g


ROOTLIBS      = $(shell root-config --libs)

LIBPATH       = $(ROOTLIBS) -L$(FASTJETDIR)/lib -L$(STARPICOPATH) -L${JETREADER}/lib
LIBPATH       += -L$(FASTJETDIR)/lib -L$(STARPICOPATH) -L${JETREADER}/lib
#LIBS          = -lMyJetlib	
LIBS         = -lfastjet -lfastjettools  -lTStarJetPico -lRecursiveTools -ljetreader -lStPicoEvent
LIBS         += -lConstituentSubtractor
#LIBS         += -lIterativeConstituentSubtractor
LIBS         += -lNsubjettiness
LIBS         += -lEnergyCorrelator
LIBS         += -lyaml-cpp

## fun with pythia :-/
## make is a horrible horrible tool. Do not touch these lines, any whitespace will make it break
#dummy := "$(shell find $(PYTHIA8DIR)/lib/ -name liblhapdfdummy\*)"
#ifneq ("",$(dummy))
#LIBS         += -llhapdfdummy
#endif

# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################
#INCS = $(SDIR)/JetAnalyzer.hh $(SDIR)/AjParameters.hh
#INCS = $(SDIR)/JetShowerParams.hh $(SDIR)/ShowerAnalyzer_Data.hh $(SDIR)/TelescopingJets.hh $(SDIR)/JetAnalyzer.hh $(SDIR)/centrality_run14.hh $(SDIR)/run14_eff.hh
#INCS = $(SDIR)/JetShowerParams.hh $(SDIR)/ShowerAnalyzer_Data_Y7.hh $(SDIR)/TelescopingJets.hh $(SDIR)/JetAnalyzer.hh $(SDIR)/centrality_run14.hh $(SDIR)/run14_eff.hh
#INCS = $(SDIR)/JetShowerParams.hh  $(SDIR)/JetAnalyzer.hh $(SDIR)/TelescopingJets.hh  $(SDIR)/ShowerAnalyzer_ppData_Embed_AuAuMinBias.hh $(SDIR)/centrality_run14.hh $(SDIR)/run14_eff.hh
#INCS = $(SDIR)/JetAnalyzer.hh $(SDIR)/JetShowerParams.hh $(SDIR)/TelescopingJets.hh  $(SDIR)/centrality_run14.hh  $(SDIR)/run14_eff.hh
#INCS = ${SDIR}/TelescopingJets.hh $(SDIR)/JetShowerParams.hh $(SDIR)/run14_eff.hh ${SDIR}/centrality_def.hh ${SDIR}/centrality_run14.hh 

###############################################################################
# standard rules $(OBJ): $(ODIR)/%.o: src/%.c $(DEPS)
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $^ $(LIBS) -o $@

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all     : ${BDIR}/dataset_QA
#	 doxy

${BDIR}/dataset_QA  : ${ODIR}/dataset_QA.o 

###############################################################################
##################################### MISC ####################################
###############################################################################


doxy: html/index.html

html/index.html : $(INCS) src/* Doxyfile
#	doxygen
	@echo 
	@echo Updating documentation
	( cat Doxyfile ; echo "QUIET=YES" ) | doxygen -

clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -rvf $(BDIR)/*dSYM
	rm -rvf lib/*dSYM	
	rm -vf $(BDIR)/*
	rm -vf lib/*
	rm -vf $(SDIR)/dict.cxx $(SDIR)/dict.h

.PHONY : clean doxy
