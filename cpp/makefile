CXX = g++
LD = $(CXX)
OPTIM = -O0 #-march=native

bold := $(shell tput bold)
sgr0 := $(shell tput sgr0)

CXXFLAGS = -ggdb -pg 
CXXFLAGS += -Wall -pedantic -std=c++17 $(OPTIM) 

OBJDIR := build
SRCDIR := src
LIBDIR := lib
BINDIR := bin

SRCS := test.cpp
OBJS := $(SRCS:%.cpp=$(OBJDIR)/%.o)

PROGS := test
PROGBIN := $(PROGS:%.cpp=%)

INCLUDE= -I$(SRCDIR) \
         -I /home/corbin/uni/empirical-thermosphere-models/deps/eigen/ \
         -I ../src 

	 
LIB= -L ../lib -lSHTOOLS -lgfortran \
     -ldl -lfftw3 -lm -llapack -lblas \


all: $(OBJDIR) $(BINDIR) $(PROGBIN)

# places dependency files into a subdirectory named
DEPDIR := .deps

# These are the special GCC-specific flags which convince the compiler to generate the dependency file
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d

COMPILE.cpp = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(INCLUDE) -c -o $@

# Delete the built-in rules for building object files from .c files, so that our rule is used instead.
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp $(DEPDIR)/%.d | $(DEPDIR)
	$(info $(bold)$<$(sgr0))
	$(COMPILE.cpp) $<
	
		
test : $(OBJDIR)/test.o 
	$(info $(bold)linking$(sgr0))
	$(LD) -pg -o $(BINDIR)/$@ $< $(LIB)
		
# Declare a rule for creating the build directory if it doesn’t exist.
$(OBJDIR): ; @mkdir -p $@
	
# Declare a rule for creating the dependency directory if it doesn’t exist.
$(DEPDIR): ; @mkdir -p $@

# Declare a rule for creating binary directory if it doesn’t exist.
$(BINDIR): ; @mkdir -p $@

# Generate a list of all the dependency files that could exist.
DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d) $(PROGS:%.cpp=$(DEPDIR)/%.d)


# Mention each dependency file as a target, so that make won’t fail if the file doesn’t exist.
$(DEPFILES):

# Include the dependency files that exist. Use wildcard to avoid failing on non-existent files.
include $(wildcard $(DEPFILES))

clean :
	rm -rf $(OBJDIR)/*.o
	rm -rf $(LIBDIR)/*.a
	rm -rf $(PROGS:%=$(BINDIR)/%)

