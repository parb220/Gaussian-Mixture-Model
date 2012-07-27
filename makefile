SOURCES = gaussian_mixture_equal_energy.cpp MultiDimSample.cpp MultiDimSampleChain.cpp EnergySet.cpp GaussianModel.cpp GaussianMixtureModel.cpp UniformModel.cpp
OBJS = $(SOURCES: .cpp = .o)
EXECUTABLE = gaussian_mixture_equal_energy

CPP = gcc
CPPFLAGS := $(CPPFLAGS) -g -Wall 
LINKFLAGS := $(LINKFLAGS) -lstdc++ -lgsl -lgslcblas -lm 
INCLUDE_DIR := $(INCLUDE_DIR) -I/usr/include/gsl
LINK_DIR := $(LINK_DIR) -L/usr/lib64 
 
all : $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJS) 
	$(CPP) $(CPPFLAGS) $(LINK_DIR) $(LINKFLAGS) $(OBJS) -o $@

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@ 

