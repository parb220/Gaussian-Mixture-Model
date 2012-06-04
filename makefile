SOURCES = gaussian_mixture_equal_energy.cpp MultiDimSample.cpp MultiDimSampleChain.cpp EnergySet.cpp GaussianModel.cpp GaussianMixtureModel.cpp UniformModel.cpp
OBJS = $(SOURCES: .cpp = .o)
CPP = g++
DEBUG = -g
CPPFLAGS = -c -Wall $(DEBUG)
LINKFLAGS = -Wall $(DEBUG)
GSL_LINKFLAGS = -lgsl -lgslcblas -lm 
EXECUTABLE = gaussian_mixture_equal_energy
 
all : $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJS) 
	$(CPP) $(LINKFLAGS) $(OBJS) $(GSL_LINKFLAGS) -o $@

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $< -o $@ 

