CXX = g++
CXXFLAGS = -O0

LIBS = -lnlopt

run_opt: check-para optimization.o model.o $(parameter).o
	$(CXX) $(CXXFLAGS)  $(filter-out $<,$^) $(LIBS)  -o $@

%.o : %.c 
	$(CXX) -c $(CXXFLAGS) $< -o $@

check-para : 
ifndef parameter
	$(error Parameter is not set)
endif

