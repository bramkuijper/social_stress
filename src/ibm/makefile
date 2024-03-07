# make the executable
EXE=stress_social.exe
CPP=stress_social.cpp main.cpp individual.cpp patch.cpp
HPP=stress_social.hpp parameters.hpp individual.hpp patch.hpp
CXX=g++
CXXFLAGS=-Wall -std=c++20 -ggdb

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP) 

clean :
	rm -rf $(EXE)
