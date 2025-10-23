CXX       := g++
CXX_FLAGS := -std=c++2b
CXX_FLAGS += -g
# CXX_FLAGS += -O2
# CXX_FLAGS += -DNDEBUG
# CXX_FLAGS += -fopenmp
CXX_FLAGS += -Wall
CXX_FLAGS += -Wno-sign-compare
CXX_FLAGS += -Wpedantic
CXX_FLAGS += -Wconversion		# int a = 0.1 <<< WARNING 

SRC     := /home/tobiasb/OneDrive/projects/PDFcode/openqcd++/src
LIBSRC	:= /home/tobiasb/OneDrive/projects/PDFcode/openqcd++/libsrc

INCLUDE :=
DIRLIBRARIES := 
LIBRARIES :=
DYNLINK :=

INCLUDE += -I/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/include
INCLUDE += -I/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/include

DIRLIBRARIES += -L/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/lib
DIRLIBRARIES += -L/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib
DIRLIBRARIES += -L/home/tobiasb/OneDrive/projects/PDFcode/openqcdrad-2.1/mycode/lib
DIRLIBRARIES += -L/home/tobiasb/OneDrive/projects/PDFcode/chaplin-1.2/lib

# LIBRARIES += -lboost_program_options
LIBRARIES += -lmycpplib
LIBRARIES += -lgsl
LIBRARIES += -lLHAPDF
LIBRARIES += -lmyqcdlib
LIBRARIES += -lchaplin
LIBRARIES += -lgfortran

DYNLINK += -Wl,-rpath=/home/tobiasb/OneDrive/projects/PDFcode/chaplin-1.2/lib

buildtest:
	$(CXX) $(CXX_FLAGS) $(INCLUDE) test.cpp $(DIRLIBRARIES) $(LIBRARIES) $(DYNLINK) -o test