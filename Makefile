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

EXE 	:= NO_EXECUTBALE_SPECIFIED
BIN     := bin
SRC     := src
LIBSRC	:= libsrc


ARGS	:=

INCLUDE :=
DIRLIBRARIES := 
LIBRARIES :=
DYNLINK :=

INCLUDE += -Iinclude
INCLUDE += -I/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/include

DIRLIBRARIES += -L/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib
DIRLIBRARIES += -L/home/tobiasb/OneDrive/projects/PDFcode/openqcdrad-2.1/mycode/lib

LIBRARIES += -lgsl
LIBRARIES += -lLHAPDF
LIBRARIES += -lboost_program_options
LIBRARIES += -lmyqcdlib
LIBRARIES += -lgfortran

# DYNLINK += -Wl,-rpath=/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib

build: $(wildcard ./$(LIBSRC)/*.cpp)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $(SRC)/$(EXE).cpp $^ $(DIRLIBRARIES) $(LIBRARIES) $(DYNLINK) -o $(BIN)/$(EXE)

checklibs:
	ldd ./bin/$(EXE)


