CXX       := g++
CXX_FLAGS := -std=c++2b
CXX_FLAGS += -g
# CXX_FLAGS += -O2
# CXX_FLAGS += -DNDEBUG
# CXX_FLAGS += -fopenmp
CXX_FLAGS += -Wall
CXX_FLAGS += -Wno-sign-compare

EXE 	:= NO_EXECUTBALE_SPECIFIED
BIN     := bin
SRC     := src

ARGS	:=

INCLUDE :=
DIRLIBRARIES := 
LIBRARIES :=
DYNLINK :=

INCLUDE += -Iinclude
INCLUDE += -I/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/include

DIRLIBRARIES += -L/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib

LIBRARIES += -lgsl
LIBRARIES += -lLHAPDF
LIBRARIES += -lboost_program_options

# DYNLINK += -Wl,-rpath=/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib

build:
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -o $(BIN)/$(EXE) $(SRC)/$(EXE).cpp $(DIRLIBRARIES) $(LIBRARIES) $(DYNLINK)

checklibs:
	ldd ./bin/$(EXE)


