cmake_minimum_required(VERSION 3.23)
project(f2heavyintegrand)

add_compile_options("-g;-Wall;-Wno-sign-compare;-Wpedantic;-Wconversion")

include(FetchContent)
FetchContent_Declare(
    openQCDpplib
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..
)
FetchContent_MakeAvailable(openQCDpplib)

include_directories(
	/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/include
	/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/include
)

add_executable(test)

set_target_properties(test PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)


target_sources(test
	PRIVATE
		test.cpp
)

target_link_libraries(test
	PRIVATE
		openQCDpplib
		/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib/libLHAPDF.so
		/home/tobiasb/OneDrive/projects/PDFcode/openqcdrad-2.1/mycode/lib/libmyqcdlib.a
)