cmake_minimum_required(VERSION 3.23)
project(integrationtest)

add_compile_options("-g;-Wall;-Wno-sign-compare;-Wpedantic;-Wconversion")

add_executable(test)

include(FetchContent)
FetchContent_Declare(
    openQCDpplib
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..
)
FetchContent_MakeAvailable(openQCDpplib)

set_target_properties(test PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

include_directories(
	/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/include
)

target_sources(test
	PRIVATE
		test.cpp
)

target_link_libraries(test
	PRIVATE
		openQCDpplib
		/home/tobiasb/OneDrive/projects/PDFcode/openqcdrad-2.1/mycode/lib/libmyqcdlib.a
		gsl
)