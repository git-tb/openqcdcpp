cmake_minimum_required(VERSION 3.23)
project(tablelookup)

add_compile_options("-g;-Wall;-Wno-sign-compare;-Wpedantic;-Wconversion")

add_executable(test)

set_target_properties(test PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_sources(test
	PRIVATE
		test.cpp
)