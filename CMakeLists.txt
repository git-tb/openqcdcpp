cmake_minimum_required(VERSION 3.23)

project(openQCDpp)

add_compile_options("-g;-Wall;-Wno-sign-compare;-Wpedantic;-Wconversion")

include_directories(
	include	
	/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/include
	/home/tobiasb/OneDrive/projects/PDFcode/hoppet/share/include
)

########################################################################################################################
########################################################################################################################

add_library(openQCDpplib)

target_link_libraries(openQCDpplib
	PRIVATE
		/home/tobiasb/OneDrive/projects/PDFcode/openqcdrad-2.1/mycode/lib/libmyqcdlib.a
		/home/tobiasb/OneDrive/projects/PDFcode/chaplin-1.2/lib/libchaplin.so
		gfortran
		gsl
)

target_sources(openQCDpplib
	PRIVATE
		libsrc/coefffunc.cpp
		libsrc/constants.cpp
		libsrc/flags.cpp
		libsrc/integrate.cpp
		libsrc/logger.cpp
		libsrc/pdf.cpp
		libsrc/structfunc.cpp
		libsrc/structfunc.cpp
		libsrc/hqcoeff.cpp

	PUBLIC
		FILE_SET HEADERS
		BASE_DIRS
    		include
		FILES
			include/coefffunc.h
			include/constants.h
			include/flags.h
			include/integrate.h
			include/logger.h
			include/pdf.h
			include/structfunc.h
			include/hqcoeff.h
)




########################################################################################################################
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)


########################################################################################################################
########################################################################################################################

add_executable(lhapdftest)

target_sources(lhapdftest
	PRIVATE
		src/lhapdftest.cpp
)

target_link_libraries(lhapdftest
	PRIVATE
		/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib/libLHAPDF.so
		boost_program_options
)

########################################################################################################################
########################################################################################################################

add_executable(hoppettest)

target_sources(hoppettest
	PRIVATE
		src/hoppettest.cpp
)

target_link_libraries(hoppettest
	PRIVATE
		/home/tobiasb/OneDrive/projects/PDFcode/hoppet/share/lib/libhoppet.a
		/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib/libLHAPDF.so
		gfortran
)

########################################################################################################################
########################################################################################################################

add_executable(teststructfunc)

target_sources(teststructfunc
	PRIVATE
		src/teststructfunc.cpp
)

target_link_libraries(teststructfunc
	PRIVATE
		openQCDpplib
		/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib/libLHAPDF.so
)

########################################################################################################################
########################################################################################################################

add_executable(testF2heavyintegrand)

target_sources(testF2heavyintegrand
	PRIVATE
		src/testF2heavyintegrand.cpp
)

target_link_libraries(testF2heavyintegrand
	PRIVATE
		openQCDpplib
		/home/tobiasb/OneDrive/projects/PDFcode/LHAPDF-6.5.5/build/lib/libLHAPDF.so
)


endif()
########################################################################################################################