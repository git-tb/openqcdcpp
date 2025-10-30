# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/../..")
  file(MAKE_DIRECTORY "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/../..")
endif()
file(MAKE_DIRECTORY
  "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/build/_deps/openqcdpplib-build"
  "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/build/_deps/openqcdpplib-subbuild/openqcdpplib-populate-prefix"
  "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/build/_deps/openqcdpplib-subbuild/openqcdpplib-populate-prefix/tmp"
  "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/build/_deps/openqcdpplib-subbuild/openqcdpplib-populate-prefix/src/openqcdpplib-populate-stamp"
  "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/build/_deps/openqcdpplib-subbuild/openqcdpplib-populate-prefix/src"
  "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/build/_deps/openqcdpplib-subbuild/openqcdpplib-populate-prefix/src/openqcdpplib-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/build/_deps/openqcdpplib-subbuild/openqcdpplib-populate-prefix/src/openqcdpplib-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/tobiasb/OneDrive/projects/PDFcode/openqcd++/tests/F2perf/build/_deps/openqcdpplib-subbuild/openqcdpplib-populate-prefix/src/openqcdpplib-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
