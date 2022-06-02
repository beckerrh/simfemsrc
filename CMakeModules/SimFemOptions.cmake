# ==============================================================================
#                              Build options
# ==============================================================================

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

OPTION(BUILD_SHARED_LIBS "Build dynamic libraries" ON)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
	add_definitions(-DCLANG)
    SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -stdlib=libc++ -std=gnu++11" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
	add_definitions(-DGNU)
endif()

#----------- build-type --------------
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#  add_definitions(-DARMA_EXTRA_DEBUG)
	add_definitions(-DDEBUG)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(-DARMA_NO_DEBUG)
endif()

#----------- matrix-storage --------------
#add_definitions(-DSPARSITYPATTERNLONG)

# SET(SIMFEMCPPLIBRARY ${SIMFEM_INSTALL_DIR}/lib/libSimFem${CMAKE_BUILD_TYPE}.dylib)
