# ==============================================================================
#                              Find external libraries
# ==============================================================================

find_package(Armadillo REQUIRED)
find_package(MPI  COMPONENTS CXX REQUIRED)
#set(MPI_CXX_COMPILER mpicxx)

FIND_LIBRARY(UMFPACK_LIBRARIES NAMES umfpack libumfpack UMFPACK libUMFPACK
               PATH_SUFFIXES ${CMAKE_LIBRARY_ARCHITECTURE}
               HINTS ${CMAKE_INSTALL_PREFIX}/External/lib)
if(NOT UMFPACK_LIBRARIES)
    message(FATAL_ERROR "*** umfpack library not find (consider installing suite-sparse)")
endif()
SET(SIMFEM_EXTERNAL_LINKLIBS ${ARMADILLO_LIBRARIES} ${UMFPACK_LIBRARIES} ${MPI_LIBRARIES})
SET(SIMFEM_EXTERNAL_INCLUDES ${ARMADILLO_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
