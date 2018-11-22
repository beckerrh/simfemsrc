#-------------------------------------------#
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(SimFemLib)
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../CMakeModules)
INCLUDE(SimFemOptions)
INCLUDE(SimFemFindExternalLibraries)

set(LIBRARYCPP SimFem${CMAKE_BUILD_TYPE})
set(LIBRARYCPPY simfempy)

ADD_SUBDIRECTORY(libcpp)
message("{WITHOUTPYTHON} ${WITHOUTPYTHON}")
if(NOT ${WITHOUTPYTHON})
    ADD_SUBDIRECTORY(libcppy)
endif()
