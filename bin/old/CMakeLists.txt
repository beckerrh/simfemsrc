SET(EXECUTABLE simfemcpp)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/libcpp ${SIMFEM_EXTERNAL_INCLUDES})
ADD_EXECUTABLE(${EXECUTABLE} main)
TARGET_LINK_LIBRARIES(${EXECUTABLE} ${SIMFEM_EXTERNAL_LINKLIBS} simfem)
INSTALL(TARGETS ${EXECUTABLE} DESTINATION bin)
