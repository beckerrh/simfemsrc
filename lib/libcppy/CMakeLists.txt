#-------------------------------------------#
INCLUDE(SimFemWrapper)


SET(DIRECTORIES
# Alat
Mesh
Solvers
)

AUX_SOURCE_DIRECTORY(. SRC)
AUX_SOURCE_DIRECTORY(Tools SRC)
INCLUDE_DIRECTORIES(.)
get_filename_component(SMFEM_SOURCE_DIR ${CMAKE_SOURCE_DIR} PATH)
FOREACH(DIRECTORY ${DIRECTORIES})
  AUX_SOURCE_DIRECTORY(${DIRECTORY}Wrap SRC)
  INCLUDE_DIRECTORIES(${SMFEM_SOURCE_DIR}/lib/libcpp/${DIRECTORY})
ENDFOREACH(DIRECTORY)
INCLUDE_DIRECTORIES(${SIMFEM_EXTERNAL_INCLUDES})
INCLUDE_DIRECTORIES(${SMFEM_SOURCE_DIR}/lib/libcpp/Alat)
INCLUDE_DIRECTORIES(${SMFEM_SOURCE_DIR}/lib/libcpp/Perulangan)

INCLUDE(SimFemSummary)

ADD_LIBRARY(${LIBRARYCPPY} ${SRC})
set_target_properties(${LIBRARYCPPY} PROPERTIES PREFIX "" )
set_target_properties(${LIBRARYCPPY} PROPERTIES SUFFIX .so)
# set(CPPLIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libSimFem${CMAKE_BUILD_TYPE}.dylib)
TARGET_LINK_LIBRARIES(${LIBRARYCPPY} ${SIMFEM_EXTERNAL_LINKLIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${LIBRARYCPP})

#-------------------------------------------#
INSTALL(TARGETS ${LIBRARYCPPY} DESTINATION lib)
FOREACH(DIRECTORY ${DIRECTORIES})
  FILE(GLOB HFILES ${DIRECTORY}/${DIRECTORY}/*.hpp)
  INSTALL(FILES ${HFILES} DESTINATION include/${DIRECTORY})
ENDFOREACH(DIRECTORY)
