
# ==============================================================================
#                              Save variables to file
# ==============================================================================

MESSAGE(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
set(CacheForExecutable ${CMAKE_INSTALL_PREFIX}/cmake/CMakeCacheForExecutable.cmake)
file(WRITE ${CacheForExecutable} "")
set(Vars
SIMFEM_EXTERNAL_LINKLIBS
SIMFEM_EXTERNAL_INCLUDES
SIMFEM_INCLUDES
)
foreach(Var ${Vars})
  file(APPEND ${CacheForExecutable} "set(${Var} \"${${Var}}\")\n")
endforeach()
