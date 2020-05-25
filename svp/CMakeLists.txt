set(libraries ${libraries} libplll-static)
set(libraries_shared ${libraries_shared} libplll-shared)

include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include ${PROJECT_BINARY_DIR}/packages/eigen)
message(STATUS "eigen inc ${EIGEN_INCLUDES}")
link_directories(${pLLL_BINARY_DIR}/plll)

add_executable(svp src/svp.cpp)
add_executable(svp2 src/svp2.cpp)

target_link_libraries(svp ${libraries})
target_link_libraries(svp2 ${libraries})
