set(libraries ${libraries} libplll-static)
set(libraries_shared ${libraries_shared} libplll-shared)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include(${PROJECT_SOURCE_DIR}/cmake/FindEigen3.cmake)
include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include)
message(STATUS "eigen inc ${EIGEN_INCLUDES}")
link_directories(${pLLL_BINARY_DIR}/plll)
add_executable(svp src/svp.cpp)
add_executable(svp2 src/svp2.cpp)

target_link_libraries(svp ${libraries})
target_link_libraries(svp2 ${libraries} Eigen3::Eigen)
