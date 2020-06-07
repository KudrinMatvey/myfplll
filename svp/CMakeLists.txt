set(libraries ${libraries} libplll-static)
set(libraries_shared ${libraries_shared} libplll-shared)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include(${PROJECT_SOURCE_DIR}/cmake/FindEigen3.cmake)
include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include)
message(STATUS "eigen inc ${EIGEN_INCLUDES}")
message(STATUS "ashk ${CMAKE_TOOLCHAIN_FILE}")
link_directories(${pLLL_BINARY_DIR}/plll)
add_executable(svp src/svp.cpp)
add_executable(svp2 src/svp2.cpp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread -o")
add_definitions(${CMAKE_REQUIRED_FLAGS})
add_executable(helloworld src/helloworld.cpp)

# RESTinio dependencies:
# 1. ASIO or Boost::ASIO (goes as headers, vcpkg knows where)
# 2. HTTP parser
find_package(unofficial-http-parser REQUIRED)
# 3. fmtlib
find_package(fmt REQUIRED)
# RESTinio itself
find_package(restinio REQUIRED)

# Make your project dependent on restinio,
# and let cmake deal with all the headers paths and linked libs.
target_link_libraries(helloworld PRIVATE ${libraries} Eigen3::Eigen restinio::restinio)

target_link_libraries(svp ${libraries})
target_link_libraries(svp2 ${libraries} Eigen3::Eigen)