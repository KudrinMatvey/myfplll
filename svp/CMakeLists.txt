set(libraries ${libraries} libplll-static)
set(libraries_shared ${libraries_shared} libplll-shared)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include(${PROJECT_SOURCE_DIR}/cmake/FindEigen3.cmake)
include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include)
message(STATUS "eigen inc ${EIGEN_INCLUDES}")
message(STATUS "ashk ${CMAKE_TOOLCHAIN_FILE}")
link_directories(${pLLL_BINARY_DIR}/plll)
link_directories(${PROJECT_SOURCE_DIR}/svp/src/denis)
message(STATUS "ashk ${PROJECT_SOURCE_DIR}/svp/src/denis")
message(STATUS "${PLLL_FOUND} ${EIGEN_FOUND}")

add_executable(svp src/svp.cpp)
add_executable(d-svp src/denis/svp.cpp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread -o")
add_definitions(${CMAKE_REQUIRED_FLAGS})
add_executable(helloworld src/helloworld.cpp)

file(GLOB SOURCES
    ${PROJECT_SOURCE_DIR}/svp/src/denis/*.h
    ${PROJECT_SOURCE_DIR}/svp/src/denis/svp.cpp
)
add_executable(d-main src/denis/main.cpp ${SOURCES})
add_executable(svp2 src/svp2.cpp ${SOURCES})
add_executable(main2 src/kudrin/main2.cpp src/kudrin/ksvp.cpp src/kudrin/ksvp.h)


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
target_link_libraries(svp ${libraries})
target_link_libraries(d-svp ${libraries})
target_link_libraries(d-main ${libraries})
target_link_libraries(svp2 ${libraries} Eigen3::Eigen)
target_link_libraries(main2 ${libraries} Eigen3::Eigen)
