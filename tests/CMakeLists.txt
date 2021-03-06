set(libraries ${libraries} libplll-static)
set(libraries_shared ${libraries_shared} libplll-shared)

include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include ${pLLL_SOURCE_DIR}/plll/include-internal)
link_directories(${pLLL_BINARY_DIR}/plll)

add_executable(test-arithmetic src/test-arithmetic.cpp)
add_executable(test-int src/test-int.cpp)
add_executable(test-nfpint src/test-nfpint.cpp)
add_executable(test-rng src/test-rng.cpp)
add_executable(test-vecmat src/test-vecmat.cpp)
add_executable(test-real src/test-real.cpp)
add_executable(test-hnf src/test-hnf.cpp)
add_executable(test-voronoi src/test-voronoi.cpp)
add_executable(test-plll src/test-plll.cpp)
add_executable(test-taskman src/test-taskman.cpp)
add_executable(test-ensurehugeexponents src/test-ensurehugeexponents.cpp)
add_executable(timing-vecmat src/timing-vecmat.cpp)

target_link_libraries(test-arithmetic ${libraries})
target_link_libraries(test-int ${libraries})
target_link_libraries(test-nfpint ${libraries})
target_link_libraries(test-rng ${libraries})
target_link_libraries(test-vecmat ${libraries})
target_link_libraries(test-real ${libraries})
target_link_libraries(test-hnf ${libraries})
target_link_libraries(test-voronoi ${libraries})
target_link_libraries(test-plll ${libraries})
target_link_libraries(test-taskman ${libraries})
target_link_libraries(test-ensurehugeexponents ${libraries})
target_link_libraries(timing-vecmat ${libraries})
