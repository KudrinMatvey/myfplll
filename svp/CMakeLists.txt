#set(libraries ${libraries} libplll-static)
#set(libraries_shared ${libraries_shared} libplll-shared)
#set (CMAKE_CXX_STANDARD 11)
#set_source_files_properties(File2.cpp -fno-rtti)
#include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include)
#link_directories(${pLLL_BINARY_DIR}/plll)

#add_executable(svp src/svp.cpp)

#target_link_libraries(svp ${libraries})

set(libraries ${libraries} libplll-static)
set(libraries_shared ${libraries_shared} libplll-shared)

include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include)
link_directories(${pLLL_BINARY_DIR}/plll)

add_executable(svp src/svp.cpp)

target_link_libraries(svp ${libraries})
