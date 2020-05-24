# Install script for directory: /home/user/plll-1.0/plll

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/user/plll-1.0/build/plll/libplll.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libplll.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libplll.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libplll.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/user/plll-1.0/build/plll/libplll.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libplll.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libplll.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libplll.so"
         OLD_RPATH "/usr/local/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libplll.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/home/user/plll-1.0/plll/include/plll.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/plll" TYPE FILE FILES
    "/home/user/plll-1.0/build/plll/include/plll/config.hpp"
    "/home/user/plll-1.0/plll/include/plll/arguments.hpp"
    "/home/user/plll-1.0/plll/include/plll/arithmetic.hpp"
    "/home/user/plll-1.0/plll/include/plll/arithmetic-expressions.hpp"
    "/home/user/plll-1.0/plll/include/plll/arithmetic-gmp.hpp"
    "/home/user/plll-1.0/plll/include/plll/arithmetic-gmp-conv.hpp"
    "/home/user/plll-1.0/plll/include/plll/arithmetic-gmp-iops.hpp"
    "/home/user/plll-1.0/plll/include/plll/arithmetic-gmp-rops.hpp"
    "/home/user/plll-1.0/plll/include/plll/arithmetic-nint.hpp"
    "/home/user/plll-1.0/plll/include/plll/arithmetic-nint-conv.hpp"
    "/home/user/plll-1.0/plll/include/plll/documentation.hpp"
    "/home/user/plll-1.0/plll/include/plll/helper.hpp"
    "/home/user/plll-1.0/plll/include/plll/linalg.hpp"
    "/home/user/plll-1.0/plll/include/plll/matrix.hpp"
    "/home/user/plll-1.0/plll/include/plll/matrix-mem.hpp"
    "/home/user/plll-1.0/plll/include/plll/matrix-ops.hpp"
    "/home/user/plll-1.0/plll/include/plll/matrix-ops2.hpp"
    "/home/user/plll-1.0/plll/include/plll/rational.hpp"
    "/home/user/plll-1.0/plll/include/plll/rational-conv.hpp"
    "/home/user/plll-1.0/plll/include/plll/rational-ops.hpp"
    )
endif()

