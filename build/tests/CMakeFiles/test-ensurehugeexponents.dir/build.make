# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /root/myfplll

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /root/myfplll/build

# Include any dependencies generated for this target.
include tests/CMakeFiles/test-ensurehugeexponents.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/test-ensurehugeexponents.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/test-ensurehugeexponents.dir/flags.make

tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o: tests/CMakeFiles/test-ensurehugeexponents.dir/flags.make
tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o: ../tests/src/test-ensurehugeexponents.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/myfplll/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o"
	cd /root/myfplll/build/tests && /usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o -c /root/myfplll/tests/src/test-ensurehugeexponents.cpp

tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.i"
	cd /root/myfplll/build/tests && /usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/myfplll/tests/src/test-ensurehugeexponents.cpp > CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.i

tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.s"
	cd /root/myfplll/build/tests && /usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/myfplll/tests/src/test-ensurehugeexponents.cpp -o CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.s

tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o.requires:

.PHONY : tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o.requires

tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o.provides: tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o.requires
	$(MAKE) -f tests/CMakeFiles/test-ensurehugeexponents.dir/build.make tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o.provides.build
.PHONY : tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o.provides

tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o.provides.build: tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o


# Object files for target test-ensurehugeexponents
test__ensurehugeexponents_OBJECTS = \
"CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o"

# External object files for target test-ensurehugeexponents
test__ensurehugeexponents_EXTERNAL_OBJECTS =

tests/test-ensurehugeexponents: tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o
tests/test-ensurehugeexponents: tests/CMakeFiles/test-ensurehugeexponents.dir/build.make
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_thread.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_system.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_chrono.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_date_time.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_atomic.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libpthread.so
tests/test-ensurehugeexponents: /usr/local/lib/libgmp.so
tests/test-ensurehugeexponents: /usr/local/lib/libmpfr.so
tests/test-ensurehugeexponents: plll/libplll.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_thread.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_system.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_chrono.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_date_time.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libboost_atomic.a
tests/test-ensurehugeexponents: /usr/lib/x86_64-linux-gnu/libpthread.so
tests/test-ensurehugeexponents: /usr/local/lib/libgmp.so
tests/test-ensurehugeexponents: /usr/local/lib/libmpfr.so
tests/test-ensurehugeexponents: tests/CMakeFiles/test-ensurehugeexponents.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/root/myfplll/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test-ensurehugeexponents"
	cd /root/myfplll/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test-ensurehugeexponents.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/test-ensurehugeexponents.dir/build: tests/test-ensurehugeexponents

.PHONY : tests/CMakeFiles/test-ensurehugeexponents.dir/build

tests/CMakeFiles/test-ensurehugeexponents.dir/requires: tests/CMakeFiles/test-ensurehugeexponents.dir/src/test-ensurehugeexponents.cpp.o.requires

.PHONY : tests/CMakeFiles/test-ensurehugeexponents.dir/requires

tests/CMakeFiles/test-ensurehugeexponents.dir/clean:
	cd /root/myfplll/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/test-ensurehugeexponents.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/test-ensurehugeexponents.dir/clean

tests/CMakeFiles/test-ensurehugeexponents.dir/depend:
	cd /root/myfplll/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /root/myfplll /root/myfplll/tests /root/myfplll/build /root/myfplll/build/tests /root/myfplll/build/tests/CMakeFiles/test-ensurehugeexponents.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/test-ensurehugeexponents.dir/depend

