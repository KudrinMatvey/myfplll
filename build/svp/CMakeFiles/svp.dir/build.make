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
include svp/CMakeFiles/svp.dir/depend.make

# Include the progress variables for this target.
include svp/CMakeFiles/svp.dir/progress.make

# Include the compile flags for this target's objects.
include svp/CMakeFiles/svp.dir/flags.make

svp/CMakeFiles/svp.dir/src/svp.cpp.o: svp/CMakeFiles/svp.dir/flags.make
svp/CMakeFiles/svp.dir/src/svp.cpp.o: ../svp/src/svp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/myfplll/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object svp/CMakeFiles/svp.dir/src/svp.cpp.o"
	cd /root/myfplll/build/svp && /usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/svp.dir/src/svp.cpp.o -c /root/myfplll/svp/src/svp.cpp

svp/CMakeFiles/svp.dir/src/svp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/svp.dir/src/svp.cpp.i"
	cd /root/myfplll/build/svp && /usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/myfplll/svp/src/svp.cpp > CMakeFiles/svp.dir/src/svp.cpp.i

svp/CMakeFiles/svp.dir/src/svp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/svp.dir/src/svp.cpp.s"
	cd /root/myfplll/build/svp && /usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/myfplll/svp/src/svp.cpp -o CMakeFiles/svp.dir/src/svp.cpp.s

svp/CMakeFiles/svp.dir/src/svp.cpp.o.requires:

.PHONY : svp/CMakeFiles/svp.dir/src/svp.cpp.o.requires

svp/CMakeFiles/svp.dir/src/svp.cpp.o.provides: svp/CMakeFiles/svp.dir/src/svp.cpp.o.requires
	$(MAKE) -f svp/CMakeFiles/svp.dir/build.make svp/CMakeFiles/svp.dir/src/svp.cpp.o.provides.build
.PHONY : svp/CMakeFiles/svp.dir/src/svp.cpp.o.provides

svp/CMakeFiles/svp.dir/src/svp.cpp.o.provides.build: svp/CMakeFiles/svp.dir/src/svp.cpp.o


# Object files for target svp
svp_OBJECTS = \
"CMakeFiles/svp.dir/src/svp.cpp.o"

# External object files for target svp
svp_EXTERNAL_OBJECTS =

svp/svp: svp/CMakeFiles/svp.dir/src/svp.cpp.o
svp/svp: svp/CMakeFiles/svp.dir/build.make
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_thread.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_system.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_chrono.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_date_time.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_atomic.a
svp/svp: /usr/lib/x86_64-linux-gnu/libpthread.so
svp/svp: /usr/local/lib/libgmp.so
svp/svp: /usr/local/lib/libmpfr.so
svp/svp: plll/libplll.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_thread.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_system.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_chrono.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_date_time.a
svp/svp: /usr/lib/x86_64-linux-gnu/libboost_atomic.a
svp/svp: /usr/lib/x86_64-linux-gnu/libpthread.so
svp/svp: /usr/local/lib/libgmp.so
svp/svp: /usr/local/lib/libmpfr.so
svp/svp: svp/CMakeFiles/svp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/root/myfplll/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable svp"
	cd /root/myfplll/build/svp && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/svp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
svp/CMakeFiles/svp.dir/build: svp/svp

.PHONY : svp/CMakeFiles/svp.dir/build

svp/CMakeFiles/svp.dir/requires: svp/CMakeFiles/svp.dir/src/svp.cpp.o.requires

.PHONY : svp/CMakeFiles/svp.dir/requires

svp/CMakeFiles/svp.dir/clean:
	cd /root/myfplll/build/svp && $(CMAKE_COMMAND) -P CMakeFiles/svp.dir/cmake_clean.cmake
.PHONY : svp/CMakeFiles/svp.dir/clean

svp/CMakeFiles/svp.dir/depend:
	cd /root/myfplll/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /root/myfplll /root/myfplll/svp /root/myfplll/build /root/myfplll/build/svp /root/myfplll/build/svp/CMakeFiles/svp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : svp/CMakeFiles/svp.dir/depend

