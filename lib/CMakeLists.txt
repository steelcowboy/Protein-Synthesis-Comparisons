cmake_minimum_required(VERSION 3.0)

project(libsynth)

#Create a library output (Instead of an executable)
add_library(libsynth ${PROJECT_SOURCE_DIR}/src/synthesis.cc)

#Specify the include directories of this project, and allow other superprojects to include these directories as well.
target_include_directories(
   libsynth PUBLIC
   src)
