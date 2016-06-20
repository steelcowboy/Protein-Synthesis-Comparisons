cmake_minimum_required(VERSION 3.0)

project(protein_synthesis)
#Add a sub-project (A library) as a dependency
add_subdirectory(lib)

#Create the client application
add_executable(protein_synthesis ${PROJECT_SOURCE_DIR}/src/protein_synthesis.cc)
#Link the library from the subproject
target_link_libraries(protein_synthesis libsynth)


