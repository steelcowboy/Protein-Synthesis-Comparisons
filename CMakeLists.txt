cmake_minimum_required(VERSION 3.0)

project(protein_synthesis)
#Add a sub-project (A library) as a dependency
add_subdirectory(lib)

#Create the client application
add_executable(protein_synthesis ${PROJECT_SOURCE_DIR}/src/protein_synthesis.cc)

#Give more warnings
if ( CMAKE_COMPILER_IS_GNUCC )
    set_property( TARGET protein_synthesis APPEND_STRING PROPERTY COMPILE_FLAGS -Werror)
endif ( CMAKE_COMPILER_IS_GNUCC )

#Link the library from the subproject
target_link_libraries(protein_synthesis libsynth)


