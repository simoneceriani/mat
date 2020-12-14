cmake_minimum_required(VERSION 3.17)
project(mat LANGUAGES CXX)

#------------------------------------------------------------------------------------------
# Declare dependencies, use as find_package
# find_package_wrapper(...)
find_package_wrapper(Eigen3 REQUIRED CONFIG)
find_package_wrapper(OpenMP REQUIRED)

#------------------------------------------------------------------------------------------
# prepare headers and sources

# define the library
add_library(mat STATIC)

# use the target_sources_wrapper, which wrap cmake native target_sources. 
# Wrapping is necessary to cover both build and installation step
# Normally headers are PUBLIC and cpp are PRIVATE, but you might have private headers too
# If the option NO_AUTO_SOURCE_GROUP is omitted the target sources and headers are 
# automatically grouped respectively in "Sources" and "Headers" folder in the IDE
# environment whereas the dependency headers (i.e. those from linked targets) are grouped
# in the "Dependency Headers" folder.
# Specify NO_AUTO_SOURCE_GROUP if you want all sources and headers grouped in the IDE
# standard folders ("Source Files" and "Header Files" in Visual Studio).
# note: you can repeat target_sources_wrapper if you want to organize your files
target_sources_wrapper(mat
	PUBLIC 
    Global.h
		DimensionDescriptor.h
		DimensionDescriptor.hpp
    DimensionDescriptorTraits.h
    DimensionDescriptorTraits.hpp
    MatrixBlockDescriptor.h
    MatrixBlockDescriptor.hpp
    MatrixBlockTypeTraits.h
    DenseMatrixBlock.h
    DenseMatrixBlock.hpp
    DiagonalMatrixBlock.h
    DiagonalMatrixBlock.hpp
    MatrixBlockBase.h
    MatrixBlockBase.hpp
    SparsityPattern.h
    SparseMatrixBlock.h
    SparseMatrixBlock.hpp
	PRIVATE
    Global.cpp
		DimensionDescriptor.cpp
    DimensionDescriptorTraits.cpp
    MatrixBlockDescriptor.cpp
    DenseMatrixBlock.cpp
    DiagonalMatrixBlock.cpp
    MatrixBlockBase.cpp
    SparsityPattern.cpp
    SparseMatrixBlock.cpp
)

# if you do not want to have matd.{dll,lib,so} in debug disable this 
set_target_properties(mat  PROPERTIES DEBUG_POSTFIX "d")

# special compiler flags for this target
# PUBLIC will propagate do dependent, PRIVATE if they are only for this project
#if(MSVC)
#	target_compile_options(mat <<PUBLIC|PRIVATE>> /bigobj)
#endif()
#target_compile_features(mat PUBLIC cxx_std_11) # cxx_std_14

#Add an alias so that library can be used inside the build tree
add_library(mat::mat ALIAS mat)

# link libraries
target_link_libraries(mat
	PUBLIC
		Eigen3::Eigen
    OpenMP::OpenMP_CXX
	PRIVATE
)

# add qt automoc per target if you need (or remove this section if qt independent)
# AUTOMOC for signal slot, AUTOUIC for User Interfaces, AUTORCC for resource files
# set_target_properties(mat PROPERTIES AUTOMOC ON)
# set_target_properties(mat PROPERTIES AUTOUIC ON)
# set_target_properties(mat PROPERTIES AUTORCC ON)

#------------------------------------------------------------------------------------------
# GLOBAL, SELF or a custom export name allowed
installAllTargets(SELF)

#------------------------------------------------------------------------------------------
# this will add the samples and tests folder (if existing) and will perfom common operations:
# - link the target to the matVersion library, which allow to have the root include directory added
# - create folder structure for MSVC IDE
# - analyze Dll dependencies per target and globally, which in turn will be used for deep installation and PATH generation
commonEndProjectCMake(matSol)
