cmake_minimum_required(VERSION 3.17)

# tell cmake where is CMakeBaseBuildSystem folder
set(matSol_CMAKEBASEBUILDSYSTEM_PATH "${CMAKE_CURRENT_LIST_DIR}/../../cmakebasebuildsystem" CACHE PATH "path to the CMakeBaseBuildSystem folder")
list(APPEND CMAKE_MODULE_PATH "${matSol_CMAKEBASEBUILDSYSTEM_PATH}")

# OK, now a normal project can start
project(matSol VERSION 1.0.0 LANGUAGES CXX)

# add the CMakeBaseBuildSystem modules, just a bunch of includes here
include(winDllRecursiveAnalysis)
include(utils)
include(CMakeBaseBuildSystem)
include(gitTracer)
# add CMakeBaseBuildSystem in source dir, allow to track version!
configureGitCMakeFile(${matSol_CMAKEBASEBUILDSYSTEM_PATH} CMakeBaseBuildSystem ${CMAKE_CURRENT_SOURCE_DIR})

if (WIN32)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

	# global definitions, edit as you need, but try as much as you can to use specific target properties!
	
	# define a single output folder for binaries (dll and exe)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
	add_definitions(-DUNICODE)
	
	# silent generally an annoying warning from visual studio about dll export/import
	add_compile_options(/wd4251)
	
	# add multicore compilation in visual studio
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
	
endif()

# message on cmake_prefix_path, it is always nice to know where I am looking too!
MESSAGE(STATUS "Packages search paths:")
MESSAGE(STATUS "ENV{CMAKE_PREFIX_PATH} = $ENV{CMAKE_PREFIX_PATH}")
MESSAGE(STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")

# enable cmake testing
enable_testing()

# add myself version, but in binary dir, otherwise each commit will dirty itself, possible future use (so far commented)
# configureGitCMakeFile(${CMAKE_CURRENT_SOURCE_DIR} matSol ${CMAKE_CURRENT_BINARY_DIR})

# add fixed for libraries without modern config files if you need, 
# divide them in common and unix/win specific 
# NOTE: content of this folder will be copied at installation step
# if you need some special cmake additional file please put it somewhere else
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
#if(UNIX)
#	list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/unix")
#endif()
#if(WIN32)
#	list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/win32")
#endif()
#installSpecialCMake("${CMAKE_CURRENT_LIST_DIR}/cmake/")


# set default option values
# set ON if your project is an application and need to be standalone, OFF if is a library to be distributed
OPTION(matSol_INSTALL_DEPENDENCIES "Install Dll dependencies" OFF)
OPTION(matSol_BUILD_SAMPLES "Global switch for samples building" ON)
OPTION(matSol_BUILD_TESTS "Global switch for tests building" ON)

# installation destination
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/installation" CACHE PATH "prefix to installation folder" FORCE)
endif()

# prepare vesion.h and version.cpp and link them in a library with name 
# NOTE: this will prepare and install a target <<MAIN_PROJECT_NAME>>Version (STATIC LIBRARY) 
prepareVersionFile(matSol)
# alternative, use prepareVersionFileInterface(<<MAIN_PROJECT_NAME>>) to create <<MAIN_PROJECT_NAME>>Version as INTERFACE LIBRARY

# add the libs (please respect the structure)
add_subdirectory(mat)
# add_subdirectory(...)

#--------------------------------------------------------------------------------------
# create the config/configVersion file, 
# first argument, name of the main project 
# use additional arguments to specify defaul list of components when this will be used in FIND_PACKAGE. 
# NONE ALL or custom list are accepted
# note: if I am not the main project skip (submodule)
if(matSol STREQUAL ${CMAKE_PROJECT_NAME})
	createGlobalConfigFile(mat ALL)
endif()

# install other cross dependent dll/so if requested
if(${matSol_INSTALL_DEPENDENCIES})
	prepareDLLInstall(DEBUG RELEASE)
endif()

# configure the bat file for VS, cmd with PATH, this will be produced in build folder
if(MSVC)
	prepareBatFile()
endif()

# only a print so far, TODO properly plugin management for installation!
# include(QtCMakeUtils)
# printQtPluginsList()
