cmake_minimum_required(VERSION 3.17)
PROJECT(matTest  LANGUAGES CXX)

include(CTest)
find_package(Catch2 REQUIRED)
include(Catch)

# note: tests are not installed by default, 

# add these lines block for each test, or if you need more complex structure, create a subfolder
#add_executable(<<TARGET_NAME>>)
#target_sources(<<TARGET_NAME>> PRIVATE <<...>>)
#if (WIN32)
#	set_target_properties(<<TARGET_NAME>> PROPERTIES FOLDER mat/test)
#endif()
#target_link_libraries(<<TARGET_NAME>> PRIVATE <<...>> Boost::unit_test_framework)
#add_test(NAME <<TARGET_NAME>> COMMAND <<TARGET_NAME>>)

# and/or compile this list to have simple 1-file test

set(TESTS_LIST 
	DimensionDescriptorCatch.cpp
  MatrixBlockDescriptorCatch.cpp
  DenseMatrixBlockCatch.cpp
  DiagonalMatrixBlockCatch.cpp
  SparsityPatternCatch.cpp
  SparseMatrixBlockCatch.cpp
  SparsityPatternBlockDescriptorCatch.cpp
)

foreach(s ${TESTS_LIST})
	get_filename_component(sn "${s}" NAME_WE)
	add_executable(matTest_${sn})
	target_sources(matTest_${sn} PRIVATE ${s})
	if (WIN32)
		set_target_properties(matTest_${sn} PROPERTIES FOLDER mat/test)
	endif()	
	target_link_libraries(matTest_${sn} PRIVATE mat::mat Catch2::Catch2WithMain )
  
  # add to tests
  catch_discover_tests(matTest_${sn})
endforeach()
