cmake_minimum_required(VERSION 3.17)
# note, suggested name is PARENT_FOLDER_samples (e.g., parent folder is myLib, this is myLib_samples
PROJECT(mat_samples LANGUAGES CXX) 

## complex case, create a dir and do it as a standard project
#add_subdirectory(...)

# note: samples are not installed by default, 
# we do not use target_sources_wrapper here 
# but the native target_sources function

## add these lines block for each sample, or if you need more complex structure, create a subfolder
# add_executable(mat_samples)
# target_sources(mat_samples PRIVATE <<...>>)
# if (WIN32)
# 	set_target_properties(mat_samples PROPERTIES FOLDER <<SAMPLE_FOLDER>>)
# endif()
# target_link_libraries(mat_samples PRIVATE <<...>>)

## and/or compile this list to have simple 1-file samples

set(SAMPLES_LIST 
	eigenPlayground.cpp
)

foreach(s ${SAMPLES_LIST})
	get_filename_component(sn "${s}" NAME_WE)
	add_executable(matSamples_${sn})
	target_sources(matSamples_${sn} PRIVATE ${s})
	if (WIN32)
		set_target_properties(matSamples_${sn} PROPERTIES FOLDER mat/samples)
	endif()	
	target_link_libraries(matSamples_${sn} PRIVATE mat)
endforeach()