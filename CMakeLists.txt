CMAKE_MINIMUM_REQUIRED(VERSION 2.8) 

PROJECT (b00203579_fft)

#
#	dir variables
#
SET(incl_dir "${CMAKE_CURRENT_SOURCE_DIR}/incl")
MESSAGE (STATUS  "** incl_dir - ${incl_dir}")
SET(src_dir "${CMAKE_CURRENT_SOURCE_DIR}/src")
MESSAGE (STATUS  "** src_dir - ${src_dir}")

SET(fftw_dir "${CMAKE_CURRENT_SOURCE_DIR}/fftw3")
MESSAGE (STATUS  "** fftw_dir - ${fftw_dir}")

SET(res_dir "${CMAKE_CURRENT_BINARY_DIR}/res")
MESSAGE (STATUS  "** res_dir - ${res_dir}")

#
#	source files
#
FILE(GLOB_RECURSE PROJECT_INCL_FILES "${incl_dir}/*.h")
FILE(GLOB PROJECT_SRC_FILES "${src_dir}/*.cpp")

#
#	register project as executable
#
ADD_EXECUTABLE( ${CMAKE_PROJECT_NAME} 
				${PROJECT_INCL_FILES} ${PROJECT_SRC_FILES} )
TARGET_LINK_LIBRARIES(	${CMAKE_PROJECT_NAME} "${fftw_dir}/libfftw3-3.lib")
SET_TARGET_PROPERTIES(	${CMAKE_PROJECT_NAME} 
						PROPERTIES INCLUDE_DIRECTORIES "${incl_dir};${fftw_dir}")

#
#	output directory and path to executable
#
SET(exec_dir $<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>) 

#
#	copy libfftw3-3.dll into the same directory as the built executable
#
ADD_CUSTOM_COMMAND(	TARGET ${CMAKE_PROJECT_NAME} POST_BUILD
               		COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   	"${fftw_dir}/libfftw3-3.dll" ${exec_dir})
