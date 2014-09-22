cmake_minimum_required (VERSION 2.6)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/subdir-header.cmake)
message(STATUS "Configuring CaLib...")

find_package(ZLIB REQUIRED)
find_package(Curses REQUIRED)

include_directories("include")



set(MOD "CaLib")
ROOT_PREPARE_LIB("src" "include" "include" "src/G__${MOD}.cxx" SRCS EXESRCS)
add_library("lib${MOD}" SHARED ${SRCS})
target_link_libraries("lib${MOD}" ${ROOT_LIBRARIES} ${CURSES_LIBRARIES} ${ZLIB_LIBRARIES} Spectrum)
set_target_properties("lib${MOD}" PROPERTIES OUTPUT_NAME ${MOD})

# we add many additional files marked as "headers"
# so they appear automatically in QtCreator
file(GLOB_RECURSE ADDITIONAL_FILES "macros/*" "config/*")
set_source_files_properties(${ADDITIONAL_FILES} PROPERTIES HEADER_FILE_ONLY TRUE)

# and we create the calib_manager
add_executable("calib_manager" "src/Main${MOD}Manager.cxx" ${ADDITIONAL_FILES})
target_link_libraries("calib_manager" "lib${MOD}")

# now track the incs/libs for acqu_user
SET_ACQU_USER_VARS("lib${MOD}")
