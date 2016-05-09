
# a2display is needed in Core,
# because AcquDAQ can't be properly disentangled
# from AcquRoot analysis part grrrr
add_subdirectory(a2display)
if(ACQU_DAQ)
  return()
endif()

add_subdirectory(Worker)
# check for external OSCAR
if(DEFINED ENV{OSCAR})
  message("-- Found external OSCAR in " $ENV{OSCAR})
  include_directories($ENV{OSCAR}/include)
  link_directories($ENV{OSCAR}/lib)
  set(ACQU_USER_LIBS ${ACQU_USER_LIBS} OSCAR)
else()
  add_subdirectory(OSCAR)
endif()

# bump the variables to parent scope set by the macro SET_ACQU_USER_VARS
set(ACQU_USER_INCS ${ACQU_USER_INCS} PARENT_SCOPE)
set(ACQU_USER_LIBS ${ACQU_USER_LIBS} PARENT_SCOPE)
