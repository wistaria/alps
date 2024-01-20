######################################################################
# Symmetries management
######################################################################
set(DMRG_BUILD_SYMMETRIES "NU1" CACHE STRING "List of symmetry to include in the build objects. U1, TwoU1, NU1, Z2 and NONE are currently available.")
set(DMRG_NUMSYMM 6 CACHE STRING "Maximum number of U1 symmetries for NU1.")
mark_as_advanced(DMRG_BUILD_SYMMETRIES DMRG_NUMSYMM)

macro(get_symm_suffix RET SYMM)
  set(${RET} ${SYMM})
  string(REGEX REPLACE "^U1$" "u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TWOU1$" "2u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TwoU1$" "2u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^NU1$"  "nu1" ${RET} ${${RET}})
  string(REGEX REPLACE "^NONE$" "none" ${RET} ${${RET}})
  string(REGEX REPLACE "^Z2$" "Ztwo" ${RET} ${${RET}})
endmacro(get_symm_suffix)

macro(get_symm_group_name RET SYMM)
  set(${RET} ${SYMM})
  string(REGEX REPLACE "^U1$" "U1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TWOU1$" "TwoU1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TwoU1$" "TwoU1" ${RET} ${${RET}})
  string(REGEX REPLACE "^NU1$"  "NU1" ${RET} ${${RET}})
  string(REGEX REPLACE "^NONE$" "TrivialGroup" ${RET} ${${RET}})
  string(REGEX REPLACE "^Z2$" "Ztwo" ${RET} ${${RET}})
endmacro(get_symm_group_name)

macro(get_symm_files TYPE RET FILEBASE)
  if(NOT ${TYPE} STREQUAL "APPEND")
    set(${RET} "")
  endif(NOT ${TYPE} STREQUAL "APPEND")
  foreach(SYMM ${DMRG_BUILD_SYMMETRIES})
    get_symm_suffix(SYMM_SUFFIX ${SYMM})
    string(REPLACE "{SYMM}" ${SYMM_SUFFIX} SYMM_FILE ${FILEBASE})
    list(APPEND ${RET} ${SYMM_FILE})
  endforeach(SYMM)
endmacro(get_symm_files)

macro(configure_symm_file INPUT OUTBASE VARNAME)
  foreach(SYMM ${DMRG_BUILD_SYMMETRIES})
    get_symm_suffix(SYMM_SUFFIX ${SYMM})
    get_symm_group_name(${VARNAME} ${SYMM})
    string(REPLACE "{SYMM}" ${SYMM_SUFFIX} SYMM_FILE ${OUTBASE})
    configure_file(${INPUT} ${SYMM_FILE})
  endforeach(SYMM)
endmacro(configure_symm_file)


macro(set_symmetry_difinitions RET)
  foreach(SYMM ${DMRG_BUILD_SYMMETRIES})
    get_symm_group_name(SYMM_NAME ${SYMM})
    message(STATUS "MPS: enabling ${SYMM_NAME} symmetry.")
    list(APPEND ${RET} -DHAVE_${SYMM_NAME})
  endforeach(SYMM)
endmacro(set_symmetry_difinitions)
