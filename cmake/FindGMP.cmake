# Locate GMP. Honours GMP_ROOT, then pkg-config, then the usual prefixes
# (including Homebrew on both Apple Silicon and Intel).
#
# Defines the imported target GMP::GMP.

find_path(GMP_INCLUDE_DIR gmp.h
  HINTS ${GMP_ROOT} ENV GMP_ROOT
  PATH_SUFFIXES include
  PATHS /opt/homebrew /usr/local /usr /opt/local)

find_library(GMP_LIBRARY NAMES gmp libgmp
  HINTS ${GMP_ROOT} ENV GMP_ROOT
  PATH_SUFFIXES lib
  PATHS /opt/homebrew /usr/local /usr /opt/local)

if(GMP_INCLUDE_DIR AND EXISTS "${GMP_INCLUDE_DIR}/gmp.h")
  file(STRINGS "${GMP_INCLUDE_DIR}/gmp.h" _gmp_ver_major
       REGEX "^#define[ \t]+__GNU_MP_VERSION[ \t]+[0-9]+")
  file(STRINGS "${GMP_INCLUDE_DIR}/gmp.h" _gmp_ver_minor
       REGEX "^#define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+[0-9]+")
  file(STRINGS "${GMP_INCLUDE_DIR}/gmp.h" _gmp_ver_patch
       REGEX "^#define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+[0-9]+")
  string(REGEX MATCH "[0-9]+" _maj "${_gmp_ver_major}")
  string(REGEX MATCH "[0-9]+" _min "${_gmp_ver_minor}")
  string(REGEX MATCH "[0-9]+" _pat "${_gmp_ver_patch}")
  set(GMP_VERSION "${_maj}.${_min}.${_pat}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  REQUIRED_VARS GMP_LIBRARY GMP_INCLUDE_DIR
  VERSION_VAR GMP_VERSION)

if(GMP_FOUND AND NOT TARGET GMP::GMP)
  add_library(GMP::GMP UNKNOWN IMPORTED)
  set_target_properties(GMP::GMP PROPERTIES
    IMPORTED_LOCATION "${GMP_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
