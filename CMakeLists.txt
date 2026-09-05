# ELiPS -- pairing-based cryptography over BN and BLS12 curves.
#
# Replaces the autotools build (issue #2). See MODERNIZATION_PLAN.md.
#
#   cmake -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j
#   ctest --test-dir build --output-on-failure
#
# Sanitizer builds:
#   cmake -B build-asan -DCMAKE_BUILD_TYPE=Asan
#
# Note on curve selection: the plan (section 1.5) calls for compile-time curve
# selection via -DELIPS_CURVE, so that FP_LIMBS becomes a constant. That only
# becomes meaningful in Phase 3, when field elements stop being heap mpz_t.
# Today the library selects its curve at runtime through init_bn() /
# bls12_inits(), so adding the switch now would be scaffolding with nothing
# behind it. Deferred to Phase 3 deliberately.

cmake_minimum_required(VERSION 3.16)
project(elips
  VERSION 1.0.1
  DESCRIPTION "Efficient Library for Pairing-based Systems (BN, BLS12)"
  LANGUAGES C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(ELIPS_BUILD_TESTS "Build the test suite" ON)
option(ELIPS_WERROR      "Treat warnings as errors" OFF)

# ---------------------------------------------------------------- sanitizers
# Exposed as build types so a sanitizer run is one cache variable, not a
# hand-assembled flag string.
set(_san_common "-O1 -g -fno-omit-frame-pointer")
set(CMAKE_C_FLAGS_ASAN  "${_san_common} -fsanitize=address,undefined -fno-sanitize-recover=all"
    CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_UBSAN "${_san_common} -fsanitize=undefined -fno-sanitize-recover=all"
    CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_MSAN  "${_san_common} -fsanitize=memory -fsanitize-memory-track-origins"
    CACHE STRING "" FORCE)
foreach(_t ASAN UBSAN MSAN)
  set(CMAKE_EXE_LINKER_FLAGS_${_t}    "${CMAKE_C_FLAGS_${_t}}" CACHE STRING "" FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_${_t} "${CMAKE_C_FLAGS_${_t}}" CACHE STRING "" FORCE)
endforeach()
mark_as_advanced(CMAKE_C_FLAGS_ASAN CMAKE_C_FLAGS_UBSAN CMAKE_C_FLAGS_MSAN)

find_package(GMP REQUIRED)

# ------------------------------------------------------------------- library
file(GLOB ELIPS_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")
file(GLOB ELIPS_HEADERS CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/include/ELiPS_bn_bls/*.h")

add_library(elips ${ELIPS_SOURCES})
add_library(ELiPS::elips ALIAS elips)

target_include_directories(elips
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(elips PUBLIC GMP::GMP)
set_target_properties(elips PROPERTIES
  VERSION   ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  C_VISIBILITY_PRESET default
  POSITION_INDEPENDENT_CODE ON)

target_compile_options(elips PRIVATE -Wall -Wextra)
if(ELIPS_WERROR)
  target_compile_options(elips PRIVATE -Werror)
endif()

# --------------------------------------------------------------------- tests
if(ELIPS_BUILD_TESTS)
  enable_testing()

  # Known-answer vectors. One process per curve: the library keeps its curve
  # parameters in globals, so init_bn() and bls12_inits() overwrite each other.
  # That constraint goes away in Phase 3.
  add_executable(kat_runner test/kat_runner.c)
  target_link_libraries(kat_runner PRIVATE elips)

  foreach(_curve bn_462 bls12_461 bls12_381)
    add_test(NAME kat.${_curve}
             COMMAND kat_runner "${CMAKE_CURRENT_SOURCE_DIR}/test/kat/${_curve}.vec")
  endforeach()

  # Negative control: the runner must actually fail on a corrupted vector.
  # A test oracle that cannot fail is not an oracle, so this is not optional.
  add_test(NAME kat.detects_corruption
           COMMAND ${CMAKE_COMMAND}
                   -DRUNNER=$<TARGET_FILE:kat_runner>
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/test/kat/bn_462.vec
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorruptCheck.cmake)

  add_executable(property_tests test/property_tests.c)
  target_link_libraries(property_tests PRIVATE elips)
  foreach(_c bn bls12)
    add_test(NAME property.${_c} COMMAND property_tests ${_c})
  endforeach()

  # The library's two final exponentiations must agree and currently do not.
  # Registered as an expected failure so the bug is visible in CI without
  # blocking the gate. Remove WILL_FAIL when issue #16 is fixed.
  add_executable(finalexp_agreement test/finalexp_agreement.c)
  target_link_libraries(finalexp_agreement PRIVATE elips)
  foreach(_c bn bls12)
    add_test(NAME finalexp.agreement.${_c} COMMAND finalexp_agreement ${_c})
    set_tests_properties(finalexp.agreement.${_c} PROPERTIES WILL_FAIL TRUE)
  endforeach()

  set_tests_properties(kat.bn_462 kat.bls12_461 kat.bls12_381
                       property.bn property.bls12
                       PROPERTIES TIMEOUT 600)
endif()

# ------------------------------------------------------------------- install
include(GNUInstallDirs)
install(TARGETS elips EXPORT elipsTargets
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES ${ELIPS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ELiPS_bn_bls)
install(EXPORT elipsTargets
        FILE ELiPSTargets.cmake NAMESPACE ELiPS::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ELiPS)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ELiPSConfigVersion.cmake"
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ELiPSConfig.cmake"
"include(CMakeFindDependencyMacro)
list(APPEND CMAKE_MODULE_PATH \"\${CMAKE_CURRENT_LIST_DIR}\")
find_dependency(GMP)
include(\"\${CMAKE_CURRENT_LIST_DIR}/ELiPSTargets.cmake\")
")
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/ELiPSConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/ELiPSConfigVersion.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGMP.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ELiPS)

message(STATUS "ELiPS ${PROJECT_VERSION}  build=${CMAKE_BUILD_TYPE}  GMP=${GMP_VERSION}")
