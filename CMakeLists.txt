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
# Curve selection: the Phase 3 arithmetic in src/arith is compiled per curve,
# because FP_LIMBS has to be a compile-time constant for the loops to unroll and
# for the Phase 6 assembly to target a known width (plan section 1.5). The
# legacy mpz layer in src/*.c still picks its curve at runtime, so the two
# coexist until the migration finishes.

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

# src/util/sysrand.c is the system CSPRNG. It carries no curve dependency, so
# both layers compile it: the legacy library needs it to seed its gmp_randstate,
# which used to come from time(NULL).
add_library(elips ${ELIPS_SOURCES} src/util/sysrand.c)
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

# ------------------------------------------------- Phase 3 arithmetic (new)
# Built once per curve. These are the fixed-width Montgomery routines that
# replace the heap mpz_t field layer; see include/elips/fp.h.
set(ELIPS_CURVES BLS12_381 BLS12_461 BN_462)
foreach(_cv IN LISTS ELIPS_CURVES)
  add_library(elips_arith_${_cv} STATIC
    src/arith/fp.c src/arith/fpx.c src/arith/ec.c
    src/arith/serialize.c src/util/random.c src/util/sysrand.c
    src/pairing/miller.c)
  target_include_directories(elips_arith_${_cv} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)
  target_compile_definitions(elips_arith_${_cv} PUBLIC ELIPS_CURVE_${_cv})
  target_link_libraries(elips_arith_${_cv} PUBLIC GMP::GMP)
  target_compile_options(elips_arith_${_cv} PRIVATE -Wall -Wextra)
  set_target_properties(elips_arith_${_cv} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endforeach()

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

  # Differential test: the new Montgomery layer against GMP's mpz, on random
  # inputs plus the boundary values that Montgomery code actually gets wrong.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(fp_difftest_${_cv} test/fp_difftest.c)
    target_link_libraries(fp_difftest_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    add_test(NAME fp.diff.${_cv} COMMAND fp_difftest_${_cv} 20000)
    set_tests_properties(fp.diff.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # The Phase 0 vectors run against the NEW layer too, so both implementations
  # are checked against the same independent oracle rather than against each
  # other. Vector file and compiled curve must match; the runner enforces it.
  set(_vec_BLS12_381 bls12_381)
  set(_vec_BLS12_461 bls12_461)
  set(_vec_BN_462    bn_462)
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(kat_new_${_cv} test/kat_runner_new.c)
    target_link_libraries(kat_new_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    add_test(NAME kat.new.${_cv}
             COMMAND kat_new_${_cv}
                     "${CMAKE_CURRENT_SOURCE_DIR}/test/kat/${_vec_${_cv}}.vec")
    set_tests_properties(kat.new.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # The new pairing against the value Phase 0 established as correct, plus
  # non-degeneracy and mu_r membership. BLS12-461 only: the reference side is
  # the legacy library, which initialises that curve at runtime.
  # The new layer on its own, no legacy library linked. All three curves,
  # including BLS12-381, which the legacy layer never supported at all.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(standalone_${_cv} test/standalone_test.c)
    target_link_libraries(standalone_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    add_test(NAME standalone.${_cv} COMMAND standalone_${_cv})
    set_tests_properties(standalone.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # Serialization, randomness and square roots. All three curves; the negative
  # cases are the point, so this runs everywhere rather than on one curve.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(serialize_test_${_cv} test/serialize_test.c)
    target_link_libraries(serialize_test_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    add_test(NAME serialize.${_cv} COMMAND serialize_test_${_cv})
    set_tests_properties(serialize.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # Timing-leakage tests. Statistical, so they are noisy by nature and their
  # sample sizes are kept small enough to run on every CI build; a developer
  # investigating a specific routine runs the binary directly with a larger n.
  # BLS12-381 only in CTest: the other two exercise the same code paths at
  # different widths, and three copies of a slow statistical test is not three
  # times the information.
  add_executable(dudect_test test/dudect_test.c)
  target_link_libraries(dudect_test PRIVATE elips_arith_BLS12_381 GMP::GMP m)

  # Built in every configuration so the harness itself is sanitized, but only
  # RUN in un-instrumented builds. Under Asan or Ubsan the measurement is of the
  # instrumentation, not of the routine: every memory access grows a shadow
  # lookup, which changes both the mean and the variance the t-test reads.
  # A "clean" result there would be evidence of nothing.
  if(CMAKE_BUILD_TYPE STREQUAL "Asan" OR CMAKE_BUILD_TYPE STREQUAL "Ubsan"
     OR CMAKE_BUILD_TYPE STREQUAL "Msan")
    message(STATUS "dudect timing tests built but not registered under ${CMAKE_BUILD_TYPE}")
  else()
    # control_vartime is the negative control: it calls fp_inv_vartime on secret
    # inputs and the test FAILS if no leak is found. A leakage detector that
    # cannot detect a leak says nothing about the targets it passes -- the same
    # argument that put kat.detects_corruption in the Phase 0 suite.
    foreach(_t fp_mul fp_add fp_cselect fp_inv ep_mul ep2_mul ep2_mul_glv
               control_vartime)
      add_test(NAME dudect.${_t} COMMAND dudect_test ${_t})
      set_tests_properties(dudect.${_t} PROPERTIES TIMEOUT 900)
    endforeach()
  endif()

  # BLS12-461 and BN-462: the reference side is the legacy library, which picks
  # its curve at runtime, so only the curves it can actually initialise. The
  # BLS12-381 field and curve layers are covered by the vector tests instead.
  foreach(_cv BLS12_461 BN_462)
    add_executable(pairing_test_${_cv} test/pairing_test.c)
    target_link_libraries(pairing_test_${_cv} PRIVATE elips_arith_${_cv} elips GMP::GMP)
    add_test(NAME pairing.matches_reference.${_cv} COMMAND pairing_test_${_cv})
    set_tests_properties(pairing.matches_reference.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # Group-law properties. Only the BLS12-461 build: the test borrows a genuine
  # G2 point from the legacy layer, which initialises that curve at runtime.
  add_executable(ec_test test/ec_test.c)
  target_link_libraries(ec_test PRIVATE elips_arith_BLS12_461 elips GMP::GMP)
  add_test(NAME ec.properties COMMAND ec_test)
  set_tests_properties(ec.properties PROPERTIES TIMEOUT 900)

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
