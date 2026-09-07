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

# ------------------------------------------------------------ curve selection
# The modern arithmetic needs its curve at compile time: FP_LIMBS has to be a
# constant for the loops to unroll and for the Phase 6 assembly to target a
# known width. So a consumer picks one curve and gets one library.
#
#   cmake -B build -DELIPS_CURVE=BLS12_381
#
# That is the target that gets installed and that find_package(ELiPS) hands
# back as ELiPS::arith. The test suite still builds all three, because the
# vectors have to cover all three.
set(ELIPS_CURVES BLS12_381 BLS12_461 BN_462)
set(ELIPS_CURVE "BLS12_381" CACHE STRING
    "Curve for the installed library: BLS12_381 (the standard one), BLS12_461 or BN_462")
set_property(CACHE ELIPS_CURVE PROPERTY STRINGS ${ELIPS_CURVES})
if(NOT ELIPS_CURVE IN_LIST ELIPS_CURVES)
  message(FATAL_ERROR
    "ELIPS_CURVE is '${ELIPS_CURVE}'. Choose one of: ${ELIPS_CURVES}")
endif()

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
# One library, built once per curve (issue #17 retired the second one).
#
# There used to be a runtime-curve mpz_t layer under src/*.c alongside this,
# carrying the public headers, defects A3 and A4, and the wrong-exponent final
# exponentiation of issue #16. It is gone; these are the only sources now.
set(ELIPS_ARITH_SOURCES
  src/arith/fp.c
  src/arith/fpx.c
  src/arith/ec.c
  src/arith/serialize.c
  src/pairing/miller.c
  src/hash/sha256.c
  src/hash/hash_to_curve.c
  src/util/random.c
  src/util/sysrand.c)

foreach(_cv IN LISTS ELIPS_CURVES)
  add_library(elips_arith_${_cv} STATIC ${ELIPS_ARITH_SOURCES})
  target_include_directories(elips_arith_${_cv} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>)
  target_compile_definitions(elips_arith_${_cv} PUBLIC ELIPS_CURVE_${_cv})
  target_link_libraries(elips_arith_${_cv} PUBLIC GMP::GMP)
  target_compile_options(elips_arith_${_cv} PRIVATE -Wall -Wextra)
  if(ELIPS_WERROR)
    target_compile_options(elips_arith_${_cv} PRIVATE -Werror)
  endif()
  set_target_properties(elips_arith_${_cv} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})
endforeach()

# The selected curve, under a stable name. This is what gets installed, so a
# downstream project writes target_link_libraries(app PRIVATE ELiPS::arith)
# and never names a curve in its own build files.
add_library(ELiPS::arith ALIAS elips_arith_${ELIPS_CURVE})
message(STATUS "ELiPS curve for the installed library: ${ELIPS_CURVE}")

# --------------------------------------------------------------------- tests
if(ELIPS_BUILD_TESTS)
  enable_testing()

  # Every test target gets the same warning treatment as the library. Without
  # this, "zero warnings" would describe half the code.
  function(elips_test_warnings tgt)
    target_compile_options(${tgt} PRIVATE -Wall -Wextra)
    if(ELIPS_WERROR)
      target_compile_options(${tgt} PRIVATE -Werror)
    endif()
  endfunction()

  # Curve-name to vector-file mapping, used by several suites below.
  set(_vec_BLS12_381 bls12_381)
  set(_vec_BLS12_461 bls12_461)
  set(_vec_BN_462    bn_462)

  # Differential test: the new Montgomery layer against GMP's mpz, on random
  # inputs plus the boundary values that Montgomery code actually gets wrong.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(fp_difftest_${_cv} test/fp_difftest.c)
    target_link_libraries(fp_difftest_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(fp_difftest_${_cv})
    add_test(NAME fp.diff.${_cv} COMMAND fp_difftest_${_cv} 20000)
    set_tests_properties(fp.diff.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # The Phase 0 field and curve vectors, from the independent Python oracle.
  # Vector file and compiled curve must match; the runner enforces it.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(kat_runner_${_cv} test/kat_runner.c)
    target_link_libraries(kat_runner_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(kat_runner_${_cv})
    add_test(NAME kat.${_cv}
             COMMAND kat_runner_${_cv}
                     "${CMAKE_CURRENT_SOURCE_DIR}/test/kat/${_vec_${_cv}}.vec")
    set_tests_properties(kat.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # Negative control: the runner must actually fail on a corrupted vector.
  # A test oracle that cannot fail is not an oracle, so this is not optional.
  add_test(NAME kat.detects_corruption
           COMMAND ${CMAKE_COMMAND}
                   -DRUNNER=$<TARGET_FILE:kat_runner_BLS12_381>
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/test/kat/bls12_381.vec
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorruptCheck.cmake)

  # All three curves, with no second implementation linked in.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(standalone_${_cv} test/standalone_test.c)
    target_link_libraries(standalone_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(standalone_${_cv})
    add_test(NAME standalone.${_cv} COMMAND standalone_${_cv})
    set_tests_properties(standalone.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # Serialization, randomness and square roots. All three curves; the negative
  # cases are the point, so this runs everywhere rather than on one curve.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(serialize_test_${_cv} test/serialize_test.c)
    target_link_libraries(serialize_test_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(serialize_test_${_cv})
    add_test(NAME serialize.${_cv} COMMAND serialize_test_${_cv})
    set_tests_properties(serialize.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # Hash-to-curve. On BLS12-381 these records use RFC 9380's own domain
  # separation tag, so they ARE the specification's published vectors and this
  # test is what pins interoperability. The other two curves have no registered
  # suite, so theirs pin agreement with the Python reference.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(h2c_test_${_cv} test/h2c_test.c)
    target_link_libraries(h2c_test_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(h2c_test_${_cv})
    add_test(NAME h2c.${_cv}
             COMMAND h2c_test_${_cv}
                     "${CMAKE_CURRENT_SOURCE_DIR}/test/kat/h2c_${_vec_${_cv}}.vec")
    set_tests_properties(h2c.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # Negative control, same argument as kat.detects_corruption: a runner that
  # cannot fail proves nothing about the vectors it passes.
  add_test(NAME h2c.detects_corruption
           COMMAND ${CMAKE_COMMAND}
                   -DRUNNER=$<TARGET_FILE:h2c_test_BLS12_381>
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/test/kat/h2c_bls12_381.vec
                   -DPATTERN=h2c_g1
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorruptCheck.cmake)

  # Edge cases and API contracts: boundaries, the identity, degenerate group-law
  # inputs, and aliasing. Aliasing has its own section because two of the
  # original audit's defects were exactly that, and no vector file would catch
  # one -- vectors never alias their inputs and outputs.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(edge_test_${_cv} test/edge_test.c)
    target_link_libraries(edge_test_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(edge_test_${_cv})
    add_test(NAME edge.${_cv} COMMAND edge_test_${_cv})
    set_tests_properties(edge.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # Timing-leakage tests. Statistical, so they are noisy by nature and their
  # sample sizes are kept small enough to run on every CI build; a developer
  # investigating a specific routine runs the binary directly with a larger n.
  # BLS12-381 only in CTest: the other two exercise the same code paths at
  # different widths, and three copies of a slow statistical test is not three
  # times the information.
  add_executable(dudect_test test/dudect_test.c)
  target_link_libraries(dudect_test PRIVATE elips_arith_BLS12_381 GMP::GMP m)
  elips_test_warnings(dudect_test)

  # BN-462 as well, and the reason is a real exception to the paragraph above.
  # BN G1 GLV is not the BLS12 code at a different width: it rounds against a
  # lattice basis, with wide multiplies, a rounded division and signed digits
  # carried in two's complement, none of which BLS12-381 ever executes. So it
  # gets its own measurement, and its own negative controls, because a clean
  # reading from a harness that has not been shown to work on this build says
  # nothing.
  add_executable(dudect_test_bn test/dudect_test.c)
  target_link_libraries(dudect_test_bn PRIVATE elips_arith_BN_462 GMP::GMP m)
  elips_test_warnings(dudect_test_bn)

  # Built in every configuration so the harness itself is sanitized, but only
  # RUN in un-instrumented builds. Under Asan or Ubsan the measurement is of the
  # instrumentation, not of the routine: every memory access grows a shadow
  # lookup, which changes both the mean and the variance the t-test reads.
  # A "clean" result there would be evidence of nothing.
  if(CMAKE_BUILD_TYPE STREQUAL "Asan" OR CMAKE_BUILD_TYPE STREQUAL "Ubsan"
     OR CMAKE_BUILD_TYPE STREQUAL "Msan")
    message(STATUS "dudect timing tests built but not registered under ${CMAKE_BUILD_TYPE}")
  else()
    # Two negative controls, at opposite ends, both of which FAIL if no leak is
    # found. A leakage detector that cannot detect a leak says nothing about the
    # targets it passes -- the same argument that put kat.detects_corruption in
    # the Phase 0 suite.
    #
    #   control_vartime  calls fp_inv_vartime on secrets; reads in the hundreds,
    #                    so it proves the harness finds an obvious leak.
    #   sensitivity      plants ONE extra field multiplication on one bit of the
    #                    secret; reads 30 to 60, so it proves the harness finds a
    #                    small one. It is also what keeps the decision threshold
    #                    honest: raise the threshold past a one-multiply leak and
    #                    this test starts failing.
    foreach(_t fp_mul fp_add fp_cselect fp_inv ep_mul ep2_mul ep2_mul_glv
               ep_mul_glv hash_to_g1 hash_to_g2 sensitivity control_vartime)
      add_test(NAME dudect.${_t} COMMAND dudect_test ${_t})
      set_tests_properties(dudect.${_t} PROPERTIES TIMEOUT 900)
    endforeach()

    # Only the routines whose code BN does not share, plus the controls.
    foreach(_t ep_mul_glv ep2_mul_glv sensitivity control_vartime)
      add_test(NAME dudect.bn.${_t} COMMAND dudect_test_bn ${_t})
      set_tests_properties(dudect.bn.${_t} PROPERTIES TIMEOUT 900)
    endforeach()
  endif()

  # The pairing, against the independent Python oracle in
  # tools/reference/pairing_ref.py rather than against a second implementation.
  # All three curves, including BLS12-381, which had no pairing reference at all
  # while the legacy layer was the reference side (issue #17).
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(pairing_test_${_cv} test/pairing_test.c)
    target_link_libraries(pairing_test_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(pairing_test_${_cv})
    add_test(NAME pairing.${_cv}
             COMMAND pairing_test_${_cv}
                     "${CMAKE_CURRENT_SOURCE_DIR}/test/kat/pairing_${_vec_${_cv}}.vec")
    set_tests_properties(pairing.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  # And the pairing runner must be able to fail too.
  add_test(NAME pairing.detects_corruption
           COMMAND ${CMAKE_COMMAND}
                   -DRUNNER=$<TARGET_FILE:pairing_test_BLS12_381>
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/test/kat/pairing_bls12_381.vec
                   -DPATTERN=pair
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorruptCheck.cmake)

  # Subgroup membership, in both directions. The rejection direction is the
  # reason this exists: a test that is too permissive accepts attacker-chosen
  # points of small order while every other test still passes.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(subgroup_test_${_cv} test/subgroup_test.c)
    target_link_libraries(subgroup_test_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(subgroup_test_${_cv})
    add_test(NAME subgroup.${_cv}
             COMMAND subgroup_test_${_cv}
                     "${CMAKE_CURRENT_SOURCE_DIR}/test/kat/subgroup_${_vec_${_cv}}.vec")
    set_tests_properties(subgroup.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

  add_test(NAME subgroup.detects_corruption
           COMMAND ${CMAKE_COMMAND}
                   -DRUNNER=$<TARGET_FILE:subgroup_test_BLS12_381>
                   -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/test/kat/subgroup_bls12_381.vec
                   -DPATTERN=outg2 -DRELABEL=ing2
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorruptCheck.cmake)

  # Benchmarks (issue #15).
  #
  # Built for every curve, but NOT registered with CTest. A benchmark that runs
  # in CI becomes a test that fails when the runner is busy, and the point of
  # this harness is to tell a real regression from a noisy machine rather than
  # to conflate the two. Run it deliberately:
  #
  #   cmake --build build --target bench_BLS12_381
  #   ./build/bench_BLS12_381 --reps 15
  #   ./build/bench_BLS12_381 --json > before.json
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(bench_${_cv} bench/bench.c)
    target_link_libraries(bench_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(bench_${_cv})
  endforeach()

  # Group-law properties, built for every curve. The GLV routines differ by family (four
  # dimensions on BLS12 G2, two on BN G2 and on BLS12 G1), so a single-curve
  # build would leave most of them untested.
  foreach(_cv IN LISTS ELIPS_CURVES)
    add_executable(ec_test_${_cv} test/ec_test.c)
    target_link_libraries(ec_test_${_cv} PRIVATE elips_arith_${_cv} GMP::GMP)
    elips_test_warnings(ec_test_${_cv})
    add_test(NAME ec.properties.${_cv} COMMAND ec_test_${_cv})
    set_tests_properties(ec.properties.${_cv} PROPERTIES TIMEOUT 900)
  endforeach()

endif()

# ------------------------------------------------------------------ examples
# Runnable demonstrations of the consumer-facing API, built for the selected
# curve. They are also tests: CI runs each one and requires exit 0, so an API
# change that breaks a caller breaks the build rather than rotting quietly in a
# README.
option(ELIPS_BUILD_EXAMPLES "Build the example programs" ON)
if(ELIPS_BUILD_EXAMPLES)
  set(ELIPS_EXAMPLES pairing hash bls)
  set(_ex_pairing 01_pairing.c)
  set(_ex_hash    02_hash_and_serialize.c)
  set(_ex_bls     03_bls_signature.c)
  foreach(_ex IN LISTS ELIPS_EXAMPLES)
    add_executable(elips_example_${_ex} examples/${_ex_${_ex}})
    target_link_libraries(elips_example_${_ex} PRIVATE elips_arith_${ELIPS_CURVE})
    target_compile_options(elips_example_${_ex} PRIVATE -Wall -Wextra)
    if(ELIPS_WERROR)
      target_compile_options(elips_example_${_ex} PRIVATE -Werror)
    endif()
    set_target_properties(elips_example_${_ex} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples")
    if(ELIPS_BUILD_TESTS)
      add_test(NAME example.${_ex} COMMAND elips_example_${_ex})
      set_tests_properties(example.${_ex} PROPERTIES TIMEOUT 900)
    endif()
  endforeach()
endif()

# ------------------------------------------------------------------- install
include(GNUInstallDirs)
install(TARGETS elips_arith_${ELIPS_CURVE} EXPORT elipsTargets
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/elips
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h")
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

# The library is built per curve, so its exported name carries one.
# ELiPS::arith is the stable alias a consumer links against, and
# ELiPS_CURVE says which curve this installation was built for.
set(ELiPS_CURVE \"${ELIPS_CURVE}\")
if(NOT TARGET ELiPS::arith)
  add_library(ELiPS::arith ALIAS ELiPS::elips_arith_${ELIPS_CURVE})
endif()
# ELiPS::elips named the runtime-curve mpz layer, retired by issue #17. Kept as
# an alias so an existing consumer still configures; it now resolves to the
# same library as ELiPS::arith.
if(NOT TARGET ELiPS::elips)
  add_library(ELiPS::elips ALIAS ELiPS::elips_arith_${ELIPS_CURVE})
endif()
")
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/ELiPSConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/ELiPSConfigVersion.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGMP.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ELiPS)

message(STATUS "ELiPS ${PROJECT_VERSION}  build=${CMAKE_BUILD_TYPE}  GMP=${GMP_VERSION}")
