# Negative control for a vector runner: corrupt one expected result and require
# a non-zero exit. A test oracle that cannot fail is not an oracle, so every
# vector file gets one of these.
#
#   RUNNER   the executable under test
#   SRC      the vector file
#   PATTERN  the record type to damage (default fp12_mul, the field runner's)
if(NOT DEFINED PATTERN)
  set(PATTERN "fp12_mul")
endif()

#   RELABEL  optional: instead of damaging an expected value, rename every
#            PATTERN record to this. Used for the subgroup file, whose records
#            carry no expected value -- the claim IS the label. Relabelling
#            "must be rejected" as "must be accepted" is the sharper control
#            there, because it exercises the one direction that can fail
#            silently.
file(READ "${SRC}" _txt)
if(DEFINED RELABEL)
  string(REGEX REPLACE "\n${PATTERN} " "\n${RELABEL} " _bad "${_txt}")
else()
  string(REGEX REPLACE "(\n${PATTERN} [^\n]*)= [0-9a-f]+" "\\1= deadbeef" _bad "${_txt}")
endif()
if(_bad STREQUAL _txt)
  message(FATAL_ERROR
    "corruption harness found no '${PATTERN}' record to alter in ${SRC}")
endif()

# The name carries the pattern so two of these can run concurrently.
set(_tmp "${CMAKE_CURRENT_BINARY_DIR}/corrupted-${PATTERN}.vec")
file(WRITE "${_tmp}" "${_bad}")
execute_process(COMMAND "${RUNNER}" "${_tmp}" RESULT_VARIABLE _rc OUTPUT_QUIET ERROR_QUIET)
file(REMOVE "${_tmp}")
if(_rc EQUAL 0)
  message(FATAL_ERROR "runner returned 0 on a corrupted vector file; the oracle cannot detect errors")
endif()
message(STATUS "runner correctly rejected corrupted vectors (exit ${_rc})")
