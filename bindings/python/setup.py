"""Everything is in pyproject.toml; cffi_modules needs a setup.py to hang off.

The .S source-extension patch that the assembly backend needs lives in
build_ffi.py, next to the source list it applies to, so that building the
extension directly works the same way pip does.
"""
from setuptools import setup

setup(cffi_modules=["build_ffi.py:ffibuilder"])
