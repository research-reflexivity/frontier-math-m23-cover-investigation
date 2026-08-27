#!/usr/bin/env sage-python
"""Build the optional FLINT Acb FFT bridge in this directory."""

from __future__ import annotations

import pkgconfig
from Cython.Build import cythonize
from setuptools import Extension, setup


flint = pkgconfig.parse("flint")
extension = Extension(
    "acb_fft",
    ["acb_fft.pyx"],
    include_dirs=flint["include_dirs"],
    library_dirs=flint["library_dirs"],
    libraries=flint["libraries"],
)

setup(
    name="m23-acb-fft",
    ext_modules=cythonize([extension], compiler_directives={"language_level": 3}),
)
