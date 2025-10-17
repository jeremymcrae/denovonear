
import sys
import os
import shutil
from pathlib import Path

from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ['-std=c++11', '-I/usr/include']
EXTRA_LINK_ARGS = []
if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += [
        "-stdlib=libc++",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include",
        ]
    EXTRA_LINK_ARGS += [
        "-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib",
        ]
elif sys.platform == "win32":
    EXTRA_COMPILE_ARGS += ['/std:c++14']

extensions = [
    Extension("denovonear.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=["src/denovonear/weights.pyx",
            "src/weighted_choice.cpp",
            "src/simulate.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("denovonear.site_specific_rates",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "src/denovonear/site_specific_rates.pyx",
            "src/weighted_choice.cpp",
            ],
        include_dirs=["src/"],
        language="c++"),
    ]

setup(
    package_dir={'': 'src'},
    package_data={"denovonear": ['data/rates.txt', 'weights.pxd']},
    ext_modules=cythonize(extensions),
    test_suite="tests")
