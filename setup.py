
import sys
import os
import shutil
from pathlib import Path

from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize
import gencodegenes

EXTRA_COMPILE_ARGS = []
EXTRA_LINK_ARGS = []
if sys.platform == "win32":
    EXTRA_COMPILE_ARGS += ['/std:c++14']
else:
    EXTRA_COMPILE_ARGS += ['-std=c++11']
    if sys.platform == "darwin":
        EXTRA_COMPILE_ARGS += ["-stdlib=libc++"]
        EXTRA_LINK_ARGS += ["-stdlib=libc++"]

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
        include_dirs=["src/", os.path.dirname(gencodegenes.__file__)],
        language="c++"),
    ]

setup(
    package_dir={'': 'src'},
    package_data={"denovonear": ['data/rates.txt', 'weights.pxd']},
    ext_modules=cythonize(extensions),
    test_suite="tests")
