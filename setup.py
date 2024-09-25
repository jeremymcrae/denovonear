
import sys
import os
import io
import glob
from pathlib import Path

from distutils.ccompiler import new_compiler
from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize
import gencodegenes

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

# Find the path to a gencodegenes code file required for site_rates.cpp
# Or we could put the file directly in this package (or via a git submodule).
tx_cpp = Path(gencodegenes.__file__).parent / 'tx.cpp'

weights = cythonize([
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
            "src/site_rates.cpp",
            str(tx_cpp),
            ],
        include_dirs=["src/"],
        language="c++"),
    ])

setup(name="denovonear",
    description='Package to examine de novo clustering',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version="0.11.1",
    author="Jeremy McRae",
    author_email="jeremy.mcrae@gmail.com",
    license="MIT",
    url='https://github.com/jeremymcrae/denovonear',
    packages=["denovonear"],
    install_requires=[
        'aiohttp >= 3.9.0b0',
        'scipy >= 0.9.0',
        'gencodegenes >= 1.0.10',
    ],
    extras_requires={
        'pysam': ['pysam'],
    },
    package_dir={'': 'src'},
    package_data={"denovonear": ['data/rates.txt', 'weights.pxd']},
    entry_points={'console_scripts': ['denovonear = denovonear.__main__:main']},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.8',
    ext_modules=weights,
    test_suite="tests")
