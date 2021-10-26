
import sys
import os
import io
from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ["-std=c++11"]
EXTRA_LINK_ARGS = []

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
    EXTRA_LINK_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

weights = cythonize([
    Extension("denovonear.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=["denovonear/weights.pyx",
            "src/weighted_choice.cpp",
            "src/simulate.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("denovonear.transcript",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "denovonear/transcript.pyx",
            "src/tx.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("denovonear.gencode",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "denovonear/gencode.pyx",
            "src/gencode.cpp",
            "src/gtf.cpp",
            "src/tx.cpp",
            "src/gzstream/gzstream.C",
        ],
        include_dirs=["src/", "src/gzstream"],
        libraries=["z"],
        language="c++"),
    Extension("denovonear.site_specific_rates",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "denovonear/site_specific_rates.pyx",
            "src/weighted_choice.cpp",
            "src/tx.cpp",
            "src/site_rates.cpp"],
        include_dirs=["src/"],
        language="c++"),
    ])

setup(name="denovonear",
    description='Package to examine de novo clustering',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version="0.9.5",
    author="Jeremy McRae",
    author_email="jeremy.mcrae@gmail.com",
    license="MIT",
    url='https://github.com/jeremymcrae/denovonear',
    packages=["denovonear"],
    install_requires=[
        'aiohttp >= 3.0',
        'scipy >= 0.9.0',
        'cython >= 0.19.0',
        'pyfaidx >= 0.5.8',
    ],
    package_data={"denovonear": ['data/rates.txt', 'weights.pxd']},
    entry_points={'console_scripts': ['denovonear = denovonear.__main__:main']},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.6',
    ext_modules=weights,
    test_suite="tests")
