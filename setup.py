
import sys
import os
import io
from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ["-std=c++11"]

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += ["-stdlib=libc++"]

weights = cythonize([
    Extension("denovonear.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=["denovonear/weights.pyx",
            "src/weighted_choice.cpp",
            "src/simulate.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("denovonear.transcript",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=[
            "denovonear/transcript.pyx",
            "src/tx.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("denovonear.site_specific_rates",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=[
            "denovonear/site_specific_rates.pyx",
            "src/weighted_choice.cpp",
            "src/tx.cpp",
            "src/site_rates.cpp"],
        include_dirs=["src/"],
        language="c++"),
    ])

setup (name="denovonear",
        description='Package to examine de novo clustering',
        long_description=io.open('README.rst', encoding='utf-8').read(),
        version="0.6.4",
        author="Jeremy McRae",
        author_email="jeremy.mcrae@sanger.ac.uk",
        license="MIT",
        url='https://github.com/jeremymcrae/denovonear',
        packages=["denovonear", "denovonear.gene_plot"],
        install_requires=['scipy >= 0.9.0',
                          'cairocffi >= 0.7.2',
                          'webcolors >= 1.5',
                          'cython >= 0.19.0'
        ],
        package_data={"denovonear": ['data/rates.txt', 'weights.pxd']},
        entry_points={'console_scripts': ['denovonear = denovonear.__main__:main']},
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "License :: OSI Approved :: MIT License",
        ],
        ext_modules=weights,
        test_suite="tests")
