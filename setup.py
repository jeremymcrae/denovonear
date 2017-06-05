
import sys
from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ["-std=c++11"]

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS = ["-stdlib=libc++"]

weights = cythonize([
    Extension("denovonear.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=["denovonear/weights.pyx",
            "denovonear/weighted_choice.cpp",
            "denovonear/simulate.cpp"],
        language="c++"),
    Extension("denovonear.transcript",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=[
            "denovonear/transcript.pyx",
            "denovonear/tx.cpp"],
        language="c++"),
    Extension("denovonear.site_specific_rates",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=[
            "denovonear/site_specific_rates.pyx",
            "denovonear/weighted_choice.cpp",
            "denovonear/tx.cpp",
            "denovonear/site_rates.cpp"],
        language="c++"),
    ])

setup (name="denovonear",
        description='Package to examine de novo clustering',
        version="0.4.0",
        author="Jeremy McRae",
        author_email="jeremy.mcrae@sanger.ac.uk",
        license="MIT",
        packages=["denovonear", "denovonear.gene_plot"],
        install_requires=['scipy >= 0.9.0',
                          'cairocffi >= 0.7.2',
                          'webcolors >= 1.5',
                          'cython >= 0.19.0'
        ],
        package_data={"denovonear": ['data/rates.txt']},
        classifiers=[
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: MIT License",
        ],
        ext_modules=weights,
        test_suite="tests")
