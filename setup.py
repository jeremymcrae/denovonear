
import sys
from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ["-std=c++0x"]

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS = ["-stdlib=libc++"]

module1 = cythonize([
    Extension("denovonear.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=["denovonear/cpp/weights.pyx",
            "denovonear/cpp/weighted_choice.cpp",
            "denovonear/cpp/simulate.cpp"],
        language="c++"),
    ])

setup (name="denovonear",
        description='Package to examine de novo clustering',
        version="0.1.1",
        author="Jeremy McRae",
        author_email="jeremy.mcrae@sanger.ac.uk",
        license="MIT",
        packages=["denovonear", "denovonear.gene_plot"],
        install_requires=['pysam >= 0.8.0',
                          'scipy >= 0.9.0',
                          'cairocffi >= 0.7.2',
                          'webcolors >= 1.5',
                          'cython >= 0.22.0'
        ],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: MIT License",
        ],
        ext_modules=module1,
        test_suite="tests")
