
import sys
from setuptools import setup
from distutils.core import Extension

EXTRA_COMPILE_ARGS = ["-std=c++0x"]

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS = ["-stdlib=libc++"]

module1 = Extension("libsimulatedenovo",
        extra_compile_args = EXTRA_COMPILE_ARGS,
        sources = ["denovonear/cpp/simulate.cpp", "denovonear/cpp/weighted_choice.cpp"])

setup (name = "denovonear",
        description = 'Package to examine de novo clustering',
        version = "0.1.1",
        author = "Jeremy McRae",
        author_email = "jeremy.mcrae@sanger.ac.uk",
        license="MIT",
        packages=["denovonear"],
        install_requires=['numpy >= 1.6.1',
                          'scipy >= 0.9.0'
        ],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: MIT License",
        ],
        ext_modules = [module1],
        test_suite="tests")
