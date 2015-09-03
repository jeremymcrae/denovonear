
import sys
from distutils.core import setup, Extension

EXTRA_COMPILE_ARGS = ["-std=c++0x"]

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS = ["-stdlib=libc++"]

module1 = Extension("libsimulatedenovo",
        extra_compile_args = EXTRA_COMPILE_ARGS,
        sources = ["src/cpp/simulate.cpp", "src/cpp/weighted_choice.cpp"])

setup (name = "De novo clustering",
        description = 'Package to examine de novo clustering',
        version = "0.1",
        author = "Jeremy McRae",
        author_email = "jeremy.mcrae@sanger.ac.uk",
        ext_modules = [module1])
