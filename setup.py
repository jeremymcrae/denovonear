
import sys
from setuptools import setup, Extension

EXTRA_COMPILE_ARGS = ["-std=c++0x"]

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS = ["-stdlib=libc++"]

class lazy_cythonize(list):
    def __init__(self, callback):
        self._list, self.callback = None, callback
    def c_list(self):
        if self._list is None: self._list = self.callback()
        return self._list
    def __iter__(self):
        for e in self.c_list(): yield e
    def __getitem__(self, ii): return self.c_list()[ii]
    def __len__(self): return len(self.c_list())

def extensions():
    from Cython.Build import cythonize
    ext = Extension("denovonear.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=["denovonear/weights.pyx",
            "denovonear/weighted_choice.cpp",
            "denovonear/simulate.cpp"],
        language="c++")
    return cythonize([ext])

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
                          'cython >= 0.19.0'
        ],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: MIT License",
        ],
        ext_modules=lazy_cythonize(extensions),
        test_suite="tests")
