
import sys
import os
import io
import glob

from distutils.ccompiler import new_compiler
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

def build_zlib():
    ''' compile zlib code to object files for linking on windows
    
    Returns:
        list of paths to compiled object code
    '''
    include_dirs = ['src/zlib/']
    sources = list(glob.glob('src/zlib/*.c'))
    extra_compile_args = []
    
    cc = new_compiler()
    return cc.compile(sources, include_dirs=include_dirs,
        extra_preargs=extra_compile_args)

def get_gzstream_path():
    ''' workaround for building gzstream on windows
    cython on windows didn't like the .C extension for gzstream. This just
    renames the file (on windows only), and returns the relative path.
    '''
    gzstream_path = 'src/gzstream/gzstream.C'
    if sys.platform == 'win32':
        gzstream_win_path = 'src/gzstream/gzstream.cpp'
        try:
            os.rename(gzstream_path, gzstream_win_path)
        except FileNotFoundError:
            pass  # avoid error on github actions
        gzstream_path = gzstream_win_path
    return gzstream_path

def scrub_gzstream():
    ''' workaround for compilation error on macos
    
    compiling gzstream requires the corresponding gzstream.h file, but if we 
    include the gzstream directory in the include dirs, then clang complains
    about the version file in the gzstream folder. If we remove the gzstream
    directory from the include dirs, then clang complains about the missing
    gzstream.h. This is because gzstream.C identifies it's header file with
    angle brackets. Replacing the angle brackets in that line seems to work.
    '''
    with open(get_gzstream_path(), 'rt') as handle:
        lines = handle.readlines()
    
    with open(get_gzstream_path(), 'wt') as handle:
        for line in lines:
            if line == '#include <gzstream.h>\n':
                line = '#include "gzstream.h"\n'
            handle.write(line)

if sys.platform == 'win32':
    zlib, libs = build_zlib(), []
else:
    zlib, libs = [], ['z']

scrub_gzstream()

weights = cythonize([
    Extension("denovonear.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=["src/denovonear/weights.pyx",
            "src/weighted_choice.cpp",
            "src/simulate.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("denovonear.transcript",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "src/denovonear/transcript.pyx",
            "src/tx.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("denovonear.gencode",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "src/denovonear/gencode.pyx",
            "src/gencode.cpp",
            "src/gtf.cpp",
            "src/tx.cpp",
            get_gzstream_path(),
        ],
        include_dirs=["src/"],
        libraries=["z"],
        language="c++"),
    Extension("denovonear.site_specific_rates",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "src/denovonear/site_specific_rates.pyx",
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
    version="0.9.15",
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
