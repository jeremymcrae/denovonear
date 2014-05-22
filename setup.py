from distutils.core import setup, Extension

module1 = Extension('libweightedchoice',
          extra_compile_args = ["-std=c++0x"],
          sources = ['./src/weighted_choice.cpp'])

setup (name = 'De novo clustering',
       description = 'Package to examine de novo clustering',
       version = "0.1",
       author = 'Jeremy McRae',
       author_email = 'jeremy.mcrae@sanger.ac.uk',
       url = 'http://www.sanger.ac.uk',
       long_description = '''
This is really just a demo package.
''',
       ext_modules = [module1])