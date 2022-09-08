from distutils.core import setup, Extension

setup(
    name='myspkmeans',
    author='Yara and Eldad',
    version='1.0',
    description='Python wrapper for custom C extension',
    ext_modules=[Extension('myspkmeans', sources=['spkmeans.c', 'spkmeansmodule.c'])])
