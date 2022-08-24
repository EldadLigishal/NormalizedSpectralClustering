from distutils.core import setup, Extension

setup(
    name='mySpkmeans',
    author='Yara and Eldad',
    version='1.0',
    ext_modules=[Extension('mySpkmeans',sources=['spkmeans.c', 'spkmeansmodule.c'])])
