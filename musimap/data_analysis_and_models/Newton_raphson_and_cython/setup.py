from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Cython codes',
    ext_modules=cythonize("cython_codes.pyx"),
    zip_safe=False,
)