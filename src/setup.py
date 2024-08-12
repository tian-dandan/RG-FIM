import pybind11
from distutils.core import setup, Extension
from distutils import ccompiler

cpp_args = [ ]
print("Compiler is: ", ccompiler.get_default_compiler())
sfc_module = Extension(
    'ImageObject',
    sources=['ImageObject.cpp','pybind11.cpp'],
    
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
    ],
    language='c++',
    extra_compile_args=cpp_args)


setup(
    name='ImageObject',
    version='1.1',
    description='Python package for image object operations',
    ext_modules=[sfc_module],
)
