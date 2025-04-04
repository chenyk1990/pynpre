#!/usr/bin/env python
# -*- encoding: utf8 -*-
import io
import os

#from setuptools import find_packages
from setuptools import setup
from distutils.core import Extension
import numpy


long_description = """
Source code: https://github.com/chenyk1990/pynpre""".strip() 


def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")).read()

from distutils.core import Extension

nprec3d_module = Extension('npre3dcfun', sources=['pynpre/src/npre3d.c',
												'pynpre/src/npre_fxynpre.c',
												'pynpre/src/npre_alloc.c',
												'pynpre/src/npre_kissfft.c',
												'pynpre/src/npre_komplex.c',
												'pynpre/src/npre_conjgrad.c',
												'pynpre/src/npre_cdivn.c',
												'pynpre/src/npre_triangle.c',
												'pynpre/src/npre_trianglen.c',
												'pynpre/src/npre_ntriangle.c',
												'pynpre/src/npre_ntrianglen.c',		
												'pynpre/src/npre_decart.c',	
												'pynpre/src/npre_win.c',	
												'pynpre/src/npre_memcpy.c',			
												'pynpre/src/npre_fft1.c'],
                                                include_dirs=[numpy.get_include()])

ftfa_module = Extension('ftfacfun', sources=['pynpre/src/tf.c',
                                                'pynpre/src/npre_fxynpre.c',
                                                'pynpre/src/npre_alloc.c',
                                                'pynpre/src/npre_kissfft.c',
                                                'pynpre/src/npre_komplex.c',
                                                'pynpre/src/npre_conjgrad.c',
                                                'pynpre/src/npre_cdivn.c',
                                                'pynpre/src/npre_triangle.c',
                                                'pynpre/src/npre_trianglen.c',
                                                'pynpre/src/npre_ntriangle.c',
                                                'pynpre/src/npre_ntrianglen.c',
                                                'pynpre/src/npre_decart.c',
                                                'pynpre/src/npre_win.c',
                                                'pynpre/src/npre_memcpy.c',
                                                'pynpre/src/npre_fft1.c'],
                                                include_dirs=[numpy.get_include()])
                                                
setup(
    name="pynpre",
    version="0.0.2",
    license='MIT License',
    description="A python package of non-stationary predictive filtering for denoising and interpolation of multi-dimensional multi-channel seismic data",
    long_description=long_description,
    author="pynpre developing team",
    author_email="chenyk2016@gmail.com",
    url="https://github.com/chenyk1990/pynpre",
    ext_modules=[nprec3d_module,ftfa_module],
    packages=['pynpre'],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    keywords=[
        "seismology", "earthquake seismology", "exploration seismology", "array seismology", "denoising", "science", "engineering", "structure", "local slope", "filtering"
    ],
    install_requires=[
        "numpy", "scipy", "matplotlib"
    ],
    extras_require={
        "docs": ["sphinx", "ipython", "runipy"]
    }
)
