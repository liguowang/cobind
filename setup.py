import sys, os, platform, glob
from distutils.core import setup
from setuptools import *

"""
Setup script for cobind  -- collocation analysis of genomics regions.
"""


def main():
    setup(  name = "cobind",
            version = "0.0.5",
            python_requires='>=3.5',
            py_modules = [ 'psyco_full' ],
            packages = find_packages( 'lib' ),
            package_dir = { '': 'lib' },
            package_data = { '': ['*.ps'] },
            scripts = glob.glob( "bin/*.py"),
            ext_modules = [],
            test_suite = 'nose.collector',
            setup_requires = ['nose>=0.10.4','cython>=0.12'],
            author = "Liguo Wang",
            author_email ="wangliguo78@gmail.com",
            platforms = ['Linux','MacOS'],
            requires = [],
            install_requires = ['scipy', 'numpy', 'pandas', 'bx-python','pyBigWig'],
            description = "collocation analysis of genomics regions",
            url = "",
            zip_safe = False,
            dependency_links = [],
            classifiers=[
                'Development Status :: 4 - Development',
                'Environment :: Console',
                'Intended Audience :: Science/Research',
                'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
                'Operating System :: MacOS :: MacOS X',
                'Operating System :: POSIX',
                'Programming Language :: Python',
                'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],

            keywords='genome regions, overlap coefficient, cooccur, pointwise mutual information',
             )


if __name__ == "__main__":
    main()
