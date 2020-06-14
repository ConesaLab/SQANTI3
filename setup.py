#!/usr/bin/env python
import sys

if sys.version_info < (3,):
    sys.exit("SQANTI3 requires Python >= 3.7")
from pathlib import Path
from setuptools import setup, find_packages

try:
    from sqanti3 import __author__, __email__
except ImportError:  # Deps not yet installed
    __author__ = __email__ = ""

setup(
    name="SQANTI3",
    version="1.0.1",
    description="Quality Control of Long-Read Defined Transcriptomes",
    long_description=Path("README.MD").read_text("utf-8"),
    url="https://github.com/ConesaLab/SQANTI3",
    author=__author__,
    author_email=__email__,
    license="GPL3",
    python_requires=">=3.7",
    install_requires=[
        l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    classifiers=[
        "Development Status :: 3 - Production",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    # extras_require=dict(doc=["sphinx", "sphinx_rtd_theme", "sphinx_autodoc_typehints"]),
    packages=find_packages(),
    include_package_data=True,
    entry_points={"console_scripts": ["sqanti3_qc = sqanti3.sqanti3_qc:main",
                                      "sqanti3_RulesFilter = sqanti3.sqanti3_RulesFilter.py:main"]},
    keywords="",
    package_dir={"sqanti3": "sqanti3"},
    package_data={"sqanti3": ["utilities/**", "GMST", "RTS"]},
)
