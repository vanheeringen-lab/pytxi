#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.util import convert_path
import re

# get the version from __init__.py
with open(convert_path("pytxi/__init__.py")) as ver_file:
    match = next(
        re.finditer('__version__ = "(.*)"', ver_file.read(), re.MULTILINE)
    )
    version = match.group(1)

# use the readme as long description
with open("README.md") as f:
    long_description = f.read()

classifiers = [
    "Development Status :: 1 - Planning",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="pytxi",
    version=version,
    description="We want tximport in Python too!",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    data_files=[("", ["LICENSE", "README.md"])],
    entry_points={"console_scripts": ["pytxi=pytxi.cli:main"]},
    install_requires=["genomepy>=0.10"],
    author="Simon van Heeringen",
    author_email="simon.vanheeringen@gmail.com",
    url="https://github.com/vanheeringen-lab/pytxi",
    license="MIT",
    classifiers=classifiers,
)
