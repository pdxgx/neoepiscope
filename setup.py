#!/usr/bin/env python
# coding=utf-8

from __future__ import absolute_import, division, print_function
from setuptools import setup, find_packages
from distutils.core import Command
from setuptools.command.install import install

# Borrowed (with revisions) from https://stackoverflow.com/questions/17001010/
# how-to-run-unittest-discover-from-python-setup-py-test/21726329#21726329
class DiscoverTest(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import os
        import sys
        import unittest

        # get setup.py directory
        setup_file = sys.modules["__main__"].__file__
        test_dir = os.path.join(os.path.abspath(os.path.dirname(setup_file)), "tests")
        # use the default shared TestLoader instance
        test_loader = unittest.defaultTestLoader
        # use the basic test runner that outputs to sys.stderr
        test_runner = unittest.TextTestRunner()
        # automatically discover all tests
        # NOTE: only works for python 2.7 and later
        test_suite = test_loader.discover(test_dir)
        print(test_suite)
        # run the test suite
        test_runner.run(test_suite)


class DownloadDependencies(Command):
    # Wrapper to accommodate old-style class Command
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        from neoepiscope.download import NeoepiscopeDownloader

        downloader = NeoepiscopeDownloader()
        downloader.run()


setup(
    name="neoepiscope",
    version="0.2.1",
    description="comprehensive neoepitope prediction software",
    long_description=(
        "neoepiscope is software for predicting neoepitopes "
        "from DNA sequencing (DNA-seq) data. Where most neoepitope "
        "prediction software confines attention to neoepitopes arising "
        "from at most one somatic mutation, often just an SNV, "
        "neoepiscope uses assembled haplotype output of HapCUT2 "
        "(Edge et al., 2017) to also enumerate neoepitopes arising "
        "from more than one somatic mutation. neoepiscope also takes "
        "into account frameshifting from indels and permits "
        "personalizing the reference transcriptome using germline variants."
    ),
    url="http://github.com/ohsu-comp-bio/neoepiscope",
    download_url="https://github.com/pdxgx/neoepiscope/tarball/v0.2.1",
    author="Mary A. Wood, Austin Nguyen, Abhinav Nellore, Reid F. Thompson",
    author_email="thompsre@ohsu.edu",
    license="MIT",
    packages=["neoepiscope"],
    include_package_data=True,
    package_data={"neoepiscope": ["*.py", "*.pickle"]},
    zip_safe=True,
    install_requires=["intervaltree", "mhcflurry", "mhcnuggets", "numpy==1.14.5"],
    entry_points={"console_scripts": ["neoepiscope=neoepiscope:main"]},
    cmdclass={"download": DownloadDependencies, "test": DiscoverTest},
    keywords=["neoepitope", "neoantigen", "cancer", "immunotherapy"],
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
