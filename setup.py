from setuptools import setup, find_packages

# Borrowed from https://stackoverflow.com/questions/17001010/
# how-to-run-unittest-discover-from-python-setup-py-test/21726329#21726329
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

def discover_and_run_tests():
    import os
    import sys
    import unittest

    # get setup.py directory
    setup_file = sys.modules['__main__'].__file__
    setup_dir = os.path.abspath(os.path.dirname(setup_file))

    # use the default shared TestLoader instance
    test_loader = unittest.defaultTestLoader

    # use the basic test runner that outputs to sys.stderr
    test_runner = unittest.TextTestRunner()

    # automatically discover all tests
    # NOTE: only works for python 2.7 and later
    test_suite = test_loader.discover(setup_dir)
    print test_suite

    # run the test suite
    test_runner.run(test_suite)

try:
    from setuptools.command.test import test

    class DiscoverTest(test):

        def finalize_options(self):
            test.finalize_options(self)
            self.test_args = []
            self.test_suite = True

        def run_tests(self):
            discover_and_run_tests()

except ImportError:
    from distutils.core import Command

    class DiscoverTest(Command):
        user_options = []

        def initialize_options(self):
                pass

        def finalize_options(self):
            pass

        def run(self):
            discover_and_run_tests()

setup(name='neoepiscope',
      version='0.1.0',
      description='comprehensive neoepitope prediction software',
      url='http://github.com/ohsu-comp-bio/neoepiscope',
      download_url = 'https://github.com/ohsu-comp-bio/neoepiscope/tarball/0.1.0',
      author='Mary A. Wood, Julianne David, Austin Nguyen, Abhinav Nellore, Reid F. Thompson',
      author_email='thompsre@ohsu.edu',
      license='MIT',
      packages=['neoepiscope'],
      package_data={'neoepiscope': ['*.py']},
      zip_safe=True,
      install_requires=[
      		'intervaltree', 'mhcflurry'
      	],
      entry_points={
        'console_scripts': [
            'neoepiscope=neoepiscope:main',
        ],},
      cmdclass={'test': DiscoverTest},
      keywords=['neoepitope', 'neoantigen', 'cancer', 'immunotherapy'],
      classifiers=[
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 3',
          'License :: OSI Approved :: MIT License',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Education',
          'Operating System :: MacOS',
          'Operating System :: Unix',
          'Operating System :: Windows',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
    )
