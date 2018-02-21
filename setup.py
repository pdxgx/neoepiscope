from setuptools import setup

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
