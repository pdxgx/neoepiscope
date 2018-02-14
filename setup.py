from setuptools import setup

setup(name='neoepiscope',
      version='0.1.0',
      description='comprehensive neoepitope prediction software',
      url='http://github.com/ohsu-comp-bio/neoepiscope',
      download_url = 'https://github.com/nellore/vickitrix/tarball/0.1.5',
      author='Mary Wood, Julianne David, Austin Nguyen, Mihir Paralkar, Mayur Paralkar, Abhinav Nellore, Reid Thompson',
      author_email='thompsre@ohsu.edu',
      license='MIT',
      packages=['neoepiscope'],
      package_data={'neoepiscope': ['*.py']},
      zip_safe=True,
      install_requires=[
      		'twython', 'gdax', 'pycrypto'
      	],
      entry_points={
        'console_scripts': [
            'vickitrix=vickitrix:go',
        ],},
      keywords=['bitcoin', 'btc', 'ethereum', 'eth', 'twitter'],
      classifiers=[
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 3',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Operating System :: MacOS',
          'Operating System :: Unix',
          'Topic :: Utilities'
        ]
    )
