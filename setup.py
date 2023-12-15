from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

exec(open("./mergem/__version.py").read())

setup(name='mergem',
      version=_version,
      description='mergem is a python package and command-line tool for merging, comparing, and translating genome-scale metabolic models.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://lobolab.umbc.edu/mergem/',
      author='Lobo Lab',
      author_email='lobolab@umbc.edu',
      maintainer='Archana Hari',
      maintainer_email="archh1@umbc.edu",
      license='GNU GPLv3',
      packages=['mergem'],
      entry_points={
            'console_scripts': [
                  'mergem = mergem.cli:main'
            ]
      },
      classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Operating System :: OS Independent",
            "Intended Audience :: Science/Research",
            "Development Status :: 4 - Beta",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
      ],
      install_requires=['cobra >= 0.15.4', 'click>=8.0.3', 'requests'],
      include_package_data=True,
      zip_safe=False)
