from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

version = '0.7.1'

setup(name='mergem',
      version=version,
      description='Merge two or more genome scale metabolic models.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/lobolab/mergem',
      author='Lobo lab',
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
