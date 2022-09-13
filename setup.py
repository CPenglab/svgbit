import re
from setuptools import setup, find_packages

VERSIONFILE = "svgbit/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." %
                       (VERSIONFILE, ))

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name="svgbit",
    version=verstr,
    author="CPenglab",
    author_email="chengpeng@ynu.edu.cn",
    discription="Find spatial variable genes for Spatial Trasncriptomics data.",
    long_description=long_description,
    url="https://github.com/CPenglab/svgbit",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "setuptools",
        "matplotlib>=3.5.1",
        "pandana>=0.6.1",
        "pysal>=2.4.0",
        "pillow>=9.2.0",
        "seaborn>=0.11.2",
    ],
    entry_points={'console_scripts': [
        'svgbit = svgbit.__main__:main',
    ]},
    python_requires=">=3.8.12",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ])
