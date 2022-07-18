from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="svgene",
    version="0.0.4",
    long_description=long_description,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "setuptools",
        "pysal>=2.4.0",
        "pandana>=0.6.1",
    ],
    entry_points={'console_scripts': [
        'svgene = svgene.__main__:main',
    ]},
    py_modules=[],
    license="MIT",
    python_requires=">=3.8.12",
)
