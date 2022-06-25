from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="SVGLib",
    version="0.0.1",
    long_description=long_description,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "setuptools",
        "pysal>=2.4.0",
        "pandana>=0.6.1",
    ],
    py_modules=[],
    license="MIT",
    python_requires=">=3.8.12",
)
