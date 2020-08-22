import pathlib
from setuptools import setup

# The directory containing this file
cwd = pathlib.Path(__file__).parent

# The text of the README file
README = (cwd / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="py2sambvca",
    version="1.0.0",
    description="Simple thin client to interface python scripts with SambVca catalytic pocket Fortran calculator.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/JacksonBurns/py2sambvca",
    author="Jackson Burns",
    license="GNU GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3"
    ],
    packages=["py2sambvca"],
    include_package_data=True
)