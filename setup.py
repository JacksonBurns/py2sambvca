import os.path
import codecs
import pathlib
from setuptools import setup

# The directory containing this file
cwd = pathlib.Path(__file__).parent

# The text of the README file
README = (cwd / "README.md").read_text()


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


# This call to setup() does all the work
setup(
    name="py2sambvca",
    version=get_version("py2sambvca/__init__.py"),
    description="Simple thin client to interface python scripts with SambVca catalytic pocket Fortran calculator.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/GwydionJon/py2sambvca",
    author="Jackson Burns",
    license="GNU GPLv3",
    classifiers=["Programming Language :: Python :: 3"],
    packages=["py2sambvca"],
    include_package_data=True,
)
