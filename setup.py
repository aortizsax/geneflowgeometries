#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re
from setuptools import setup, find_packages


def _read(path_components, **kwargs):
    path = os.path.join(os.path.dirname(__file__), *path_components)
    if sys.version_info.major < 3:
        return open(path, "rU").read()
    else:
        with open(path, encoding=kwargs.get("encoding", "utf8")) as src:
            s = src.read()
        return s


def _read_requirements(path):
    return [
        line.strip()
        for line in _read([path]).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


project_init = _read(["src", "geneflowgeometries", "__init__.py"])
__version__ = re.match(
    r".*^__version__\s*=\s*['\"](.*?)['\"]\s*$.*", project_init, re.S | re.M
).group(1)
__project__ = re.match(
    r".*^__project__\s*=\s*['\"](.*?)['\"]\s*$.*", project_init, re.S | re.M
).group(1)

setup(
    name=__project__,
    version=__version__,
    author="Adrian Ortiz-Velez",
    authors=["Adrian Ortiz-Velez", "Jeet Sukumaran"],
    author_email="aortizsax@gmail.com",
    packages=find_packages("src"),
    package_dir={"": "src"},
    entry_points={
        "console_scripts": [
            # "name-of-executable = module.with:function_to_execute"
            "geneflowgeometries = geneflowgeometries.application.geneflowgeometries:main",
        ]
    },
    include_package_data=True,
    # MANIFEST.in: only used in source distribution packaging.
    # ``package_data``: only used in binary distribution packaging.
    package_data={
        "": ["*.txt", "*.md", "*.rst",],
        "geneflowgeometries": [
            # For files in this package's direct namespace
            # (e.g., "src/{normalized_project_name}/*.json")
            # "*.json",
            # For files in a (non-subpackage) subdirectory direct namespace
            # (e.g., "src/{normalized_project_name}/resources/config/*.json")
            # "resources/config/*.json",
            # For files located in 'src/geneflowgeometries-data/'
            # "../geneflowgeometries-data/*.json",
            # For files located in 'resources'/'
            # "../../resources/*.json",
        ],
    },
    test_suite="tests",
    # url="http://pypi.python.org/pypi/geneflowgeometries",
    url="https://github.com/aortizsax/geneflowgeometries",
    license="LICENSE",
    description="A Project",
    long_description=_read(["README.md"]),
    long_description_content_type="text/markdown",
    # long_description_content_type="text/x-rst",
    install_requires=_read_requirements("requirements.txt"),
    extras_require={"test": _read_requirements("requirements-test.txt")},
)
