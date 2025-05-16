#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

requirements = [
    "cftime",
    "fv3config",
    "xarray",
]

setup_requirements = []

test_requirements = ["pytest"]

setup(
    author="Allen Institute of Artificial Intelligence",
    author_email="spencerc@allenai.org",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
    description="shield_run is used when running simulations on Gaea.",
    install_requires=requirements,
    license="Apache 2.0 license",
    include_package_data=True,
    keywords="shield_run",
    name="shield_run",
    packages=find_packages(include=["shield_run", "shield_run.*"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    version="0.1.0",
    zip_safe=False,
)
