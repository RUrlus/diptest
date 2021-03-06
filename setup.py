"""Project: Diptest

Created: 2021/09/24

Description:
    setup script to install diptest package.

Authors:
    Ralph Urlus [rurlus.dev@gmail.com]

Redistribution and use in source and binary forms, with or without
modification, are permitted according to the terms listed in the file
LICENSE.
"""

import os
import pybind11
from setuptools import find_packages
from skbuild import setup

NAME = 'diptest'

MAJOR = 0
REVISION = 5
PATCH = 1
DEV = False

VERSION = '{major}.{revision}.{patch}'.format(major=MAJOR, revision=REVISION, patch=PATCH)
FULL_VERSION = VERSION
if DEV:
    FULL_VERSION += '.dev'

# read the contents of readme file
with open("README.md", encoding="utf-8") as f:
    long_description = f.read()


def write_version_py(filename: str = 'diptest/version.py') -> None:
    """Write package version to version.py.

    This will ensure that the version in version.py is in sync with us.

    Parameters
    ----------
    filename : str
        the path the file to write the version.py

    """
    # Do not modify the indentation of version_str!
    version_str = """\"\"\"THIS FILE IS AUTO-GENERATED BY diptest SETUP.PY.\"\"\"

name = '{name!s}'
version = '{version!s}'
full_version = '{full_version!s}'
release = {is_release!s}
"""

    with open(filename, 'w') as version_file:
        version_file.write(
            version_str.format(name=NAME.lower(), version=VERSION, full_version=FULL_VERSION, is_release=not DEV)
        )


if __name__ == '__main__':
    write_version_py()
    if 'DIPTEST_MANUAL_BUILD' in os.environ:
        from setuptools import setup
        print('Diptest: running pip install without extension')

    setup(
        name=NAME,
        packages=find_packages(),
        version=FULL_VERSION,
        cmake_args=[
            f"-DDIPTEST_VERSION_INFO:STRING={VERSION}",
            f"-Dpybind11_DIR:STRING={pybind11.get_cmake_dir()}",
            "-DDIPTEST_ENABLE_ARCH_FLAGS:BOOL=ON",
        ]
    )
