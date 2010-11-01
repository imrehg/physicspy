#!/usr/bin/env python
"""Physicspy: physics calculations in Python

"""

DOCLINES = __doc__.split("\n")

import os
from distutils.core import setup

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def setup_package():
    from distutils.core import setup

    try:
        setup(name='physicspy',
            author = 'Gergely Imreh',
            author_email = 'imrehg@gmail.com',
            version = '0.1.0',
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            url = 'http://github.com/imrehg/physicspy/',
            license = 'MIT',
            platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
            packages = ['physicspy','physicspy.optics','physicspy.quantum'])
    finally:
        pass

    return


if __name__ == '__main__':
    setup_package()

