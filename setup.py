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
            version='0.0.1',
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            license='MIT',
            platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
            packages=['physicspy','physicspy.optics'])
    finally:
        pass

    return


if __name__ == '__main__':
    setup_package()

