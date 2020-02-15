# Git interaction code adapted from https://github.com/desihub/desiutil.
# The desiutil project is distributed under a 3-clause BSD style license:
#
# Copyright (c) 2014-2017, DESI Collaboration <desi-data@desi.lbl.gov>
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of the DESI Collaboration nor the names of its
#   contributors may be used to endorse or promote products derived from this
#   software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
"""Some code for interacting with git.
"""

import re
from os.path import abspath, exists, isdir, isfile, join
from setuptools import Command
from distutils.log import INFO, WARN, ERROR


def get_version():
    """Get the value of ``__version__`` without having to import the module.

    Returns
    -------
    :class:`str`
        The value of ``__version__``.
    """
    ver = 'unknown'
    try:
        version_dir = find_version_directory()
    except IOError:
        return ver
    version_file = join(version_dir, '_version.py')
    if not isfile(version_file):
        update_version()
    with open(version_file, "r") as f:
        for line in f.readlines():
            mo = re.match("__version__ = '(.*)'", line)
            if mo:
                ver = mo.group(1)
    return ver


def get_git_version(git='git'):
    """Use ``git describe`` to generate a version string.

    Parameters
    ----------
    git : :class:`str`, optional
        Path to the git executable, if not in :envvar:`PATH`.

    Returns
    -------
    :class:`str`
        A :pep:`386`-compatible version string.

    Notes
    -----
    The version string should be compatible with :pep:`386` and
    :pep:`440`.
    """
    from subprocess import Popen, PIPE
    myversion = '0.0.1.dev0'
    try:
        p = Popen([git, "describe", "--tags", "--dirty", "--always"],
                  universal_newlines=True, stdout=PIPE, stderr=PIPE)
    except OSError:
        return myversion
    out, err = p.communicate()
    if p.returncode != 0:
        return myversion
    ver = out.rstrip().split('-')[0]+'.dev'
    try:
        p = Popen([git, "rev-list", "--count", "HEAD"],
                  universal_newlines=True, stdout=PIPE, stderr=PIPE)
    except OSError:
        return myversion
    out, err = p.communicate()
    if p.returncode != 0:
        return myversion
    ver += out.rstrip()
    return ver


def update_version(tag=None):
    """Update the _version.py file.

    Parameters
    ----------
    tag : :class:`str`, optional
        Set the version to this string, unconditionally.

    Raises
    ------
    IOError
        If the repository type could not be determined.
    """
    version_dir = find_version_directory()
    if tag is not None:
        ver = tag
    else:
        if isdir(".git"):
            ver = get_git_version()
        else:
            raise IOError("Repository type is not git.")
    version_file = join(version_dir, '_version.py')
    with open(version_file, "w") as f:
        f.write("__version__ = '{}'\n".format(ver))
    return


def find_version_directory():
    """Return the name of a directory containing version information.

    Looks for files in the following places:

    * python/snewpy/_version.py
    * snewpy/_version.py

    Returns
    -------
    :class:`str`
        Name of a directory that can or does contain version information.

    Raises
    ------
    IOError
        If no valid directory can be found.
    """
    packagename='snewpy'
    setup_dir = abspath('.')
    if isdir(join(setup_dir, 'python', packagename)):
        version_dir = join(setup_dir, 'python', packagename)
    elif isdir(join(setup_dir, packagename)):
        version_dir = join(setup_dir, packagename)
    else:
        raise IOError("Could not find a directory containing version information!")
    return version_dir


class SetVersion(Command):
    """Allow users to easily update the package version with
    ``python setup.py version``.
    """
    description = "update _version.py from git repo"
    user_options = [('tag=', 't',
                     'Set the version to a name in preparation for tagging.'),
                    ]
    boolean_options = []

    def initialize_options(self):
        self.tag = None

    def finalize_options(self):
        pass

    def run(self):
        meta = self.distribution.metadata
        update_version(tag=self.tag)
        ver = get_version()
        self.announce("Version is now {}.".format(ver), level=INFO)
