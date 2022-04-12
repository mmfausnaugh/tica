#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__package__ = "tica"
__author__ = 'Michael Fausnaugh'

import logging
import platform
import locale
import sys
import os.path

fstem = os.path.abspath(os.path.dirname(__file__) + '/../')
with open(os.path.join(fstem, 'VERSION'),'r') as infile:
    version = infile.read()

def platform_info():
    lines = []
    lines.append('TICA version: {}'.format(version) )
    lines.append('Python version: {}'.format( sys.version.replace('\n', ' ') ))
    lines.append('Host: {}'.format( platform.node() ) )
    lines.append('Platform: {}'.format( platform.platform()) )
    lines.append('Locale: {}'.format( locale.setlocale(locale.LC_ALL)) )
    return lines

def setup_logging(debug=False, filename=None):
    """Provide a sane set of defaults for logging."""
    level = logging.DEBUG if debug else logging.INFO
    if filename is None:
        filename = '-'
    if filename == '-':
        hand = logging.StreamHandler()
    else:
        hand = logging.FileHandler(filename)
    root_logger = logging.getLogger()

    fmt = '%(asctime)s %(levelname)s %(funcName)s: %(message)s' if level == logging.DEBUG else '%(asctime)s %(message)s'
    datefmt = '%Y.%m.%d %H:%M:%S'
    hand.setFormatter(logging.Formatter(fmt, datefmt))

    root_logger.setLevel(level)
    root_logger.handlers = []
    root_logger.addHandler(hand)

    for line in platform_info():
        logging.debug(line)
