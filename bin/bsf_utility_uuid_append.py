#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2013 - 2022 Michael K. Schuster
#
#  Biomedical Sequencing Facility (BSF), part of the genomics core facility
#  of the Research Center for Molecular Medicine (CeMM) of the
#  Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
#  This file is part of BSF Python.
#
#  BSF Python is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BSF Python is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.
#
"""The :py:mod:`bin.bsf_utility_uuid_append` module is a script to append a :emphasis:`UUID` to a string.
"""

import sys
import uuid
from argparse import ArgumentParser
from typing import Optional


def run(input_value: Optional[str] = None) -> int:
    """Run function.

    :param input_value: An input value.
    :type input_value: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if input_value:
        print('_'.join((input_value, uuid.uuid4().hex)))
    else:
        print(uuid.uuid4().hex)

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Append a UUID to a string.')

    argument_parser.add_argument(
        'input_value',
        help="input string, can be empty by supplying ''")

    name_space = argument_parser.parse_args()

    return run(input_value=name_space.input_value)


if __name__ == '__main__':
    sys.exit(main())
