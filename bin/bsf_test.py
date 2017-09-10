#! /usr/bin/env python
#
# BSF Python script to test library functions.
#
#
# Copyright 2013 - 2016 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility
# of the Research Center for Molecular Medicine (CeMM) of the
# Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
# This file is part of BSF Python.
#
# BSF Python is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BSF Python is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.


from bsf.argument import Argument

for key, value in (
        ('-switch_short', None),
        ('--switch_long', None),
        ('-option', 'short'),
        ('--option', 'long'),
        ('option', 'pair'),
        ('option', 'pair space'),
        ('-option_pair=short', None),
        ('--option_pair=long', None),
        ('--option=pair space', None)):
    test_argument = Argument.from_key_value(key=key, value=value)

    print "Argument: {!r}".format(test_argument)
    print "Argument.get_str(): {!r}".format(test_argument.get_str())
    print "Argument.get_list(): {!r}".format(test_argument.get_list())
