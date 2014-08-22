"""Bio

The Bio package aligns the BSF Python project to the
Open Bioinformatics Foundation (OBF) BioPython project.

http://www.biopython.org/
http://www.open-bio.org/
"""

#
# Copyright 2013 Michael K. Schuster
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

from pkgutil import extend_path

# For compatibility with the BioPython project, extend the path.
# Importantly, in the PYTHONPATH environment variable, the bsfpython module
# needs to be specified after BioPython.

__path__ = extend_path(__path__, __name__)
