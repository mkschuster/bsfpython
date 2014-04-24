"""Bio.BSF.DRMS.SLURM

A package of methods supporting the Simple Linux Utility for Resource Management (SLURM) system.
"""

#
# Copyright 2014 Michael K. Schuster
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


import errno
import os.path
import re
import string
import subprocess
import warnings

# Global dictionary to map from Executable.name entries to SLURM job identifiers,
# which are required for specifying job dependencies. Isn't that what a DRMS should do for us?

executable_name_dict = dict()

# TODO: This module could create a file that records SLURM Process identifiers, which could be used by
# further scripts to query the state of jobs.
# Since SLURM only reports job identifiers upon submission via sbatch, it is no longer possible to submit
# analysis jobs in stages. The only way to record job identifiers of previous analyses is via a file. Sigh.


def submit(self, debug=0):

    """Submit BSF Executable objects into the Simple Linux Utility for Resource Management (SLURM)
     Distributed Resource Management System.

    :param self: BSF DRMS
    :type self: DRMS
    :param debug: Debug level
    :type debug: int
    :return: Nothing
    :rtype: None
    """

    output = str()
    output += "#! /bin/bash\n"
    output += "\n"

    if debug > 0:
        output += "# BSF-Python debug mode: {}\n".format(debug)
        output += "\n"

    for executable in self.executables:

        command = list()

        command.append('sbatch')

        # Add DRMS-specific options.

        # Binary or script

        # Job resource string ...

        # SLURM-specific sanity checks ...

        # If a hard memory limit has been set, use it as the minimum free required.

        if self.memory_limit_hard:
            if not self.memory_free_virtual:
                self.memory_free_virtual = self.memory_limit_hard

        # TODO: The memory must be specified in MB.
        # Maybe it would be worth having a routine that converts suffixes into MB.
        # This should use memory_limit_soft or memory_limit_hard.

        if self.memory_limit_hard:
            command.append('--mem')
            command.append(self.memory_limit_hard)

        # Parallel environment

        if self.parallel_environment:
            command.append('--distribution')
            command.append(self.parallel_environment)
            command.append('--ntasks')
            command.append('1')
            command.append('--cpus-per-task')
            command.append(str(self.threads))
            command.append('--share')

        # Queue name

        if self.queue:
            command.append('--partition')
            command.append(self.queue)

        # Working directory, standard output and standard error streams.

        if self.work_directory:
            command.append('--workdir')
            command.append(self.work_directory)

            # Write standard output and standard error streams into a
            # 'bsfpython_slurm_output' directory under the 'working_directory'.

            # TODO: Use slurm_output_name to keep --error and --output relative to the --workdir and
            # slurm_output_path to create the directory.
            slurm_output_directory = os.path.join(self.work_directory, 'bsfpython_slurm_output')

            if not os.path.isdir(slurm_output_directory):
                # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.makedirs(slurm_output_directory)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            command.append('--error')
            command.append(os.path.join(slurm_output_directory, '{}_%j.err'.format(executable.name)))

            command.append('--output')
            command.append(os.path.join(slurm_output_directory, '{}_%j.out'.format(executable.name)))

        # Job name

        if executable.name:
            command.append('--job-name')
            command.append(executable.name)

        # Job hold conditions

        if len(executable.dependencies):
            identifier_list = list()
            for executable_name in executable.dependencies:
                if executable_name in executable_name_dict:
                    identifier_list.append(executable_name_dict[executable_name])
                elif debug == 0:
                    # Dependencies can only be calculated if jobs have been submitted.
                    # TODO: This means, it is no longer possible to submit jobs in analysis stages.
                    # It is thus necessary to write sample annotation sheets that contain SLURM job identifiers.
                    message = "While submitting Executable with name {!r}, " \
                              "Executable with name {!r} that it depends on, has not been submitted before.". \
                        format(executable.name, executable_name)
                    warnings.warn(message, UserWarning)
            if len(identifier_list):
                # If no jobs have been submitted before, the identifier list is empty.
                command.append('--dependency')
                command.append(string.join(map(lambda x: 'afterok:' + x, identifier_list), ','))

        command.extend(executable.command_list())

        if executable.stdout_path:
            command.append("1>{}".format(executable.stdout_path))
        if executable.stderr_path:
            command.append("2>{}".format(executable.stderr_path))

        # Finally, submit this command if not in debug mode.

        if debug == 0:

            child_process = subprocess.Popen(args=command,
                                             bufsize=4096,
                                             stdin=subprocess.PIPE,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE,
                                             shell=False,
                                             close_fds=True)

            # Although subprocess.communicate() may block when memory buffers
            # have been filled up, not much STDOUT and STDERR is expected from
            # SGE qsub.

            (child_stdout, child_stderr) = child_process.communicate(input=None)

            child_return_code = child_process.returncode

            if child_return_code:
                message = "SLURM sbatch exit code {!r}\n".format(child_return_code)
                message += "STDOUT: {}\n".format(child_stdout)
                message += "STDERR: {}\n".format(child_stderr)
                message += "Command list representation: {!r}".format(command)
                raise Exception(message)

            # Parse the multi-line STDOUT string to get the SGE process identifier and name.
            # The response to the SGE qsub command looks like:
            # Your job 137657 ("ls") has been submitted

            for line in child_stdout.splitlines(False):
                match = re.search(pattern=r'Submitted batch job (\d+)', string=line)
                if match:
                    executable.process_identifier = match.group(1)
                    # Correlate Executable.name and Executable.process_identifier information.
                    if executable.name in executable_name_dict:
                        message = "Overwriting Executable with name {!r} and process identifier {!r} " \
                                  "that has been submitted to SLURM before.". \
                            format(executable.name, executable.process_identifier)
                        warnings.warn(message, UserWarning)
                    executable_name_dict[executable.name] = executable.process_identifier
                else:
                    print('Could not parse SLURM sbatch response line {}'.format(line))

        # Copy the SLURM command line to the Bash script.

        output += string.join(words=command, sep=' ') + "\n"
        output += "\n"

    script_path = os.path.join(self.work_directory, 'bsfpython_slurm_{}.sh'.format(self.name))
    script_file = open(name=script_path, mode='w')
    script_file.write(output)
    script_file.close()
