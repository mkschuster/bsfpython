"""bsf.drms.sge

A package of methods supporting the Son of Grid Engine (SGE) system.
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


import errno
import os.path
import re
import subprocess

# TODO: This module could create a file that records SGE Process identifiers, which could be used by
# further scripts to query the state of jobs.


output_directory = 'bsfpython_sge_output'


def submit(drms, debug=0):

    """Submit C{Executable} objects into the Son of Grid Engine (SGE) Distributed Resource Management System (DRMS).

    @param drms: Distributed Resource Management System (C{DRMS})
    @type drms: DRMS
    @param debug: Debug level
    @type debug: int
    """

    output = str()
    output += "#! /bin/bash\n"
    output += "\n"

    if debug > 0:
        output += "# BSF-Python debug mode: {}\n".format(debug)
        output += "\n"

    for executable in drms.executables:

        command = list()

        command.append('qsub')

        # Add DRMS-specific options.

        # Clear settings inherited from $SGE_ROOT/$SGE_CELL/common/sge_request,
        # which currently specifies the current working directory (-cwd), which is
        # explicitly set by this code base and all environment variables (-V), which
        # should be defined by the Bash startup files.
        # In this case the shell (-S) needs specifying explicitly.

        # TODO: Clearing the environment seems linked to Python (numpy?) libimf.so failures.
        # TODO: Test this again for the new Python installation at some stage ...
        # command.append('-clear')
        # command.append('-S')
        # command.append('/bin/bash')

        # Binary or script

        command.append('-b')
        if drms.is_script:
            command.append('no')
        else:
            command.append('yes')

        # Job resource string ...

        # SGE-specific sanity checks ...

        # If a hard memory limit has been set, use it as the minimum free required.

        if drms.memory_limit_hard:
            if not drms.memory_free_virtual:
                drms.memory_free_virtual = drms.memory_limit_hard

        resource_list = list()

        # Require physical memory to be free ...
        if drms.memory_free_mem:
            resource_list.append('mem_free={}'.format(drms.memory_free_mem))

        # Require swap memory to be free ...
        if drms.memory_free_swap:
            resource_list.append('swap_free={}'.format(drms.memory_free_swap))

        # Require virtual memory to be free ...
        if drms.memory_free_virtual:
            resource_list.append('virtual_free={}'.format(drms.memory_free_virtual))

        # Set hard virtual memory limit ...
        if drms.memory_limit_hard:
            resource_list.append('h_vmem={}'.format(drms.memory_limit_hard))

        # Set soft virtual memory limit ...
        if drms.memory_limit_soft:
            resource_list.append('s_vmem={}'.format(drms.memory_limit_soft))

        if len(resource_list):

            command.append('-l')
            command.append(','.join(resource_list))

        # Parallel environment

        if drms.parallel_environment:
            command.append('-pe')
            command.append(drms.parallel_environment)
            command.append(str(drms.threads))

        # Queue name

        if drms.queue:
            command.append('-q')
            command.append(drms.queue)

        # Working directory, standard output and standard error streams.

        if drms.working_directory:
            command.append('-wd')
            command.append(drms.working_directory)

            # Write standard output and standard error streams into a
            # 'bsfpython_sge_output' directory under the 'working_directory'.

            output_directory_path = os.path.join(drms.working_directory, output_directory)

            if not os.path.isdir(output_directory_path):
                # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.makedirs(output_directory_path)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            command.append('-e')
            command.append(output_directory)

            command.append('-o')
            command.append(output_directory)

        # Add Executable-specific options.

        if executable.hold:
            command.append('-h')
            # The SGE qsub command can use -h to place a user hold.
            # The second form -h {u|s|o|n|U|O|S}... is only for the SGE qalter command.

        # Job name

        if executable.name:
            command.append('-N')
            command.append(executable.name)

        # Job hold conditions

        if len(executable.dependencies):
            command.append('-hold_jid')
            command.append(','.join(executable.dependencies))

        command.extend(executable.command_list())

        if executable.stdout_path:
            command.append("1>{}".format(executable.stdout_path))
        if executable.stderr_path:
            command.append("2>{}".format(executable.stderr_path))

        # Finally, submit this command if requested and not in debug mode.

        if executable.submit and debug == 0:

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
                raise Exception(
                    "SGE qsub returned exit code {!r}\n"
                    "STDOUT: {}\n"
                    "STDERR: {}\n"
                    "Command list representation: {!r}".
                    format(child_return_code, child_stdout, child_stderr, command))

            # Parse the multi-line STDOUT string to get the SGE process identifier and name.
            # The response to the SGE qsub command looks like:
            # Your job 137657 ("ls") has been submitted

            for line in child_stdout.splitlines(False):
                match = re.search(pattern=r'Your job (\d+) \("([^"]+)"\) has been submitted',
                                  string=line)

                if match:
                    executable.process_identifier = match.group(1)
                    executable.process_name = match.group(2)
                else:
                    print('Could not parse SGE qsub response line {}'.format(line))

        # Copy the SGE command line to the Bash script.

        output += ' '.join(command) + "\n"
        output += "\n"

    script_path = os.path.join(drms.working_directory, 'bsfpython_sge_{}.bash'.format(drms.name))
    script_file = open(name=script_path, mode='w')
    script_file.write(output)
    script_file.close()
