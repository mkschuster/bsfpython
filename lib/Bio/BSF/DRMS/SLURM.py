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

from Bio.BSF.Database import DatabaseConnection, \
    JobSubmission, JobSubmissionAdaptor, \
    ProcessSLURM, ProcessSLURMAdaptor


output_directory = 'bsfpython_slurm_output'


def submit(drms, debug=0):
    """Submit C{Executable} objects into the Simple Linux Utility for Resource Management (SLURM)
     Distributed Resource Management System (DRMS).

    @param drms: Distributed Resource Management System (C{DRMS})
    @type drms: DRMS
    @param debug: Debug level
    @type debug: int
    """

    # Open or create a database.

    database_path = os.path.join(drms.work_directory, 'bsfpython_slurm_jobs.db')

    database_connection = DatabaseConnection(file_path=database_path)
    database_connection.create_schema()
    # TODO: Not sure it is a good idea to require this after every call?
    # Should this be part of the DatabaseConnection method?
    job_submission_adaptor = JobSubmissionAdaptor(database_connection=database_connection)
    process_slurm_adaptor = ProcessSLURMAdaptor(database_connection=database_connection)

    output = str()
    output += "#! /bin/bash\n"
    output += "\n"

    if debug > 0:
        output += "# BSF-Python debug mode: {}\n".format(debug)
        output += "\n"

    for executable in drms.executables:

        command = list()

        command.append('sbatch')

        # Add DRMS-specific options.

        # Binary or script

        # Job resource string ...

        # SLURM-specific sanity checks ...

        # If a hard memory limit has been set, use it as the minimum free required.

        if drms.memory_limit_hard:
            if not drms.memory_free_virtual:
                drms.memory_free_virtual = drms.memory_limit_hard

        # TODO: The memory must be specified in MB.
        # Maybe it would be worth having a routine that converts suffixes into MB.
        # This should use memory_limit_soft or memory_limit_hard.
        # TODO: Not sure how to use this in a situation like BWA where a multi-threaded application does not use
        # memory for each thread.

        if drms.memory_limit_hard:
            # command.append('--mem-per-cpu')
            command.append('--mem')
            command.append(drms.memory_limit_hard)

        command.append('--time')
        command.append(drms.time_limit)

        # Propagate none of the environment variables.

        command.append('--export')
        command.append('NONE')

        # Get the user environment resembling a login shell.

        command.append('--get-user-env=L')
        # command.append('L')

        # Parallel environment

        if drms.parallel_environment:
            command.append('--distribution')
            command.append(drms.parallel_environment)
            command.append('--ntasks')
            command.append('1')
            command.append('--cpus-per-task')
            command.append(str(drms.threads))

        command.append('--requeue')

        # The --share option may no longer be needed.
        # command.append('--share')

        # Queue name

        if drms.queue:
            command.append('--partition')
            command.append(drms.queue)

        # Working directory, standard output and standard error streams.

        if drms.work_directory:
            command.append('--workdir')
            command.append(drms.work_directory)

            # Write standard output and standard error streams into a
            # 'bsfpython_slurm_output' directory under the 'working_directory'.

            # TODO: Use slurm_output_name to keep --error and --output relative to the --workdir and
            # slurm_output_path to create the directory.
            output_directory_path = os.path.join(drms.work_directory, output_directory)

            if not os.path.isdir(output_directory_path):
                # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.makedirs(output_directory_path)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            command.append('--error')
            command.append(os.path.join(output_directory, string.join(words=(executable.name, '%j.err'), sep='_')))

            command.append('--output')
            command.append(os.path.join(output_directory, string.join(words=(executable.name, '%j.out'), sep='_')))

        # Job name

        if executable.name:
            command.append('--job-name')
            command.append(executable.name)

        # Job hold conditions
        # A particular feature of SLURM is its inability to set process dependencies on process names.
        # Rather, dependencies need setting on the process identifier, which is only obtained after
        # submitting the process. Isn't that exactly what we have a scheduler for? Sigh.
        # Consequently, SLURM process identifiers need to be tracked here, by means of an SQLite database.

        process_identifier_list = list()
        for executable_name in executable.dependencies:
            process_slurm_list = process_slurm_adaptor.select_all_by_job_name(name=executable_name)
            if len(process_slurm_list):
                # This Executable has been submitted at least once before.
                # For the moment, set the dependency on the last submission.
                process_identifier_list.append(process_slurm_list[-1].job_id)
            elif debug == 0:
                warnings.warn(
                    "While submitting Executable with name {!r}, "
                    "Executable with name {!r} that it depends on, "
                    "has not been submitted before.".
                    format(executable.name, executable_name),
                    UserWarning)
        if len(process_identifier_list):
            # Only set the dependency option if there are some on the process identifier list.
            # The identifier list may be empty if no dependencies exist or no Executable has been submitted before.
            command.append('--dependency')
            command.append(string.join(words=map(lambda x: 'afterok:' + x, process_identifier_list), sep=','))

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
            # SLURM sbatch.

            (child_stdout, child_stderr) = child_process.communicate(input=None)

            child_return_code = child_process.returncode

            if child_return_code:
                raise Exception(
                    "SLURM sbatch returned exit code {!r}\n"
                    "STDOUT: {}\n"
                    "STDERR: {}\n"
                    "Command list representation: {!r}".
                    format(child_return_code, child_stdout, child_stderr, command))

            # Parse the multi-line STDOUT string to get the SLURM process identifier and name.
            # The response to the SLURM sbatch command looks like:
            # Submitted batch job 137657
            # Set the result in the Executable.process_identifier instance variable.

            for line in child_stdout.splitlines(False):
                match = re.search(pattern=r'Submitted batch job (\d+)', string=line)
                if match:
                    executable.process_identifier = match.group(1)
                else:
                    print('Could not parse the process identifier from the SLURM sbatch response line {}'.format(line))

        # Copy the SLURM command line to the Bash script.

        output += string.join(words=command, sep=' ') + "\n"
        output += "\n"

        # Regardless of an actual Executable submission, UPDATE it in or INSERT it into the SQLite database.

        job_submission = job_submission_adaptor.select_by_name(name=executable.name)
        if job_submission:
            job_submission.command = executable.command_str()
            job_submission_adaptor.update(data_object=job_submission)
        else:
            job_submission = JobSubmission(
                executable_id=0,
                name=executable.name,
                command=executable.command_str())
            job_submission_adaptor.insert(data_object=job_submission)

        # Only store a ProcessSLURM object, if an Executable has been submitted into SLURM.

        if executable.process_identifier:
            process_slurm = process_slurm_adaptor.select_by_job_id(job_id=executable.process_identifier)
            if not process_slurm:
                process_slurm = ProcessSLURM(job_id=executable.process_identifier, job_name=executable.name)
                process_slurm_adaptor.insert(data_object=process_slurm)

        # The commit statement should affect both insert statements above.
        job_submission_adaptor.database_connection.connection.commit()

    script_path = os.path.join(drms.work_directory, 'bsfpython_slurm_{}.bash'.format(drms.name))
    script_file = open(name=script_path, mode='w')
    script_file.write(output)
    script_file.close()
