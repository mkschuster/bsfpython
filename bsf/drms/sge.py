# -*- coding: utf-8 -*-
"""Son of Grid Engine (SGE) DRMS module.

A package of methods supporting the Son of Grid Engine (SGE) system as
Distributed Resource Management System (DRMS) module
"""
#  Copyright 2013 - 2021 Michael K. Schuster
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
import errno
import os
import re
from typing import List

from bsf.connector import StandardOutputStream
from bsf.database import DatabaseAdaptor, DatabaseConnection
from bsf.process import Executable

output_directory_name = 'bsfpython_sge_output'


class ProcessSGE(object):
    """C{bsf.drms.sge.ProcessSGE} class representing one Son of Grid Engine (SGE) process.

    The instance variable names result from the SGE accounting file. See man 5 accounting.

    @ivar process_sge_id: Primary key
    @type process_sge_id: int | None
    @ivar qname: Name of the cluster queue in which the job has run
    @type qname: str | None
    @ivar hostname: Name of the execution host
    @type hostname: str | None
    @ivar sge_group: The effective group id of the job owner when executing the job
    @type sge_group: str | None
    @ivar owner: Owner of the Grid Engine job
    @type owner: str | None
    @ivar job_name: Job name
    @type job_name: str | None
    @ivar job_number: Job identifier (job number)
    @type job_number: str | None
    @ivar account: An account string as specified by the qsub(1) or qalter(1) -A option
    @type account: str | None
    @ivar priority: Priority value assigned to the job, corresponding to the priority parameter in the
        queue configuration (see queue_conf(5))
    @type priority: str | None
    @ivar submission_time: Submission time
    @type submission_time: str | None
    @ivar start_date: Start time
    @type start_date: str | None
    @ivar end_time: End time
    @type end_time: str | None
    @ivar failed: Indicates the problem which occurred in case a job could not be started on the execution host
    @type failed: str | None
    @ivar exit_status: Exit status of the job script (or Grid Engine-specific status in case of certain error
        conditions). The exit status is determined by following the normal shell conventions. If the command
        terminates normally, the value of the command is its exit status. However, in the case that the command
        exits abnormally, a value of 0200 (octal), 128 (decimal) is added to the value of the command to make up
        the exit status.
    @type exit_status: str | None
    @ivar ru_wallclock: Difference between end_time and start_time (see above), except that if the job fails,
        it is zero
    @type ru_wallclock: str | None
    @ivar project: The department which was assigned to the job
    @type project: str | None
    @ivar department: The parallel environment which was selected for the job
    @type department: str | None
    @ivar granted_pe: The number of slots which were dispatched to the job by the scheduler
    @type granted_pe: str | None
    @ivar slots: The number of slots which were dispatched to the job by the scheduler
    @type slots: str | None
    @ivar task_number: Array job task index number
    @type task_number: str | None
    @ivar cpu: The CPU time usage in seconds
    @type cpu: str | None
    @ivar mem: The integral memory usage in Gbytes seconds
    @type mem: str | None
    @ivar io: The amount of data transferred in input/output operations in GB (if available, otherwise 0)
    @type io: str | None
    @ivar category: A string specifying the job category
    @type category: str | None
    @ivar iow: The input/output wait time in seconds (if available, otherwise 0)
    @type iow: str | None
    @ivar pe_taskid: If this identifier is set, the task was part of a parallel job, and was passed to Grid Engine
        via the qrsh -inherit interface.
    @type pe_taskid: str | None
    @ivar maxvmem: The maximum vmem size in bytes
    @type maxvmem: str | None
    @ivar arid: Advance reservation identifier
    @type arid: str | None
    """

    def __init__(
            self,
            process_sge_id=None,
            qname=None,
            hostname=None,
            sge_group=None,
            owner=None,
            job_name=None,
            job_number=None,
            account=None,
            priority=None,
            submission_time=None,
            start_date=None,
            end_time=None,
            failed=None,
            exit_status=None,
            ru_wallclock=None,
            project=None,
            department=None,
            granted_pe=None,
            slots=None,
            task_number=None,
            cpu=None,
            mem=None,
            io=None,
            category=None,
            iow=None,
            pe_taskid=None,
            maxvmem=None,
            arid=None):
        """Initialise a C{bsf.drms.sge.ProcessSGE}.

        @param process_sge_id: Primary key
        @type process_sge_id: int | None
        @param qname: Name of the cluster queue in which the job has run
        @type qname: str | None
        @param hostname: Name of the execution host
        @type hostname: str | None
        @param sge_group: The effective group id of the job owner when executing the job
        @type sge_group: str | None
        @param owner: Owner of the Grid Engine job
        @type owner: str | None
        @param job_name: Job name
        @type job_name: str | None
        @param job_number: Job identifier (job number)
        @type job_number: str | None
        @param account: An account string as specified by the qsub(1) or qalter(1) -A option
        @type account: str | None
        @param priority: Priority value assigned to the job, corresponding to the priority parameter in the
            queue configuration (see queue_conf(5))
        @type priority: str | None
        @param submission_time: Submission time
        @type submission_time: str | None
        @param start_date: Start time
        @type start_date: str | None
        @param end_time: End time
        @type end_time: str | None
        @param failed: Indicates the problem which occurred in case a job could not be started on the execution host
        @type failed: str | None
        @param exit_status: Exit status of the job script (or Grid Engine-specific status in case of certain error
            conditions). The exit status is determined by following the normal shell conventions. If the command
            terminates normally, the value of the command is its exit status. However, in the case that the command
            exits abnormally, a value of 0200 (octal), 128 (decimal) is added to the value of the command to make up
            the exit status.
        @type exit_status: str | None
        @param ru_wallclock: Difference between end_time and start_time (see above), except that if the job fails,
            it is zero
        @type ru_wallclock: str | None
        @param project: The department which was assigned to the job
        @type project: str | None
        @param department: The parallel environment which was selected for the job
        @type department: str | None
        @param granted_pe: The number of slots which were dispatched to the job by the scheduler
        @type granted_pe: str | None
        @param slots: The number of slots which were dispatched to the job by the scheduler
        @type slots: str | None
        @param task_number: Array job task index number
        @type task_number: str | None
        @param cpu: The CPU time usage in seconds
        @type cpu: str | None
        @param mem: The integral memory usage in Gbytes seconds
        @type mem: str | None
        @param io: The amount of data transferred in input/output operations in GB (if available, otherwise 0)
        @type io: str | None
        @param category: A string specifying the job category
        @type category: str | None
        @param iow: The input/output wait time in seconds (if available, otherwise 0)
        @type iow: str | None
        @param pe_taskid: If this identifier is set, the task was part of a parallel job, and was passed to Grid Engine
            via the qrsh -inherit interface.
        @type pe_taskid: str | None
        @param maxvmem: The maximum vmem size in bytes
        @type maxvmem: str | None
        @param arid: Advance reservation identifier
        @type arid: str | None
        """
        super(ProcessSGE, self).__init__()

        self.process_sge_id = process_sge_id
        self.qname = qname
        self.hostname = hostname
        self.sge_group = sge_group
        self.owner = owner
        self.job_name = job_name
        self.job_number = job_number
        self.account = account
        self.priority = priority
        self.submission_time = submission_time
        self.start_date = start_date
        self.end_time = end_time
        self.failed = failed
        self.exit_status = exit_status
        self.ru_wallclock = ru_wallclock
        self.project = project
        self.department = department
        self.granted_pe = granted_pe
        self.slots = slots
        self.task_number = task_number
        self.cpu = cpu
        self.mem = mem
        self.io = io
        self.category = category
        self.iow = iow
        self.pe_taskid = pe_taskid
        self.maxvmem = maxvmem
        self.arid = arid

        return


class ProcessSGEAdaptor(DatabaseAdaptor):
    """C{bsf.drms.sge.ProcessSGEAdaptor} class providing database access for the C{bsf.drms.sge.ProcessSGE} class.

    The SQL column names result from the SGE accounting file. See man 5 accounting.
    """

    def __init__(self, database_connection):
        """Initialise a C{bsf.drms.sge.ProcessSGEAdaptor}.

        @param database_connection: C{bsf.database.DatabaseConnection}
        @type database_connection: DatabaseConnection
        """
        super(ProcessSGEAdaptor, self).__init__(
            database_connection=database_connection,
            object_type=ProcessSGE,
            table_name='process_sge',
            column_definition=[
                # Primary key
                ('process_sge_id', 'INTEGER PRIMARY KEY ASC AUTOINCREMENT'),
                # qname
                # Name of the cluster queue in which the job has run.
                ('qname', 'TEXT'),
                # hostname
                # Name of the execution host.
                ('hostname', 'TEXT'),
                # group
                # The effective group id of the job owner when executing the job.
                # Since group is a reserved word in SQL this had to be renamed to sge_group.
                ('sge_group', 'TEXT'),
                # owner
                # The owner of the Grid Engine job.
                ('owner', 'TEXT'),
                # job_name
                # Job name.
                ('job_name', 'TEXT'),
                # job_number
                # Job identifier (job number).
                ('job_number', 'TEXT'),
                # account
                # An account string as specified by the qsub(1) or qalter(1) -A option.
                ('account', 'TEXT'),
                # priority
                # The priority value assigned to the job,
                # corresponding to the priority parameter in the queue configuration (see queue_conf(5)).
                ('priority', 'TEXT'),
                # submission_time
                # Submission time.
                ('submission_time', 'TEXT'),
                # start_time
                # Start time.
                ('start_date', 'TEXT'),
                # end_time
                # End time.
                ('end_time', 'TEXT'),
                # failed
                # Indicates the problem which occurred in case a job could not be started on the execution host
                # (e.g. because the owner of the job did not have a valid account on that machine).
                # If Grid Engine tries to start a job multiple times, this may lead to multiple entries in the
                # reporting file corresponding to the same job ID.
                ('failed', 'TEXT'),
                # exit_status
                # Exit status of the job script (or Grid Engine-specific status in case of certain error conditions).
                # The exit status is determined by following the normal shell conventions. If the command terminates
                # normally, the value of the command is its exit status. However, in the case that the command exits
                # abnormally, a value of 0200 (octal), 128 (decimal) is added to the value of the command to make up
                # the exit status.
                #
                # For example: If a job dies through signal 9 (SIGKILL) then the exit status becomes 128 + 9 = 137.
                ('exit_status', 'TEXT'),
                # ru_wallclock
                # Difference between end_time and start_time (see above), except that if the job fails, it is zero.
                ('ru_wallclock', 'TEXT'),
                # ru_utime
                # ru_stime
                # ru_maxrss
                # ru_ixrss
                # ru_ismrss
                # ru_idrss
                # ru_isrss
                # ru_minflt
                # ru_majflt
                # ru_nswap
                # ru_inblock
                # ru_oublock
                # ru_msgsnd
                # ru_msgrcv
                # ru_nsignals
                # ru_nvcsw
                # ru_nivcsw
                # These entries follow the contents of the standard Unix rusage structure as described in getrusage(2).
                # Depending on the operating system where the job was executed, some fields may be 0.
                #
                # project
                # The project which was assigned to the job.
                ('project', 'TEXT'),
                # department
                # The department which was assigned to the job.
                ('department', 'TEXT'),
                # granted_pe
                # The parallel environment which was selected for the job.
                ('granted_pe', 'TEXT'),
                # slots
                # The number of slots which were dispatched to the job by the scheduler.
                ('slots', 'TEXT'),
                # task_number
                # Array job task index number.
                ('task_number', 'TEXT'),
                # cpu
                # The CPU time usage in seconds.
                # The value may be affected by the ACCT_RESERVED_USAGE execd parameter (see sge_conf(5)).
                ('cpu', 'TEXT'),
                # mem
                # The integral memory usage in Gbytes seconds.
                # The value may be affected by the ACCT_RESERVED_USAGE execd parameter (see sge_conf(5)).
                ('mem', 'TEXT'),
                # io
                # The amount of data transferred in input/output operations (if available, otherwise 0).
                ('io', 'TEXT'),
                # category
                # A string specifying the job category.
                # This contains a space-separated pseudo options list for the sub, with components as follows:
                #
                #   -U user_list
                #       An owner/group ACL list composed of host_conf(5), sge_pe(5),
                #       And queue_conf(5) user_lists/xuser_lists entries.
                #       Entries from sge_conf(5) are not considered since they can
                #       only cause a job to be accepted/rejected at submit time.
                #       Omitted if there are no such configuration entries.
                #
                #   -P project_list
                #       Like -U, but for project/xproject entries.
                #
                #   -u owner
                #       The owner's username, if it was referenced in any RQS (see sge_resource_quota(5)).
                #       Omitted if there was no such reference.
                #
                #   -q queue_list
                #       The hard queue list (only if one was specified).
                #
                #   -masterq queue_list
                #       The master queue list (only if one was specified).
                #
                #   -l resource_list
                #       The hard resource list (only if hard resources were specified).
                #
                #   -soft -l resource_list
                #       The soft resource list (only if soft resources were specified).
                #
                #   -pe pe_name pe_range
                #       The parallel environment specified for the job (only for parallel jobs).
                #
                #   -ckpt ckpt_name
                #   The job's checkpointing environment (only if one was specified).
                #
                #   -I y
                #       Present only for interactive jobs.
                #
                #   -ar ar_id
                #       The advance reservation into which the job was submitted (only if one was specified).
                ('category', 'TEXT'),
                # iow
                # The input/output wait time in seconds (if available, otherwise 0).
                ('iow', 'TEXT'),
                # pe_taskid
                # If this identifier is set, the task was part of a parallel job, and was passed to Grid Engine
                # via the qrsh -inherit interface.
                ('pe_taskid', 'TEXT'),
                # maxvmem
                # The maximum vmem size in bytes.
                # The value may be affected by the ACCT_RESERVED_USAGE execd parameter (see sge_conf(5)).
                ('maxvmem', 'TEXT'),
                # arid
                # Advance reservation identifier. If the job used the resources of an advance reservation,
                # then this field contains a positive integer identifier; otherwise the value is '0' .
                ('arid', 'TEXT'),
            ])

        return


def submit(stage, debug=0):
    """Submit each C{bsf.process.Executable} of a C{bsf.analysis.Stage}.

    Submits each C{bsf.process.Executable} into the
    Son of Grid Engine (SGE)
    Distributed Resource Management System (DRMS).

    @param stage: C{bsf.analysis.Stage}
    @type stage: bsf.analysis.Stage
    @param debug: Debug level
    @type debug: int
    """

    def submit_qsub_stdout(_file_handle, _thread_lock, _debug, _executable):
        """Thread callable to process the SGE qsub STDOUT stream.

        Parses the process identifier returned by SGE qsub and sets it as
        C{bsf.process.Executable.process_identifier}.
        The response to the SGE qsub command looks like:
        Your job 137657 ("ls") has been submitted
        @param _file_handle: File handle (i.e. pipe)
        @type _file_handle: io.TextIOWrapper
        @param _thread_lock: Thread lock
        @type _thread_lock: threading.Lock
        @param _debug: Debug level
        @type _debug: int
        @param _executable: C{bsf.process.Executable}
        @type _executable: Executable
        """
        for _line in _file_handle:
            if _debug > 0:
                _thread_lock.acquire(True)
                print('Line:', _line)
                _thread_lock.release()

            _match = re.search(pattern=r'Your job (\d+) \("([^"]+)"\) has been submitted', string=_line)

            if _match:
                _executable.process_identifier = _match.group(1)
                # _executable.process_name = _match.group(2)
            else:
                _thread_lock.acquire(True)
                print('Could not parse SGE qsub response line:', repr(_line))
                _thread_lock.release()

        return

    output_list: List[str] = list()

    output_list.append('#!/usr/bin/env bash\n')
    output_list.append('\n')

    if debug > 0:
        output_list.append('# BSF-Python debug mode: ' + repr(debug) + '\n')
        output_list.append('\n')

    for executable in stage.executable_list:
        executable_drms = Executable(
            name=executable.name,
            program='qsub',
            sub_command=executable,
            stdout=StandardOutputStream(
                thread_callable=submit_qsub_stdout,
                thread_kwargs={'_executable': executable}))

        # Add Stage-specific options.

        # Clear settings inherited from $SGE_ROOT/$SGE_CELL/common/sge_request,
        # which currently specifies the current working directory (-cwd), which is
        # explicitly set by this code base and all environment variables (-V), which
        # should be defined by the Bash startup files.
        # In this case the shell (-S) needs specifying explicitly.

        executable_drms.add_switch_short(key='clear')
        executable_drms.add_option_short(key='S', value='/bin/bash')

        # Binary or script

        if stage.is_script:
            executable_drms.add_option_short(key='b', value='no')
        else:
            executable_drms.add_option_short(key='b', value='yes')

        # Job resource string ...

        # SGE-specific sanity checks ...

        # If a hard memory limit has been set, use it as the minimum free required.

        if stage.memory_limit_hard and not stage.memory_free_virtual:
            stage.memory_free_virtual = stage.memory_limit_hard

        resource_list: List[str] = list()

        # Require physical memory to be free ...
        if stage.memory_free_mem:
            resource_list.append('mem_free=' + stage.memory_free_mem)

        # Require swap memory to be free ...
        if stage.memory_free_swap:
            resource_list.append('swap_free=' + stage.memory_free_swap)

        # Require virtual memory to be free ...
        if stage.memory_free_virtual:
            resource_list.append('virtual_free=' + stage.memory_free_virtual)

        # Set hard virtual memory limit ...
        if stage.memory_limit_hard:
            resource_list.append('h_vmem=' + stage.memory_limit_hard)

        # Set soft virtual memory limit ...
        if stage.memory_limit_soft:
            resource_list.append('s_vmem=' + stage.memory_limit_soft)

        # Set the host names ...
        if stage.node_list_include:
            for node_name in stage.node_list_include:
                resource_list.append('hostname=' + node_name)

        if len(resource_list):
            executable_drms.add_option_short(key='l', value=','.join(resource_list))

        # Parallel environment

        if stage.parallel_environment:
            # Parallel environment format: -pe pe_name pe_min-pe_max
            # Here, pe_max is not specified, but defaults to 9999999. See qsub (1).
            executable_drms.add_option_multi_short(
                key='pe',
                value=' '.join((stage.parallel_environment, str(stage.threads))))

        # Queue name

        if stage.queue:
            executable_drms.add_option_short(key='q', value=stage.queue)

        # Working directory, standard output and standard error streams.

        if stage.working_directory:
            executable_drms.add_option_short(key='wd', value=stage.working_directory)

            # Write standard output and standard error streams into an
            # output directory under the working directory.
            # Create the output directory first via its absolute path,
            # then set standard error and standard output relative to it.

            output_directory_path = os.path.join(stage.working_directory, output_directory_name)

            if not os.path.isdir(output_directory_path):
                # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.makedirs(output_directory_path)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            executable_drms.add_option_short(key='e', value=output_directory_name)
            executable_drms.add_option_short(key='o', value=output_directory_name)

        # Add bsf.process.Executable-specific options.

        if executable.hold:
            executable_drms.add_switch_short(key='h')
            # The SGE qsub command can use -h to place a user hold.
            # The second form -h {u|s|o|n|U|O|S}... is only for the SGE qalter command.

        # Job name

        if executable.name:
            executable_drms.add_option_short(key='N', value=executable.name)

        # Job hold conditions

        if len(executable.dependencies):
            executable_drms.add_option_short(key='hold_jid', value=','.join(executable.dependencies))

        # Finally, submit this command if requested and not in debug mode.

        if executable.submit and debug == 0:
            exception_str_list = executable_drms.run(debug=debug)

            if exception_str_list:
                exception_str_list.append('Command list representation: ' + repr(executable_drms.command_list()))
                raise Exception('\n'.join(exception_str_list))

        # Copy the SGE command line to the Bash script.

        output_list.append(executable_drms.command_str() + '\n')
        output_list.append('\n')

    script_path = os.path.join(stage.working_directory, 'bsfpython_sge_' + stage.name + '.bash')
    with open(file=script_path, mode='wt') as script_file:
        script_file.writelines(output_list)

    return


def check_state(stage, debug=0):
    """Check the state of each C{bsf.process.Executable} of a C{bsf.analysis.Stage}.

    @param stage: C{bsf.analysis.Stage}
    @type stage: bsf.analysis.Stage
    @param debug: Debug level
    @type debug: int
    """
    if stage:
        pass

    if debug:
        pass

    return
