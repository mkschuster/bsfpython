# -*- coding: utf-8 -*-
"""Simple Linux Utility for Resource Management (SLURM) DRMS module

A package of methods supporting the Simple Linux Utility for Resource Management (SLURM) system as
Distributed Resource Management System (DRMS) module.
"""
#  Copyright 2013 - 2019 Michael K. Schuster
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

from __future__ import print_function

import csv
import datetime
import errno
import math
import os
import re
import threading
import warnings

import bsf.database
import bsf.process

output_directory_name = 'bsfpython_slurm_output'
database_file_name = 'bsfpython_slurm_jobs.db'


class ProcessSLURM(object):
    """C{bsf.drms.slurm.ProcessSLURM} class representing a Simple Linux Utility for Resource Management (SLURM) process.

    The instance variable names result from the SLURM command C{sacct --parsable --long}

    @ivar process_slurm_id: Primary key
    @type process_slurm_id: int | None
    @ivar job_id: The number of the job or job step. It is in the form: job.jobstep
    @type job_id: str | None
    @ivar job_id_raw: The number of the job or job step. It is in the form: I{job.jobstep}
    @type job_id_raw: str | None
    @ivar job_name: The name of the job or job step
    @type job_name: str | None
    @ivar partition: Identifies the partition on which the job ran
    @type partition: str | None
    @ivar max_vm_size: Maximum virtual memory size of all tasks in job
    @type max_vm_size: str | None
    @ivar max_vm_size_node: The node on which the maximum virtual memory size occurred
    @type max_vm_size_node: str | None
    @ivar max_vm_size_task: The task identifier where the maximum virtual memory size occurred
    @type max_vm_size_task: str | None
    @ivar average_vm_size: Average virtual memory size of all tasks in job
    @type average_vm_size: str | None
    @ivar max_rss: Maximum resident set size of all tasks in job
    @type max_rss: str | None
    @ivar max_rss_node: The node on which the maximum resident set size occurred
    @type max_rss_node: str | None
    @ivar max_rss_task: The task identifier where the maximum resident set size occurred
    @type max_rss_task: str | None
    @ivar average_rss: Average resident set size of all tasks in job
    @type average_rss: str | None
    @ivar max_pages: Maximum number of page faults of all tasks in job
    @type max_pages: str | None
    @ivar max_pages_node: The node on which the maximum number of page faults occurred
    @type max_pages_node: str | None
    @ivar max_pages_task: The task identifier where the maximum number of page faults occurred
    @type max_pages_task: str | None
    @ivar average_pages: Average number of page faults of all tasks in job
    @type average_pages: str | None
    @ivar min_cpu: Minimum (system + user) CPU time of all tasks in job
    @type min_cpu: str | None
    @ivar min_cpu_node: The node on which the minimum CPU time occurred
    @type min_cpu_node: str | None
    @ivar min_cpu_task: The task identifier where the minimum CPU time occurred
    @type min_cpu_task: str | None
    @ivar average_cpu: Average (system + user) CPU time of all tasks in job
    @type average_cpu: str | None
    @ivar number_tasks: Total number of tasks in a job or step
    @type number_tasks: str | None
    @ivar allocated_cpus: Count of allocated CPUs
    @type allocated_cpus: str | None
    @ivar elapsed: The jobs elapsed time
    @type elapsed: str | None
    @ivar state: Displays the job status, or state
        Value can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or NODE_FAIL
    @type state: str | None
    @ivar exit_code: The exit code returned by the job script or salloc, typically as set by the exit() function.
        Following the colon is the signal that caused the process to terminate if it was terminated by a signal.
    @type exit_code: str | None
    @ivar average_cpu_frequency: Average weighted CPU frequency of all tasks in job, in kHz
    @type average_cpu_frequency: str | None
    @ivar requested_cpu_frequency_min: Requested minimum CPU frequency for the step, in kHz
    @type requested_cpu_frequency_min: str | None
    @ivar requested_cpu_frequency_max: Requested maximum CPU frequency for the step, in kHz
    @type requested_cpu_frequency_max: str | None
    @ivar requested_cpu_frequency_gov: Requested CPU governor
    @type requested_cpu_frequency_gov: str | None
    @ivar requested_memory: Minimum required memory for the job, in MB
    @type requested_memory: str | None
    @ivar consumed_energy: Total energy consumed by all tasks in job, in joules
    @type consumed_energy: str | None
    @ivar max_disk_read: Maximum number of bytes read by all tasks in job
    @type max_disk_read: str | None
    @ivar max_disk_read_node: The node on which the maximum number of bytes read occurred
    @type max_disk_read_node: str | None
    @ivar max_disk_read_task: The task identifier where the maximum number of bytes read occurred
    @type max_disk_read_task: str | None
    @ivar average_disk_read: Average number of bytes read by all tasks in job
    @type average_disk_read: str | None
    @ivar max_disk_write: Maximum number of bytes written by all tasks in job
    @type max_disk_write: str | None
    @ivar max_disk_write_node: The node on which the maximum number of bytes written occurred
    @type max_disk_write_node: str | None
    @ivar max_disk_write_task: The task identifier where the maximum number of bytes written occurred
    @type max_disk_write_task: str | None
    @ivar average_disk_write: Average number of bytes written by all tasks in job
    @type average_disk_write: str | None
    @ivar allocated_gres: Allocated generic consumable resources
    @type allocated_gres: str | None
    @ivar requested_gres: Requested generic consumable resources
    @type requested_gres: str | None
    @ivar allocated_tres: Allocated trackable resources
    @type allocated_tres: str | None
    @ivar requested_tres: Requested trackable resources
    @type requested_tres: str | None
    """

    def __init__(
            self,
            process_slurm_id=None,
            job_id=None,
            job_id_raw=None,
            job_name=None,
            partition=None,
            max_vm_size=None,
            max_vm_size_node=None,
            max_vm_size_task=None,
            average_vm_size=None,
            max_rss=None,
            max_rss_node=None,
            max_rss_task=None,
            average_rss=None,
            max_pages=None,
            max_pages_node=None,
            max_pages_task=None,
            average_pages=None,
            min_cpu=None,
            min_cpu_node=None,
            min_cpu_task=None,
            average_cpu=None,
            number_tasks=None,
            allocated_cpus=None,
            elapsed=None,
            state=None,
            exit_code=None,
            average_cpu_frequency=None,
            requested_cpu_frequency_min=None,
            requested_cpu_frequency_max=None,
            requested_cpu_frequency_gov=None,
            requested_memory=None,
            consumed_energy=None,
            max_disk_read=None,
            max_disk_read_node=None,
            max_disk_read_task=None,
            average_disk_read=None,
            max_disk_write=None,
            max_disk_write_node=None,
            max_disk_write_task=None,
            average_disk_write=None,
            allocated_gres=None,
            requested_gres=None,
            allocated_tres=None,
            requested_tres=None):
        """Initialise a C{bsf.drms.slurm.ProcessSLURM}.

        @param process_slurm_id:
        @type process_slurm_id: int | None
        @param job_id: The number of the job or job step. It is in the form: I{job.jobstep}
        @type job_id: str | None
        @param job_id_raw: The number of the job or job step. It is in the form: I{job.jobstep}
        @type job_id_raw: str | None
        @param job_name: The name of the job or job step
        @type job_name: str | None
        @param partition: Identifies the partition on which the job ran
        @type partition: str | None
        @param max_vm_size: Maximum virtual memory size of all tasks in job
        @type max_vm_size: str | None
        @param max_vm_size_node: The node on which the maximum virtual memory size occurred
        @type max_vm_size_node: str | None
        @param max_vm_size_task: The task identifier where the maximum virtual memory size occurred
        @type max_vm_size_task: str | None
        @param average_vm_size: Average virtual memory size of all tasks in job
        @type average_vm_size: str | None
        @param max_rss: Maximum resident set size of all tasks in job
        @type max_rss: str | None
        @param max_rss_node: The node on which the maximum resident set size occurred
        @type max_rss_node: str | None
        @param max_rss_task: The task identifier where the maximum resident set size occurred
        @type max_rss_task: str | None
        @param average_rss: Average resident set size of all tasks in job
        @type average_rss: str | None
        @param max_pages: Maximum number of page faults of all tasks in job
        @type max_pages: str | None
        @param max_pages_node: The node on which the maximum number of page faults occurred
        @type max_pages_node: str | None
        @param max_pages_task: The task identifier where the maximum number of page faults occurred
        @type max_pages_task: str | None
        @param average_pages: Average number of page faults of all tasks in job
        @type average_pages: str | None
        @param min_cpu: Minimum (system + user) CPU time of all tasks in job
        @type min_cpu: str | None
        @param min_cpu_node: The node on which the minimum CPU time occurred
        @type min_cpu_node: str | None
        @param min_cpu_task: The task identifier where the minimum CPU time occurred
        @type min_cpu_task: str | None
        @param average_cpu: Average (system + user) CPU time of all tasks in job
        @type average_cpu: str | None
        @param number_tasks: Total number of tasks in a job or step
        @type number_tasks: str | None
        @param allocated_cpus: Count of allocated CPUs
        @type allocated_cpus: str | None
        @param elapsed: The jobs elapsed time
        @type elapsed: str | None
        @param state: Displays the job status, or state.
            Value can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or NODE_FAIL
        @type state: str | None
        @param exit_code: The exit code returned by the job script or salloc, typically as set by the exit() function.
            Following the colon is the signal that caused the process to  terminate if it was terminated by a signal.
        @type exit_code: str | None
        @param average_cpu_frequency: Average weighted CPU frequency of all tasks in job, in kHz
        @type average_cpu_frequency: str | None
        @param requested_cpu_frequency_min: Requested minimum CPU frequency for the step, in kHz
        @type requested_cpu_frequency_min: str | None
        @param requested_cpu_frequency_max: Requested maximum CPU frequency for the step, in kHz
        @type requested_cpu_frequency_max: str | None
        @param requested_cpu_frequency_gov: Requested CPU governor
        @type requested_cpu_frequency_gov: str | None
        @param requested_memory: Minimum required memory for the job, in MB
        @type requested_memory: str | None
        @param consumed_energy: Total energy consumed by all tasks in job, in joules
        @type consumed_energy: str | None
        @param max_disk_read: Maximum number of bytes read by all tasks in job
        @type max_disk_read: str | None
        @param max_disk_read_node: The node on which the maximum number of bytes read occurred
        @type max_disk_read_node: str | None
        @param max_disk_read_task: The task identifier where the maximum number of bytes read occurred
        @type max_disk_read_task: str | None
        @param average_disk_read: Average number of bytes read by all tasks in job
        @type average_disk_read: str | None
        @param max_disk_write: Maximum number of bytes written by all tasks in job
        @type max_disk_write: str | None
        @param max_disk_write_node: The node on which the maximum number of bytes written occurred
        @type max_disk_write_node: str | None
        @param max_disk_write_task: The task identifier where the maximum number of bytes written occurred
        @type max_disk_write_task: str | None
        @param average_disk_write: Average number of bytes written by all tasks in job
        @type average_disk_write: str | None
        @param allocated_gres: Allocated generic consumable resources
        @type allocated_gres: str | None
        @param requested_gres: Requested generic consumable resources
        @type requested_gres: str | None
        @param allocated_tres: Allocated trackable resources
        @type allocated_tres: str | None
        @param requested_tres: Requested trackable resources
        @type requested_tres: str | None
        @return:
        @rtype:
        """
        super(ProcessSLURM, self).__init__()

        self.process_slurm_id = process_slurm_id
        self.job_id = job_id
        self.job_id_raw = job_id_raw
        self.job_name = job_name
        self.partition = partition
        self.max_vm_size = max_vm_size
        self.max_vm_size_node = max_vm_size_node
        self.max_vm_size_task = max_vm_size_task
        self.average_vm_size = average_vm_size
        self.max_rss = max_rss
        self.max_rss_node = max_rss_node
        self.max_rss_task = max_rss_task
        self.average_rss = average_rss
        self.max_pages = max_pages
        self.max_pages_node = max_pages_node
        self.max_pages_task = max_pages_task
        self.average_pages = average_pages
        self.min_cpu = min_cpu
        self.min_cpu_node = min_cpu_node
        self.min_cpu_task = min_cpu_task
        self.average_cpu = average_cpu
        self.number_tasks = number_tasks
        self.allocated_cpus = allocated_cpus
        self.elapsed = elapsed
        self.state = state
        self.exit_code = exit_code
        self.average_cpu_frequency = average_cpu_frequency
        self.requested_cpu_frequency_min = requested_cpu_frequency_min
        self.requested_cpu_frequency_max = requested_cpu_frequency_max
        self.requested_cpu_frequency_gov = requested_cpu_frequency_gov
        self.requested_memory = requested_memory
        self.consumed_energy = consumed_energy
        self.max_disk_read = max_disk_read
        self.max_disk_read_node = max_disk_read_node
        self.max_disk_read_task = max_disk_read_task
        self.average_disk_read = average_disk_read
        self.max_disk_write = max_disk_write
        self.max_disk_write_node = max_disk_write_node
        self.max_disk_write_task = max_disk_write_task
        self.average_disk_write = average_disk_write
        self.allocated_gres = allocated_gres
        self.requested_gres = requested_gres
        self.allocated_tres = allocated_tres
        self.requested_tres = requested_tres

        return


class ProcessSLURMAdaptor(bsf.database.DatabaseAdaptor):
    """C{bsf.drms.slurm.ProcessSLURMAdaptor} class providing database access for the
    C{bsf.drms.slurm.ProcessSLURM} class.

    The SQL column names result from SLURM command sacct --parsable --long
    """

    def __init__(
            self,
            database_connection):
        """Initialise a C{bsf.drms.slurm.ProcessSLURMAdaptor}.

        @param database_connection: C{bsf.database.DatabaseConnection}
        @type database_connection: bsf.database.DatabaseConnection
        @return:
        @rtype:
        """

        super(ProcessSLURMAdaptor, self).__init__(
            database_connection=database_connection,
            object_type=ProcessSLURM,
            table_name='process_slurm',
            column_definition=[
                # Primary key
                ('process_slurm_id', 'INTEGER PRIMARY KEY ASC AUTOINCREMENT'),
                # JobID
                # The number of the job or job step. It is in the form: job.jobstep.
                ('job_id', 'TEXT UNIQUE'),
                # JobIDRaw
                #
                ('job_id_raw', 'TEXT'),
                # JobName
                # The name of the job or job step.
                ('job_name', 'TEXT'),
                # Partition
                # Identifies the partition on which the job ran.
                ('partition', 'TEXT'),
                # MaxVMSize
                # Maximum virtual memory size of all tasks in job.
                ('max_vm_size', 'TEXT'),
                # MaxVMSizeNode
                # The node on which the maximum virtual memory size occurred.
                ('max_vm_size_node', 'TEXT'),
                # MaxVMSizeTask
                # The task identifier where the maximum virtual memory size occurred.
                ('max_vm_size_task', 'TEXT'),
                # AveVMSize
                # Average virtual memory size of all tasks in job.
                ('average_vm_size', 'TEXT'),
                # MaxRSS
                # Maximum resident set size of all tasks in job.
                ('max_rss', 'TEXT'),
                # MaxRSSNode
                # The node on which the maximum resident set size occurred.
                ('max_rss_node', 'TEXT'),
                # MaxRSSTask
                # The task identifier where the maximum resident set size occurred.
                ('max_rss_task', 'TEXT'),
                # AveRSS
                # Average resident set size of all tasks in job.
                ('average_rss', 'TEXT'),
                # MaxPages
                # Maximum number of page faults of all tasks in job.
                ('max_pages', 'TEXT'),
                # MaxPagesNode
                # The node on which the maximum number of page faults occurred.
                ('max_pages_node', 'TEXT'),
                # MaxPagesTask
                # The task identifier where the maximum number of page faults occurred.
                ('max_pages_task', 'TEXT'),
                # AvePages
                # Average number of page faults of all tasks in job.
                ('average_pages', 'TEXT'),
                # MinCPU
                # Minimum (system + user) CPU time of all tasks in job.
                ('min_cpu', 'TEXT'),
                # MinCPUNode
                # The node on which the minimum CPU time occurred.
                ('min_cpu_node', 'TEXT'),
                # MinCPUTask
                # The task identifier where the minimum CPU time occurred.
                ('min_cpu_task', 'TEXT'),
                # AveCPU
                # Average (system + user) CPU time of all tasks in job.
                ('average_cpu', 'TEXT'),
                # NTasks
                # Total number of tasks in a job or step.
                ('number_tasks', 'TEXT'),
                # AllocCPUS
                # Count of allocated CPUs. Equivalent to NCPUS.
                ('allocated_cpus', 'TEXT'),
                # Elapsed
                # The jobs elapsed time.
                ('elapsed', 'TEXT'),
                # State
                # Displays the job status, or state.
                # Value can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or
                # NODE_FAIL.
                ('state', 'TEXT'),
                # ExitCode
                # The exit code returned by the job script or salloc, typically as set by the exit() function.
                # Following the colon is the signal that caused the process to  terminate if it was terminated by
                # a signal.
                ('exit_code', 'TEXT'),
                # AveCPUFreq
                # Average weighted CPU frequency of all tasks in job, in kHz.
                ('average_cpu_frequency', 'TEXT'),
                # ReqCPUFreqMin
                # Requested minimum CPU frequency for the step, in kHz.
                ('requested_cpu_frequency_min', 'TEXT'),
                # ReqCPUFreqMax
                # Requested maximum CPU frequency for the step, in kHz.
                ('requested_cpu_frequency_max', 'TEXT'),
                # ReqCPUFreqGov
                # Requested CPU frequency governor.
                ('requested_cpu_frequency_gov', 'TEXT'),
                # ReqMem
                # Minimum required memory for the job, in MB.
                ('requested_memory', 'TEXT'),
                # ConsumedEnergy
                # Total energy consumed by all tasks in job, in joules.
                ('consumed_energy', 'TEXT'),
                # MaxDiskRead
                # Maximum number of bytes read by all tasks in job.
                ('max_disk_read', 'TEXT'),
                # MaxDiskReadNode
                # The node on which the maximum number of bytes read occurred.
                ('max_disk_read_node', 'TEXT'),
                # MaxDiskReadTask
                # The task identifier where the maximum number of bytes read occurred.
                ('max_disk_read_task', 'TEXT'),
                # AveDiskRead
                # Average number of bytes read by all tasks in job.
                ('average_disk_read', 'TEXT'),
                # MaxDiskWrite
                # Maximum number of bytes written by all tasks in job.
                ('max_disk_write', 'TEXT'),
                # MaxDiskWriteNode
                # The node on which the maximum number of bytes written occurred.
                ('max_disk_write_node', 'TEXT'),
                # MaxDiskWriteTask
                # The task identifier where the maximum number of bytes written occurred.
                ('max_disk_write_task', 'TEXT'),
                # AveDiskWrite
                # Average number of bytes written by all tasks in job.
                ('average_disk_write', 'TEXT'),
                # AllocGRES
                # Names and counts of generic resources allocated.
                ('allocated_gres', 'TEXT'),
                # ReqGRES
                ('requested_gres', 'TEXT'),
                # AllocTRES
                ('allocated_tres', 'TEXT'),
                # ReqTRES
                ('requested_tres', 'TEXT'),
            ])

        # NOTE: Experimentally patch the table definition for this DatabaseAdaptor.
        self.patch_table_definition()

        return

    def patch_table_definition(self):
        """Patch the SQL table definition.

        Re-synchronise the SQLite table with the current BSF Python table definition.

        @return:
        @rtype:
        """
        column_dict_old, column_dict_new = self.compare_table_definitions()

        # print('Column dict old:', column_dict_old)
        # print('Column dict new:', column_dict_new)

        if len(column_dict_new) or len(column_dict_old):
            print('Attempting to patch SQLite table', self.table_name)

            table_name_old = '_'.join((self.table_name, 'old'))

            # Rename the old table to move it sideways.
            statement = self.statement_alter_table_rename(table_name_new=table_name_old)
            # print('Statement:', statement)
            self.get_cursor().execute(statement)

            # Create a new table
            statement = self.statement_create_table()
            # print('Statement:', statement)
            self.get_cursor().execute(statement)

            # Call INSERT INTO ...
            statement_list = list()
            statement_list.append('INSERT')
            statement_list.append('INTO')
            statement_list.append(self.table_name)
            statement_list.append('SELECT')

            select_list = list()
            select_list.append('process_slurm_id')
            select_list.append('job_id')
            select_list.append("''")  # job_id_raw
            select_list.append('job_name')
            select_list.append('partition')

            select_list.append('max_vm_size')
            select_list.append('max_vm_size_node')
            select_list.append('max_vm_size_task')
            select_list.append('average_vm_size')

            select_list.append('max_rss')
            select_list.append('max_rss_node')
            select_list.append('max_rss_task')
            select_list.append('average_rss')

            select_list.append('max_pages')
            select_list.append('max_pages_node')
            select_list.append('max_pages_task')
            select_list.append('average_pages')

            select_list.append('min_cpu')
            select_list.append('min_cpu_node')
            select_list.append('min_cpu_task')
            select_list.append('average_cpu')

            select_list.append('number_tasks')
            select_list.append('allocated_cpus')
            select_list.append('elapsed')
            select_list.append('state')
            select_list.append('exit_code')

            select_list.append('average_cpu_frequency')
            select_list.append('requested_cpu_frequency')  # requested_cpu_frequency_min
            select_list.append("''")  # requested_cpu_frequency_max
            select_list.append("''")  # requested_cpu_frequency_gov
            select_list.append('requested_memory')
            select_list.append('consumed_energy')

            select_list.append('max_disk_read')
            select_list.append('max_disk_read_node')
            select_list.append('max_disk_read_task')
            select_list.append('average_disk_read')

            select_list.append('max_disk_write')
            select_list.append('max_disk_write_node')
            select_list.append('max_disk_write_task')
            select_list.append('average_disk_write')

            select_list.append("''")  # allocated_gres
            select_list.append("''")  # requested_gres
            select_list.append("''")  # allocated_tres
            select_list.append("''")  # requested_tres

            statement_list.append(', '.join(select_list))
            statement_list.append('FROM')
            statement_list.append(table_name_old)

            statement = ' '.join(statement_list)
            # print('Statement:', statement)
            self.get_cursor().execute(statement)
            # By default, the sqlite3 Python module opens transactions implicitly before a
            # Data Modification Language (DML) statement (i.e. INSERT, UPDATE, DELETE or REPLACE)
            # and commits transactions implicitly before a non-DML, non-query statement (i.e. SELECT)
            self.database_connection.commit()

            # Drop the old table
            statement = self.statement_drop_table(table_name=table_name_old)
            # print('Statement:', statement)
            self.get_cursor().execute(statement)

        return

    def select_all_by_job_name(self, name):
        """Select all C{bsf.drms.slurm.ProcessSLURM} objects by I{job_name}.

        The same C{bsf.process.Executable} can be submitted more than once into the C{bsf.Stage}.

        @param name: Job name
        @type name: str
        @return: Python C{list} of C{bsf.drms.slurm.ProcessSLURM} objects
        @rtype: list[bsf.drms.slurm.ProcessSLURM]
        """
        return self.select(statement=self.statement_select(where_clause='job_name = ?'), parameters=[name])

    def select_all_by_state(self, state=None):
        """Select all C{bsf.drms.slurm.ProcessSLURM} objects by I{state}.

        @param state: State
        @type state: str | None
        @return: Python C{list} of C{bsf.drms.slurm.ProcessSLURM} objects
        @rtype: list[bsf.drms.slurm.ProcessSLURM]
        """
        parameters = list()

        if state is None:
            statement = self.statement_select(where_clause='state IS NULL')
        else:
            statement = self.statement_select(where_clause='state = ?')
            parameters.append(state)

        return self.select(statement=statement, parameters=parameters)

    def select_all_by_states(self, state_list, negation=False):
        """Select all C{bsf.drms.slurm.ProcessSLURM} objects by a list of states.

        @param state_list: State
        @type state_list: list[str | unicode]
        @param negation: Negation i.e. SQL NOT IN
        @type negation: bool
        @return: Python C{list} of C{bsf.drms.slurm.ProcessSLURM} objects
        @rtype: list[bsf.drms.slurm.ProcessSLURM]
        """
        if negation:
            statement = self.statement_select(where_clause='state NOT IN (' + ','.join('?' * len(state_list)) + ')')
        else:
            statement = self.statement_select(where_clause='state IN (' + ','.join('?' * len(state_list)) + ')')

        return self.select(statement=statement, parameters=state_list)

    def select_by_job_id(self, job_id):
        """Select one C{bsf.drms.slurm.ProcessSLURM} object by job_id.

        @param job_id: Job identifier
        @type job_id: str
        @return: Python C{list} of C{bsf.drms.slurm.ProcessSLURM} objects
        @rtype: bsf.drms.slurm.ProcessSLURM
        """
        object_list = self.select(statement=self.statement_select(where_clause='job_id = ?'), parameters=[job_id])
        object_length = len(object_list)

        if object_length > 1:
            raise Exception("SQL database returned more than one row for unique field 'job_id'.")
        elif object_length == 1:
            return object_list[0]
        else:
            return


def _recalculate_memory(memory):
    """Recalculate the memory string.

    Multiplier suffixes K, M, G, T, P, E, Z and Y are based on 1024, while
    multiplier suffixes k, m, g, t, p, e, z and y are based on 1000 in line
    with Sun Grid Engine (SGE) conventions.
    https://en.wikipedia.org/wiki/Binary_prefix
    https://en.wikipedia.org/wiki/Metric_prefix

    @param memory: Memory specification string
    @type memory: str
    @return: Memory specification string
    @rtype: str
    """

    assert isinstance(memory, str)

    if memory[-1:] == 'K':  # kibi
        result = float(memory[:-1]) * 1024 ** 1
    elif memory[-1:] == 'k':  # kilo
        result = float(memory[:-1]) * 1000 ** 1
    elif memory[-1:] == 'M':  # mebi
        result = float(memory[:-1]) * 1024 ** 2
    elif memory[-1:] == 'm':  # mega
        result = float(memory[:-1]) * 1000 ** 2
    elif memory[-1:] == 'G':  # gibi
        result = float(memory[:-1]) * 1024 ** 3
    elif memory[-1:] == 'g':  # giga
        result = float(memory[:-1]) * 1000 ** 3
    elif memory[-1:] == 'T':  # tebi
        result = float(memory[:-1]) * 1024 ** 4
    elif memory[-1:] == 't':  # terra
        result = float(memory[:-1]) * 1000 ** 4
    elif memory[-1:] == 'P':  # pebi
        result = float(memory[:-1]) * 1024 ** 5
    elif memory[-1:] == 'p':  # peta
        result = float(memory[:-1]) * 1000 ** 5
    elif memory[-1:] == 'E':  # exbi
        result = float(memory[:-1]) * 1024 ** 6
    elif memory[-1:] == 'e':  # exa
        result = float(memory[:-1]) * 1000 ** 6
    elif memory[-1:] == 'Z':  # zebi
        result = float(memory[:-1]) * 1024 ** 7
    elif memory[-1:] == 'z':  # zetta
        result = float(memory[:-1]) * 1000 ** 7
    elif memory[-1:] == 'Y':  # yobi
        result = float(memory[:-1]) * 1024 ** 8
    elif memory[-1:] == 'y':  # yotta
        result = float(memory[:-1]) * 1000 ** 8
    else:
        result = float(memory)

    # Slurm needs memory specification in MegaBytes.

    result = int(math.ceil(result / 1024 ** 2))

    return str(result)


def submit(stage, debug=0):
    """Submit each C{bsf.process.Executable} of a C{bsf.Stage}.

    Submits each C{bsf.process.Executable} into the
    Simple Linux Utility for Resource Management (SLURM)
    Distributed Resource Management System (DRMS).
    @param stage: C{bsf.Stage}
    @type stage: bsf.Stage
    @param debug: Debug level
    @type debug: int
    @return:
    @rtype:
    """

    def submit_sbatch_stdout(_file_handle, _thread_lock, _debug, _executable):
        """Thread callable to process the SLURM I{sbatch} I{STDOUT} stream.

        Parses the process identifier returned by SLURM sbatch and sets it as
        C{bsf.process.Executable.process_identifier}.
        The response to the SLURM sbatch command looks like:
        Submitted batch job 1234567
        @param _file_handle: File handle (i.e. pipe)
        @type _file_handle: file
        @param _thread_lock: Thread lock
        @type _thread_lock: threading.lock
        @param _debug: Debug level
        @type _debug: int
        @param _executable: C{bsf.process.Executable}
        @type _executable: bsf.process.Executable
        @return:
        @rtype:
        """
        for _line in _file_handle:
            if _debug > 0:
                _thread_lock.acquire(True)
                print('Line:', _line)
                _thread_lock.release()

            _match = re.search(pattern=r'Submitted batch job (\d+)', string=_line)

            if _match:
                _executable.process_identifier = _match.group(1)
            else:
                _thread_lock.acquire(True)
                print('Could not parse the process identifier from the SLURM sbatch response line', _line)
                _thread_lock.release()

        return

    # Open or create a database.
    file_path = os.path.join(stage.working_directory, database_file_name)
    database_connection = bsf.database.DatabaseConnection(file_path=file_path)
    job_submission_adaptor = bsf.database.JobSubmissionAdaptor(database_connection=database_connection)
    process_slurm_adaptor = ProcessSLURMAdaptor(database_connection=database_connection)

    # Abort submission if the database file is not writable.
    # While this test is still prone to a race condition, an eventual database write error can only be caught
    # after a SLURM process has already been submitted. This would lead to an inconsistent state since the first
    # SLURM process would have been submitted and its identifier could no longer be updated. Hence, the test here.

    if not os.access(file_path, os.W_OK):
        raise Exception('Cannot write to SQLite database file path:', repr(file_path))

    output_list = list()
    """ @type output_list: list[str | unicode] """

    output_list.append('#!/usr/bin/env bash\n')
    output_list.append('\n')

    if debug > 0:
        output_list.append('# BSF-Python debug mode: ' + repr(debug) + '\n')
        output_list.append('\n')

    for executable in stage.executable_list:
        executable_drms = bsf.process.Executable(
            name=executable.name,
            program='sbatch',
            sub_command=executable,
            stdout_callable=submit_sbatch_stdout,
            stdout_kwargs={'_executable': executable})

        # Add Stage-specific options.

        # Binary or script

        # Job resource string ...

        # SLURM-specific sanity checks ...

        if stage.memory_limit_hard:
            executable_drms.add_option_pair_long(key='mem', value=_recalculate_memory(stage.memory_limit_hard))
        elif stage.memory_limit_soft:
            executable_drms.add_option_pair_long(key='mem', value=_recalculate_memory(stage.memory_limit_soft))

        if stage.node_list_exclude:
            executable_drms.add_option_pair_long(key='exclude', value=','.join(stage.node_list_exclude))

        if stage.node_list_include:
            executable_drms.add_option_pair_long(key='nodelist', value=','.join(stage.node_list_include))

        if stage.time_limit:
            executable_drms.add_option_pair_long(key='time', value=stage.time_limit)

        # Propagate none of the environment variables.

        executable_drms.add_option_pair_long(key='export', value='NONE')

        # Get the user environment resembling a login shell.

        executable_drms.add_option_pair_long(key='get-user-env', value='L')

        # Parallel environment

        if stage.parallel_environment:
            executable_drms.add_option_pair_long(key='distribution', value=stage.parallel_environment)
            executable_drms.add_option_pair_long(key='ntasks', value='1')
            executable_drms.add_option_pair_long(key='cpus-per-task', value=str(stage.threads))

        executable_drms.add_switch_long(key='requeue')
        executable_drms.add_switch_long(key='share')

        # Queue name

        if stage.queue:
            executable_drms.add_option_pair_long(key='partition', value=stage.queue)

        # Working directory, standard output and standard error streams.

        if stage.working_directory:
            executable_drms.add_option_pair_long(key='workdir', value=stage.working_directory)

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

            executable_drms.add_option_pair_long(
                key='error',
                value=os.path.join(output_directory_name, '_'.join((executable.name, '%j.err'))))

            executable_drms.add_option_pair_long(
                key='output',
                value=os.path.join(output_directory_name, '_'.join((executable.name, '%j.out'))))

        # Job name

        if executable.name:
            executable_drms.add_option_pair_long(key='job-name', value=executable.name)

        # Job hold conditions
        # A particular feature of SLURM is its inability to set process dependencies on process names.
        # Rather, dependencies need setting on the process identifier, which is only obtained after
        # submitting the process. Isn't that exactly what we have a scheduler for? Sigh.
        # Consequently, SLURM process identifiers need to be tracked here, by means of a SQLite database.

        process_identifier_list = list()
        """ @type process_identifier_list: list[str | unicode] """
        for executable_name in executable.dependencies:
            process_slurm_list = process_slurm_adaptor.select_all_by_job_name(name=executable_name)
            if len(process_slurm_list):
                # This bsf.process.Executable has been submitted at least once before.
                # For the moment, set the dependency on the last submission.
                process_identifier_list.append(process_slurm_list[-1].job_id)
            elif debug == 0:
                warnings.warn(
                    'While submitting Executable with name ' + repr(executable.name) + ', ' +
                    'Executable with name ' + repr(executable_name) + ' that it depends on, ' +
                    'has not been submitted before.',
                    UserWarning)
        if len(process_identifier_list):
            # Only set the dependency option if there are some on the process identifier list.
            # The identifier list may be empty if no dependencies exist or no Executable has been submitted before.
            executable_drms.add_option_pair_long(
                key='dependency',
                value=','.join(map(lambda x: 'afterok:' + x, process_identifier_list)))

        # If a process in state 'PENDING' or 'RUNNING' exists already, the Executable does not need submitting.

        process_slurm_list = process_slurm_adaptor.select_all_by_job_name(name=executable.name)
        if len(process_slurm_list) and process_slurm_list[-1].state in ('PENDING', 'RUNNING'):
            executable.submit = False

        # Finally, submit this command if requested and not in debug mode.

        if executable.submit and debug == 0:
            child_return_code = executable_drms.run(debug=debug)

            if child_return_code:
                raise Exception(
                    'SLURM sbatch returned exit code ' + repr(child_return_code) + '\n' +
                    'Command list representation: ' + repr(executable_drms.command_list()))

        # Copy the SLURM command line to the Bash script.

        output_list.append(executable_drms.command_str() + '\n')
        output_list.append('\n')

        # Regardless of an actual bsf.process.Executable submission, UPDATE it in or INSERT it into the SQLite database.

        job_submission = job_submission_adaptor.select_by_name(name=executable.name)
        if job_submission:
            job_submission.command = executable.command_str()
            job_submission_adaptor.update(object_instance=job_submission)
        else:
            job_submission = bsf.database.JobSubmission(
                executable_id=0,
                name=executable.name,
                command=executable.command_str())
            job_submission_adaptor.insert(object_instance=job_submission)

        # Only store a ProcessSLURM, if an bsf.process.Executable has been submitted into SLURM.

        if executable.process_identifier:
            process_slurm = process_slurm_adaptor.select_by_job_id(job_id=executable.process_identifier)
            if not process_slurm:
                process_slurm = ProcessSLURM(job_id=executable.process_identifier, job_name=executable.name)
                process_slurm_adaptor.insert(object_instance=process_slurm)

        # The commit statement should affect both insert statements above.
        database_connection.commit()

    script_path = os.path.join(stage.working_directory, 'bsfpython_slurm_' + stage.name + '.bash')
    with open(script_path, 'wt') as script_file:
        script_file.writelines(output_list)

    return


def check_state(stage, debug=0):
    """Check the state of each C{bsf.process.Executable} of a C{bsf.Stage}.

    @param stage: C{bsf.Stage}
    @type stage: bsf.Stage
    @param debug: Debug level
    @type debug: int
    @return:
    @rtype:
    """

    def check_state_stdout(_stdout_handle, _thread_lock, _debug, _process_slurm_adaptor, _stdout_path=None):
        """Thread callable to process the SLURM I{sacct} I{STDOUT} stream.

        @param _stdout_handle: The I{STDOUT} or I{STDERR} file handle
        @type _stdout_handle: file
        @param _thread_lock: Python C{threading.Lock}
        @type _thread_lock: threading.Lock
        @param _debug: Debug level
        @type _debug: int
        @param _process_slurm_adaptor: C{bsf.drms.slurm.ProcessSLURMAdaptor}
        @type _process_slurm_adaptor: bsf.drms.slurm.ProcessSLURMAdaptor
        @param _stdout_path: I{STDOUT} file path
        @type _stdout_path: str | unicode | None
        @return:
        @rtype:
        """
        if _debug > 0:
            _thread_lock.acquire(True)
            print('[' + datetime.datetime.now().isoformat() + ']',
                  "Started Runner 'STDOUT' processor in module " + repr(__name__) + '.')
            _thread_lock.release()

        if _stdout_path:
            output_file = open(_stdout_path, 'wt')
            if _debug > 0:
                _thread_lock.acquire(True)
                print('[' + datetime.datetime.now().isoformat() + ']',
                      "Opened 'STDOUT' file " + repr(_stdout_path) + '.')
                _thread_lock.release()
        else:
            output_file = None

        dict_reader = csv.DictReader(f=_stdout_handle, delimiter='|')

        for row_dict in dict_reader:
            new_process_slurm = ProcessSLURM(
                job_id=row_dict['JobID'],
                job_id_raw=row_dict['JobIDRaw'],
                job_name=row_dict['JobName'],
                partition=row_dict['Partition'],
                max_vm_size=row_dict['MaxVMSize'],
                max_vm_size_node=row_dict['MaxVMSizeNode'],
                max_vm_size_task=row_dict['MaxVMSizeTask'],
                average_vm_size=row_dict['AveVMSize'],
                max_rss=row_dict['MaxRSS'],
                max_rss_node=row_dict['MaxRSSNode'],
                max_rss_task=row_dict['MaxRSSTask'],
                average_rss=row_dict['AveRSS'],
                max_pages=row_dict['MaxPages'],
                max_pages_node=row_dict['MaxPagesNode'],
                max_pages_task=row_dict['MaxPagesTask'],
                average_pages=row_dict['AvePages'],
                min_cpu=row_dict['MinCPU'],
                min_cpu_node=row_dict['MinCPUNode'],
                min_cpu_task=row_dict['MinCPUTask'],
                average_cpu=row_dict['AveCPU'],
                number_tasks=row_dict['NTasks'],
                allocated_cpus=row_dict['AllocCPUS'],
                elapsed=row_dict['Elapsed'],
                state=row_dict['State'],
                exit_code=row_dict['ExitCode'],
                average_cpu_frequency=row_dict['AveCPUFreq'],
                requested_cpu_frequency_min=row_dict['ReqCPUFreqMin'],
                requested_cpu_frequency_max=row_dict['ReqCPUFreqMax'],
                requested_memory=row_dict['ReqMem'],
                consumed_energy=row_dict['ConsumedEnergy'],
                max_disk_read=row_dict['MaxDiskRead'],
                max_disk_read_node=row_dict['MaxDiskReadNode'],
                max_disk_read_task=row_dict['MaxDiskReadTask'],
                average_disk_read=row_dict['AveDiskRead'],
                max_disk_write=row_dict['MaxDiskWrite'],
                max_disk_write_node=row_dict['MaxDiskWriteNode'],
                max_disk_write_task=row_dict['MaxDiskWriteTask'],
                average_disk_write=row_dict['AveDiskWrite'],
                allocated_gres=row_dict['AllocGRES'],
                requested_gres=row_dict['ReqGRES'],
                allocated_tres=row_dict['AllocTRES'],
                requested_tres=row_dict['ReqTRES'])

            # Check if the ProcessSLURM already exists.
            old_process_slurm = _process_slurm_adaptor.select_by_job_id(job_id=new_process_slurm.job_id)

            if old_process_slurm is None:
                # The JobID is not in the database, which is caused by job_id.batch entries.
                # Insert them afterwards.
                _thread_lock.acquire(True)
                _process_slurm_adaptor.insert(object_instance=new_process_slurm)
                _thread_lock.release()
            else:
                new_process_slurm.process_slurm_id = old_process_slurm.process_slurm_id
                _thread_lock.acquire(True)
                _process_slurm_adaptor.update(object_instance=new_process_slurm)
                _thread_lock.release()

        if _debug > 0:
            _thread_lock.acquire(True)
            print('[' + datetime.datetime.now().isoformat() + ']',
                  "Received EOF on 'STDOUT' pipe.")
            _thread_lock.release()

        if output_file:
            output_file.close()

            if _debug > 0:
                _thread_lock.acquire(True)
                print('[' + datetime.datetime.now().isoformat() + ']',
                      "Closed 'STDOUT' file " + repr(_stdout_path) + '.')
                _thread_lock.release()

        # Commit changes to the database and explicitly disconnect so that other threads have access.
        _process_slurm_adaptor.commit()
        _process_slurm_adaptor.disconnect()

        return

    # Open or create a database.
    database_connection = bsf.database.DatabaseConnection(
        file_path=os.path.join(stage.working_directory, database_file_name))
    process_slurm_adaptor = ProcessSLURMAdaptor(database_connection=database_connection)

    process_slurm_list = list()
    """ @type process_slurm_list: list[ProcessSLURM] """
    # Get all Processes from the database that have no state set.
    process_slurm_list.extend(process_slurm_adaptor.select_all_by_state(state=None))
    # Get all processes from the database that have a state of 'PENDING' or 'RUNNING' as those need checking.
    process_slurm_list.extend(process_slurm_adaptor.select_all_by_states(state_list=['PENDING', 'RUNNING']))

    # Explicitly disconnect here so that the Thread can access the database.

    process_slurm_adaptor.disconnect()

    # Return if no ProcessSLURM objects need updating.
    if not process_slurm_list:
        return

    if debug > 0:
        print('Number of ProcessSLURM objects to check:', len(process_slurm_list))

    executable_drms = bsf.process.Executable(
        name='sacct',
        program='sacct',
        stdout_callable=check_state_stdout,
        stdout_kwargs={'_process_slurm_adaptor': process_slurm_adaptor})
    executable_drms.add_option_long(key='jobs', value=','.join(map(lambda x: x.job_id, process_slurm_list)))
    executable_drms.add_switch_long(key='long')
    executable_drms.add_switch_long(key='parsable')

    child_return_code = executable_drms.run(debug=debug)

    if child_return_code:
        raise Exception(
            'SLURM sacct returned exit code ' + repr(child_return_code) + '\n' +
            'Command list representation: ' + repr(executable_drms.command_list()))

    return
