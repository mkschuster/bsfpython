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
"""The :py:mod:`bsf.drms.slurm` module supports the Simple Linux Utility for Resource Management (SLURM) system as
Distributed Resource Management System (DRMS) module.
"""
import errno
import logging
import math
import os
import re
from csv import DictReader
from threading import Lock
from typing import Optional, TextIO

from bsf.analysis import Stage
from bsf.connector import StandardOutputStream
from bsf.database import DatabaseAdaptor, DatabaseConnection, JobSubmission, JobSubmissionAdaptor, \
    SQLiteTableInfoAdaptor
from bsf.process import Executable

database_file_name = 'bsfpython_slurm_jobs.db'

module_logger = logging.getLogger(name=__name__)

output_directory_name = 'bsfpython_slurm_output'


class ProcessSLURM(object):
    """The :py:class:`bsf.drms.slurm.ProcessSLURM` class represents a
    :literal:`Simple Linux Utility for Resource Management (SLURM)` process.

    The instance variable names result from the SLURM command :literal:`sacct --parsable --long`.

    :ivar process_slurm_id: Primary key
    :type process_slurm_id: int | None
    :ivar job_id: The number of the job or job step. It is in the form: :literal:`job.jobstep`.
    :type job_id: str | None
    :ivar job_id_raw: The number of the job or job step. It is in the form: :literal:`job.jobstep`.
    :type job_id_raw: str | None
    :ivar job_name: The name of the job or job step
    :type job_name: str | None
    :ivar partition: Identifies the partition on which the job ran
    :type partition: str | None
    :ivar max_vm_size: Maximum virtual memory size of all tasks in job
    :type max_vm_size: str | None
    :ivar max_vm_size_node: The node on which the maximum virtual memory size occurred
    :type max_vm_size_node: str | None
    :ivar max_vm_size_task: The task identifier where the maximum virtual memory size occurred
    :type max_vm_size_task: str | None
    :ivar average_vm_size: Average virtual memory size of all tasks in job
    :type average_vm_size: str | None
    :ivar max_rss: Maximum resident set size of all tasks in job
    :type max_rss: str | None
    :ivar max_rss_node: The node on which the maximum resident set size occurred
    :type max_rss_node: str | None
    :ivar max_rss_task: The task identifier where the maximum resident set size occurred
    :type max_rss_task: str | None
    :ivar average_rss: Average resident set size of all tasks in job
    :type average_rss: str | None
    :ivar max_pages: Maximum number of page faults of all tasks in job
    :type max_pages: str | None
    :ivar max_pages_node: The node on which the maximum number of page faults occurred
    :type max_pages_node: str | None
    :ivar max_pages_task: The task identifier where the maximum number of page faults occurred
    :type max_pages_task: str | None
    :ivar average_pages: Average number of page faults of all tasks in job
    :type average_pages: str | None
    :ivar min_cpu: Minimum (system + user) CPU time of all tasks in job
    :type min_cpu: str | None
    :ivar min_cpu_node: The node on which the minimum CPU time occurred
    :type min_cpu_node: str | None
    :ivar min_cpu_task: The task identifier where the minimum CPU time occurred
    :type min_cpu_task: str | None
    :ivar average_cpu: Average (system + user) CPU time of all tasks in job
    :type average_cpu: str | None
    :ivar number_tasks: Total number of tasks in a job or step
    :type number_tasks: str | None
    :ivar allocated_cpus: Count of allocated CPUs
    :type allocated_cpus: str | None
    :ivar elapsed: The jobs elapsed time
    :type elapsed: str | None
    :ivar state: Displays the job status, or state
        Value can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or NODE_FAIL
    :type state: str | None
    :ivar exit_code: The exit code returned by the job script or salloc, typically as set by the exit() function.
        Following the colon is the signal that caused the process to terminate if it was terminated by a signal.
    :type exit_code: str | None
    :ivar average_cpu_frequency: Average weighted CPU frequency of all tasks in job, in kHz
    :type average_cpu_frequency: str | None
    :ivar requested_cpu_frequency_min: Requested minimum CPU frequency for the step, in kHz
    :type requested_cpu_frequency_min: str | None
    :ivar requested_cpu_frequency_max: Requested maximum CPU frequency for the step, in kHz
    :type requested_cpu_frequency_max: str | None
    :ivar requested_cpu_frequency_gov: Requested CPU governor
    :type requested_cpu_frequency_gov: str | None
    :ivar requested_memory: Minimum required memory for the job, in MB
    :type requested_memory: str | None
    :ivar consumed_energy: Total energy consumed by all tasks in job, in joules
    :type consumed_energy: str | None
    :ivar max_disk_read: Maximum number of bytes read by all tasks in job
    :type max_disk_read: str | None
    :ivar max_disk_read_node: The node on which the maximum number of bytes read occurred
    :type max_disk_read_node: str | None
    :ivar max_disk_read_task: The task identifier where the maximum number of bytes read occurred
    :type max_disk_read_task: str | None
    :ivar average_disk_read: Average number of bytes read by all tasks in job
    :type average_disk_read: str | None
    :ivar max_disk_write: Maximum number of bytes written by all tasks in job
    :type max_disk_write: str | None
    :ivar max_disk_write_node: The node on which the maximum number of bytes written occurred
    :type max_disk_write_node: str | None
    :ivar max_disk_write_task: The task identifier where the maximum number of bytes written occurred
    :type max_disk_write_task: str | None
    :ivar average_disk_write: Average number of bytes written by all tasks in job
    :type average_disk_write: str | None
    :ivar allocated_gres: Allocated generic consumable resources
    :type allocated_gres: str | None
    :ivar requested_gres: Requested generic consumable resources
    :type requested_gres: str | None
    :ivar allocated_tres: Allocated trackable resources
    :type allocated_tres: str | None
    :ivar requested_tres: Requested trackable resources
    :type requested_tres: str | None
    :ivar tres_usage_in_average: Tres average usage in by all tasks in job.
    :type tres_usage_in_average: str | None
    :ivar tres_usage_in_maximum: Tres maximum usage in by all tasks in job.
    :type tres_usage_in_maximum: str | None
    :ivar tres_usage_in_maximum_node: Node for which each maximum TRES usage out occurred.
    :type tres_usage_in_maximum_node: str | None
    :ivar tres_usage_in_maximum_task: Task for which each maximum TRES usage out occurred.
    :type tres_usage_in_maximum_task: str | None
    :ivar tres_usage_in_minimum: Tres minimum usage in by all tasks in job.
    :type tres_usage_in_minimum: str | None
    :ivar tres_usage_in_minimum_node: Node for which each minimum TRES usage out occurred.
    :type tres_usage_in_minimum_node: str | None
    :ivar tres_usage_in_minimum_task: Task for which each minimum TRES usage out occurred.
    :type tres_usage_in_minimum_task: str | None
    :ivar tres_usage_in_total: Tres total usage in by all tasks in job.
    :type tres_usage_in_total: str | None
    :ivar tres_usage_out_maximum: Tres maximum usage out by all tasks in job.
    :type tres_usage_out_maximum: str | None
    :ivar tres_usage_out_maximum_node: Node for which each maximum TRES usage out occurred.
    :type tres_usage_out_maximum_node: str | None
    :ivar tres_usage_out_maximum_task: Task for which each maximum TRES usage out occurred.
    :type tres_usage_out_maximum_task: str | None
    :ivar tres_usage_out_average: Tres average usage out by all tasks in job.
    :type tres_usage_out_average: str | None
    :ivar tres_usage_out_total: Tres total usage out by all tasks in job.
    :type tres_usage_out_total: str | None
    """

    def __init__(
            self,
            process_slurm_id: Optional[int] = None,
            job_id: Optional[str] = None,
            job_id_raw: Optional[str] = None,
            job_name: Optional[str] = None,
            partition: Optional[str] = None,
            max_vm_size: Optional[str] = None,
            max_vm_size_node: Optional[str] = None,
            max_vm_size_task: Optional[str] = None,
            average_vm_size: Optional[str] = None,
            max_rss: Optional[str] = None,
            max_rss_node: Optional[str] = None,
            max_rss_task: Optional[str] = None,
            average_rss: Optional[str] = None,
            max_pages: Optional[str] = None,
            max_pages_node: Optional[str] = None,
            max_pages_task: Optional[str] = None,
            average_pages: Optional[str] = None,
            min_cpu: Optional[str] = None,
            min_cpu_node: Optional[str] = None,
            min_cpu_task: Optional[str] = None,
            average_cpu: Optional[str] = None,
            number_tasks: Optional[str] = None,
            allocated_cpus: Optional[str] = None,
            elapsed: Optional[str] = None,
            state: Optional[str] = None,
            exit_code: Optional[str] = None,
            average_cpu_frequency: Optional[str] = None,
            requested_cpu_frequency_min: Optional[str] = None,
            requested_cpu_frequency_max: Optional[str] = None,
            requested_cpu_frequency_gov: Optional[str] = None,
            requested_memory: Optional[str] = None,
            consumed_energy: Optional[str] = None,
            max_disk_read: Optional[str] = None,
            max_disk_read_node: Optional[str] = None,
            max_disk_read_task: Optional[str] = None,
            average_disk_read: Optional[str] = None,
            max_disk_write: Optional[str] = None,
            max_disk_write_node: Optional[str] = None,
            max_disk_write_task: Optional[str] = None,
            average_disk_write: Optional[str] = None,
            allocated_gres: Optional[str] = None,
            requested_gres: Optional[str] = None,
            allocated_tres: Optional[str] = None,
            requested_tres: Optional[str] = None,
            tres_usage_in_average: Optional[str] = None,
            tres_usage_in_maximum: Optional[str] = None,
            tres_usage_in_maximum_node: Optional[str] = None,
            tres_usage_in_maximum_task: Optional[str] = None,
            tres_usage_in_minimum: Optional[str] = None,
            tres_usage_in_minimum_node: Optional[str] = None,
            tres_usage_in_minimum_task: Optional[str] = None,
            tres_usage_in_total: Optional[str] = None,
            tres_usage_out_maximum: Optional[str] = None,
            tres_usage_out_maximum_node: Optional[str] = None,
            tres_usage_out_maximum_task: Optional[str] = None,
            tres_usage_out_average: Optional[str] = None,
            tres_usage_out_total: Optional[str] = None) -> None:
        """Initialise a :py:class:`bsf.drms.slurm.ProcessSLURM` object.

        :param process_slurm_id:
        :type process_slurm_id: int | None
        :param job_id: The number of the job or job step. It is in the form: :literal:`job.jobstep`.
        :type job_id: str | None
        :param job_id_raw: The number of the job or job step. It is in the form: :literal:`job.jobstep`.
        :type job_id_raw: str | None
        :param job_name: The name of the job or job step
        :type job_name: str | None
        :param partition: Identifies the partition on which the job ran
        :type partition: str | None
        :param max_vm_size: Maximum virtual memory size of all tasks in job
        :type max_vm_size: str | None
        :param max_vm_size_node: The node on which the maximum virtual memory size occurred
        :type max_vm_size_node: str | None
        :param max_vm_size_task: The task identifier where the maximum virtual memory size occurred
        :type max_vm_size_task: str | None
        :param average_vm_size: Average virtual memory size of all tasks in job
        :type average_vm_size: str | None
        :param max_rss: Maximum resident set size of all tasks in job
        :type max_rss: str | None
        :param max_rss_node: The node on which the maximum resident set size occurred
        :type max_rss_node: str | None
        :param max_rss_task: The task identifier where the maximum resident set size occurred
        :type max_rss_task: str | None
        :param average_rss: Average resident set size of all tasks in job
        :type average_rss: str | None
        :param max_pages: Maximum number of page faults of all tasks in job
        :type max_pages: str | None
        :param max_pages_node: The node on which the maximum number of page faults occurred
        :type max_pages_node: str | None
        :param max_pages_task: The task identifier where the maximum number of page faults occurred
        :type max_pages_task: str | None
        :param average_pages: Average number of page faults of all tasks in job
        :type average_pages: str | None
        :param min_cpu: Minimum (system + user) CPU time of all tasks in job
        :type min_cpu: str | None
        :param min_cpu_node: The node on which the minimum CPU time occurred
        :type min_cpu_node: str | None
        :param min_cpu_task: The task identifier where the minimum CPU time occurred
        :type min_cpu_task: str | None
        :param average_cpu: Average (system + user) CPU time of all tasks in job
        :type average_cpu: str | None
        :param number_tasks: Total number of tasks in a job or step
        :type number_tasks: str | None
        :param allocated_cpus: Count of allocated CPUs
        :type allocated_cpus: str | None
        :param elapsed: The jobs elapsed time
        :type elapsed: str | None
        :param state: Displays the job status, or state.
            Value can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or NODE_FAIL
        :type state: str | None
        :param exit_code: The exit code returned by the job script or salloc, typically as set by the exit() function.
            Following the colon is the signal that caused the process to  terminate if it was terminated by a signal.
        :type exit_code: str | None
        :param average_cpu_frequency: Average weighted CPU frequency of all tasks in job, in kHz
        :type average_cpu_frequency: str | None
        :param requested_cpu_frequency_min: Requested minimum CPU frequency for the step, in kHz
        :type requested_cpu_frequency_min: str | None
        :param requested_cpu_frequency_max: Requested maximum CPU frequency for the step, in kHz
        :type requested_cpu_frequency_max: str | None
        :param requested_cpu_frequency_gov: Requested CPU governor
        :type requested_cpu_frequency_gov: str | None
        :param requested_memory: Minimum required memory for the job, in MB
        :type requested_memory: str | None
        :param consumed_energy: Total energy consumed by all tasks in job, in joules
        :type consumed_energy: str | None
        :param max_disk_read: Maximum number of bytes read by all tasks in job
        :type max_disk_read: str | None
        :param max_disk_read_node: The node on which the maximum number of bytes read occurred
        :type max_disk_read_node: str | None
        :param max_disk_read_task: The task identifier where the maximum number of bytes read occurred
        :type max_disk_read_task: str | None
        :param average_disk_read: Average number of bytes read by all tasks in job
        :type average_disk_read: str | None
        :param max_disk_write: Maximum number of bytes written by all tasks in job
        :type max_disk_write: str | None
        :param max_disk_write_node: The node on which the maximum number of bytes written occurred
        :type max_disk_write_node: str | None
        :param max_disk_write_task: The task identifier where the maximum number of bytes written occurred
        :type max_disk_write_task: str | None
        :param average_disk_write: Average number of bytes written by all tasks in job
        :type average_disk_write: str | None
        :param allocated_gres: Allocated generic consumable resources
        :type allocated_gres: str | None
        :param requested_gres: Requested generic consumable resources
        :type requested_gres: str | None
        :param allocated_tres: Allocated trackable resources
        :type allocated_tres: str | None
        :param requested_tres: Requested trackable resources
        :type requested_tres: str | None
        :param tres_usage_in_average: Tres average usage in by all tasks in job.
        :type tres_usage_in_average: str | None
        :param tres_usage_in_maximum: Tres maximum usage in by all tasks in job.
        :type tres_usage_in_maximum: str | None
        :param tres_usage_in_maximum_node: Node for which each maximum TRES usage out occurred.
        :type tres_usage_in_maximum_node: str | None
        :param tres_usage_in_maximum_task: Task for which each maximum TRES usage out occurred.
        :type tres_usage_in_maximum_task: str | None
        :param tres_usage_in_minimum: Tres minimum usage in by all tasks in job.
        :type tres_usage_in_minimum: str | None
        :param tres_usage_in_minimum_node: Node for which each minimum TRES usage out occurred.
        :type tres_usage_in_minimum_node: str | None
        :param tres_usage_in_minimum_task: Task for which each minimum TRES usage out occurred.
        :type tres_usage_in_minimum_task: str | None
        :param tres_usage_in_total: Tres total usage in by all tasks in job.
        :type tres_usage_in_total: str | None
        :param tres_usage_out_maximum: Tres maximum usage out by all tasks in job.
        :type tres_usage_out_maximum: str | None
        :param tres_usage_out_maximum_node: Node for which each maximum TRES usage out occurred.
        :type tres_usage_out_maximum_node: str | None
        :param tres_usage_out_maximum_task: Task for which each maximum TRES usage out occurred.
        :type tres_usage_out_maximum_task: str | None
        :param tres_usage_out_average: Tres average usage out by all tasks in job.
        :type tres_usage_out_average: str | None
        :param tres_usage_out_total: Tres total usage out by all tasks in job.
        :type tres_usage_out_total: str | None
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
        self.tres_usage_in_average = tres_usage_in_average
        self.tres_usage_in_maximum = tres_usage_in_maximum
        self.tres_usage_in_maximum_node = tres_usage_in_maximum_node
        self.tres_usage_in_maximum_task = tres_usage_in_maximum_task
        self.tres_usage_in_minimum = tres_usage_in_minimum
        self.tres_usage_in_minimum_node = tres_usage_in_minimum_node
        self.tres_usage_in_minimum_task = tres_usage_in_minimum_task
        self.tres_usage_in_total = tres_usage_in_total
        self.tres_usage_out_maximum = tres_usage_out_maximum
        self.tres_usage_out_maximum_node = tres_usage_out_maximum_node
        self.tres_usage_out_maximum_task = tres_usage_out_maximum_task
        self.tres_usage_out_average = tres_usage_out_average
        self.tres_usage_out_total = tres_usage_out_total

        return


class ProcessSLURMAdaptor(DatabaseAdaptor):
    """The :py:class:`bsf.drms.slurm.ProcessSLURMAdaptor` class provides database access for the
    :py:class:`bsf.drms.slurm.ProcessSLURM` class.

    The SQL column names result from SLURM command :literal:`sacct --parsable --long`
    """

    def __init__(
            self,
            database_connection: DatabaseConnection) -> None:
        """Initialise a :py:class:`bsf.drms.slurm.ProcessSLURMAdaptor` object.

        :param database_connection: A :py:class:`bsf.database.DatabaseConnection` object.
        :type database_connection: DatabaseConnection
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
                # The identification number of the job or job step.
                # Prints the JobID in the form JobID[.JobStep] for regular, heterogeneous and array jobs.
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
                # NOTE: In SLURM 19, but no longer in SLURM 20.
                ('allocated_gres', 'TEXT'),
                # ReqGRES
                # Names and counts of generic resources requested.
                # NOTE: In SLURM 19, but no longer in SLURM 20.
                ('requested_gres', 'TEXT'),
                # AllocTRES
                # Trackable resources.
                # These are the resources allocated to the job/step after the job started running.
                ('allocated_tres', 'TEXT'),
                # ReqTRES
                # Trackable resources.
                # These are the minimum resource counts requested by the job/step at submission time.
                ('requested_tres', 'TEXT'),
                # TresUsageInAve
                # Tres average usage in by all tasks in job.
                ('tres_usage_in_average', 'TEXT'),
                # TresUsageInMax
                # Tres maximum usage in by all tasks in job.
                ('tres_usage_in_maximum', 'TEXT'),
                # TresUsageInMaxNode
                # Node for which each maximum TRES usage out occurred.
                ('tres_usage_in_maximum_node', 'TEXT'),
                # TresUsageInMaxTask
                # Task for which each maximum TRES usage out occurred.
                ('tres_usage_in_maximum_task', 'TEXT'),
                # TresUsageInMin
                # Tres minimum usage in by all tasks in job.
                ('tres_usage_in_minimum', 'TEXT'),
                # TresUsageInMinNode
                # Node for which each minimum TRES usage out occurred.
                ('tres_usage_in_minimum_node', 'TEXT'),
                # TresUsageInMinTask
                # Task for which each minimum TRES usage out occurred.
                ('tres_usage_in_minimum_task', 'TEXT'),
                # TresUsageInTot
                # Tres total usage in by all tasks in job.
                ('tres_usage_in_total', 'TEXT'),
                # TresUsageOutMax
                # Tres maximum usage out by all tasks in job.
                ('tres_usage_out_maximum', 'TEXT'),
                # TresUsageOutMaxNode
                # Node for which each maximum TRES usage out occurred.
                ('tres_usage_out_maximum_node', 'TEXT'),
                # TresUsageOutMaxTask
                # Task for which each maximum TRES usage out occurred.
                ('tres_usage_out_maximum_task', 'TEXT'),
                # TresUsageOutAve,
                # Tres average usage out by all tasks in job.
                ('tres_usage_out_average', 'TEXT'),
                # TresUsageOutTot
                # Tres total usage out by all tasks in job.
                ('tres_usage_out_total', 'TEXT'),
            ])

        # NOTE: Experimentally patch the table definition for this DatabaseAdaptor.
        self.patch_table_definition()

        return

    def patch_table_definition(self) -> None:
        """Patch the SQL table definition.

        Re-synchronise the SQLite table with the current BSF Python table definition.
        """
        column_dict_old, column_dict_new = self.compare_table_definitions()

        if len(column_dict_new) or len(column_dict_old):
            module_logger.info('Attempting to patch SQLite table: %r.', self.table_name)

            # Get the column definitions from the original table, which is needed for the
            # INSERT INTO (columns) expression.
            pragma_table_info_adaptor = SQLiteTableInfoAdaptor(
                database_connection=self.database_connection)

            pragma_table_info_list = pragma_table_info_adaptor.select_all_by_table_name(
                table_name=self.table_name)

            column_expression_old = ', '.join(map(lambda x: x.column_name, pragma_table_info_list))

            table_name_old = '_'.join((self.table_name, 'old'))

            # Rename the old table to move it sideways.
            statement = self.statement_alter_table_rename(table_name_new=table_name_old)
            module_logger.log(logging.DEBUG - 1, 'SQL statement: %r', statement)
            self.get_cursor().execute(statement)

            # Create a new table
            statement = self.statement_create_table()
            module_logger.log(logging.DEBUG - 1, 'SQL statement: %r', statement)
            self.get_cursor().execute(statement)

            # Call INSERT INTO ...
            statement_list = list()
            statement_list.append('INSERT')
            statement_list.append('INTO')
            statement_list.append(self.table_name)
            statement_list.append('(' + column_expression_old + ')')
            statement_list.append('SELECT')
            statement_list.append(column_expression_old)
            statement_list.append('FROM')
            statement_list.append(table_name_old)

            statement = ' '.join(statement_list)
            module_logger.log(logging.DEBUG - 1, 'SQL statement: %r', statement)
            self.get_cursor().execute(statement)
            # By default, the sqlite3 Python module opens transactions implicitly before a
            # Data Modification Language (DML) statement (i.e., INSERT, UPDATE, DELETE or REPLACE)
            # and commits transactions implicitly before a non-DML, non-query statement (i.e., SELECT)
            self.database_connection.commit()

            # Drop the old table
            statement = self.statement_drop_table(table_name=table_name_old)
            module_logger.log(logging.DEBUG - 1, 'SQL statement: %r', statement)
            self.get_cursor().execute(statement)

        return

    def select_all_by_job_name(self, name: str) -> list[ProcessSLURM]:
        """Select all :py:class:`bsf.drms.slurm.ProcessSLURM` objects by :literal:`job_name`.

        The same :py:class:`bsf.process.Executable` object can be submitted more than once into the
        :py:class:`bsf.analysis.Stage` object.

        :param name: Job name.
        :type name: str
        :return: Python :py:class:`list` object of :py:class:`bsf.drms.slurm.ProcessSLURM` objects.
        :rtype: list[ProcessSLURM]
        """
        return self.select(statement=self.statement_select(where_clause='job_name = ?'), parameters=[name])

    def select_all_by_state(self, state: Optional[str] = None) -> list[ProcessSLURM]:
        """Select all :py:class:`bsf.drms.slurm.ProcessSLURM` objects by :literal:`state`.

        :param state: State
        :type state: str | None
        :return: A Python :py:class:`list` object of :py:class:`bsf.drms.slurm.ProcessSLURM` objects.
        :rtype: list[ProcessSLURM]
        """
        parameters = list()

        if state is None:
            statement = self.statement_select(where_clause='state IS NULL')
        else:
            statement = self.statement_select(where_clause='state = ?')
            parameters.append(state)

        return self.select(statement=statement, parameters=parameters)

    def select_all_by_states(self, state_list: list[str], negation: bool = False) -> list[ProcessSLURM]:
        """Select all :py:class:`bsf.drms.slurm.ProcessSLURM` objects by a list of :literal:`states`.

        :param state_list: A Python :class:`list` object of Python :py:class:`str` (state) objects.
        :type state_list: list[str]
        :param negation: Negation i.e., SQL NOT IN
        :type negation: bool
        :return: A Python :py:class:`list` object of :py:class:`bsf.drms.slurm.ProcessSLURM` objects.
        :rtype: list[ProcessSLURM]
        """
        if negation:
            statement = self.statement_select(where_clause='state NOT IN (' + ','.join('?' * len(state_list)) + ')')
        else:
            statement = self.statement_select(where_clause='state IN (' + ','.join('?' * len(state_list)) + ')')

        return self.select(statement=statement, parameters=state_list)

    def select_by_job_id(self, job_id: str) -> Optional[ProcessSLURM]:
        """Select one :py:class:`bsf.drms.slurm.ProcessSLURM` object by :literal:`job_id`.

        :param job_id: Job identifier
        :type job_id: str
        :return: A Python :py:class:`list` object of :py:class:`bsf.drms.slurm.ProcessSLURM` objects.
        :rtype: ProcessSLURM | None
        """
        object_list = self.select(statement=self.statement_select(where_clause='job_id = ?'), parameters=[job_id])
        object_length = len(object_list)

        if object_length > 1:
            raise Exception(f"The SQL database returned more than one row for unique field 'job_id' "
                            f"in table {self.table_name!r}.")
        elif object_length == 1:
            return object_list[0]
        else:
            return


def _recalculate_memory(memory: str) -> str:
    """Recalculate a memory string.

    Multiplier suffixes K, M, G, T, P, E, Z and Y are based on 1024, while
    multiplier suffixes k, m, g, t, p, e, z and y are based on 1000 in line
    with Sun Grid Engine (SGE) conventions.
    https://en.wikipedia.org/wiki/Binary_prefix
    https://en.wikipedia.org/wiki/Metric_prefix

    :param memory: Memory specification string
    :type memory: str
    :return: Memory specification string
    :rtype: str
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


def submit(stage: Stage, drms_submit: Optional[bool] = None, write_script: Optional[bool] = None) -> None:
    """Submit each :py:class:`bsf.process.Executable` object of a :py:class:`bsf.analysis.Stage` object.

    Submits each :py:class:`bsf.process.Executable` object into the
    :emphasis:`Simple Linux Utility for Resource Management` (SLURM)
    :emphasis:`Distributed Resource Management System` (DRMS).

    :param stage: A :py:class:`bsf.analysis.Stage` object.
    :type stage: Stage
    :param drms_submit: Submit to the :emphasis:`Distributed Resource Management System` (DRMS).
    :type drms_submit: bool | None
    :param write_script: Write a :py:class:`bsf.analysis.Stage`-specific GNU Bash script.
    :type write_script: bool | None
    """

    def submit_sbatch_stdout(_file_handle: TextIO, _thread_lock: Lock, _executable: Executable):
        """Thread callable to process the SLURM :manpage:`sbatch(1)` :literal:`STDOUT` stream.

        Parses the process identifier returned by SLURM :literal:`sbatch` and sets it as
        :py:attr:`bsf.process.Executable.process_identifier`.
        The response to the SLURM :literal:`sbatch` command looks like:

        :literal:`Submitted batch job 1234567`

        :param _file_handle: File handle (i.e. pipe)
        :type _file_handle: TextIO
        :param _thread_lock: A Python :py:class:`threading.Lock` object.
        :type _thread_lock: Lock
        :param _executable: A :py:class:`bsf.process.Executable` object.
        :type _executable: Executable
        """
        for _line in _file_handle:
            _line = _line.rstrip()

            module_logger.debug('Line: %r', _line)

            _match = re.search(pattern=r'Submitted batch job (\d+)', string=_line)

            if _match:
                _executable.process_identifier = _match.group(1)
            else:
                module_logger.warning('Could not parse the SLURM sbatch response line: %r', _line)

        return

    # Open or create a database.
    file_path = os.path.join(stage.working_directory, database_file_name)
    database_connection = DatabaseConnection(file_path=file_path)
    job_submission_adaptor = JobSubmissionAdaptor(database_connection=database_connection)
    process_slurm_adaptor = ProcessSLURMAdaptor(database_connection=database_connection)

    # Abort submission if the database file is not writable.
    # While this test is still prone to a race condition, an eventual database write error can only be caught
    # after a SLURM process has already been submitted. This would lead to an inconsistent state since the first
    # SLURM process would have been submitted and its identifier could no longer be updated. Hence, the test here.

    if not os.access(file_path, os.W_OK):
        raise Exception(f'Cannot write to SQLite database file path: {file_path!r}')

    output_list: list[str] = list()

    output_list.append('#!/usr/bin/env bash\n')
    output_list.append('\n')

    for executable in stage.executable_list:
        executable_drms = Executable(
            name=executable.name,
            program='sbatch',
            sub_command=executable,
            stdout=StandardOutputStream(
                thread_callable=submit_sbatch_stdout,
                thread_kwargs={'_executable': executable}))

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

        # Parallel environment

        if stage.parallel_environment:
            executable_drms.add_option_pair_long(key='distribution', value=stage.parallel_environment)
            executable_drms.add_option_pair_long(key='ntasks', value='1')
            executable_drms.add_option_pair_long(key='cpus-per-task', value=str(stage.threads))

        executable_drms.add_switch_long(key='requeue')
        executable_drms.add_switch_long(key='oversubscribe')

        # Queue name

        if stage.queue:
            executable_drms.add_option_pair_long(key='partition', value=stage.queue)
            # FIXME: This is a hack to cope with the new cluster configuration.
            # TODO: This should be modelled via a DMRS-specific class.
            executable_drms.add_option_pair_long(key='qos', value=stage.queue)

        if stage.reservation:
            executable_drms.add_option_pair_long(key='reservation', value=stage.reservation)

        # Working directory, standard output and standard error streams.

        if stage.working_directory:
            executable_drms.add_option_pair_long(key='chdir', value=stage.working_directory)

            # Write standard output and standard error streams into an
            # output directory under the working directory.
            # Create the output directory first via its absolute path,
            # then set standard error and standard output relative to it.

            output_directory_path = os.path.join(stage.working_directory, output_directory_name)

            if not os.path.isdir(output_directory_path):
                try:
                    os.makedirs(output_directory_path)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise exception

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

        process_identifier_list: list[str] = list()

        for executable_name in executable.dependencies:
            process_slurm_list = process_slurm_adaptor.select_all_by_job_name(name=executable_name)
            if len(process_slurm_list):
                # This bsf.process.Executable has been submitted at least once before.
                # For the moment, set the dependency on the last submission.
                process_identifier_list.append(process_slurm_list[-1].job_id)
            else:
                module_logger.debug(
                    'While submitting Executable.name %r, Executable.name %r that it depends on, '
                    'has not been submitted before.',
                    executable.name,
                    executable_name)

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

        # Finally, submit this command if requested.

        if executable.submit and drms_submit:
            exception_str_list = executable_drms.run()

            if exception_str_list:
                exception_str_list.append('Command list representation: ' + repr(executable_drms.command_list()))
                raise Exception('\n'.join(exception_str_list))

        # Copy the SLURM command line to the Bash script.

        output_list.append(executable_drms.command_str() + '\n')
        output_list.append('\n')

        # Regardless of an actual bsf.process.Executable submission, UPDATE it in or INSERT it into the SQLite database.

        job_submission = job_submission_adaptor.select_by_name(name=executable.name)
        if job_submission:
            job_submission.command = executable.command_str()
            job_submission_adaptor.update(object_instance=job_submission)
        else:
            job_submission = JobSubmission(
                executable_id=0,
                name=executable.name,
                command=executable.command_str())
            job_submission_adaptor.insert(object_instance=job_submission)

        # Only store a ProcessSLURM, if a bsf.process.Executable has been submitted into SLURM.

        if executable.process_identifier:
            process_slurm = process_slurm_adaptor.select_by_job_id(job_id=executable.process_identifier)
            if not process_slurm:
                process_slurm = ProcessSLURM(job_id=executable.process_identifier, job_name=executable.name)
                process_slurm_adaptor.insert(object_instance=process_slurm)

        # The commit statement should affect both insert statements above.
        database_connection.commit()

    if write_script:
        with open(
                file=os.path.join(stage.working_directory, 'bsfpython_slurm_' + stage.name + '.bash'),
                mode='wt') as output_text_io:
            output_text_io.writelines(output_list)

    return


def check_state(stage: Stage) -> None:
    """Check the state of each :py:class:`bsf.process.Executable` object of a :py:class:`bsf.analysis.Stage` object.

    :param stage: A :py:class:`bsf.analysis.Stage` object.
    :type stage: Stage
    """

    def check_state_stdout(
            _stdout_handle: TextIO,
            _thread_lock: Lock,
            _process_slurm_adaptor: ProcessSLURMAdaptor) -> None:
        """Thread callable to process the SLURM :manpage:`sacct(1)` :literal:`STDOUT` stream.

        :param _stdout_handle: The :literal:`STDOUT` or :literal:`STDERR` file handle.
        :type _stdout_handle: TextIO
        :param _thread_lock: A Python :py:class:`threading.Lock` object.
        :type _thread_lock: Lock
        :param _process_slurm_adaptor: A :py:class:`bsf.drms.slurm.ProcessSLURMAdaptor` object.
        :type _process_slurm_adaptor: ProcessSLURMAdaptor
        """
        module_logger.debug("Started Runner 'STDOUT' processor.")

        dict_reader = DictReader(f=_stdout_handle, delimiter='|')

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
                requested_tres=row_dict['ReqTRES'],
                tres_usage_in_average=row_dict['TRESUsageInAve'],
                tres_usage_in_maximum=row_dict['TRESUsageInMax'],
                tres_usage_in_maximum_node=row_dict['TRESUsageInMaxNode'],
                tres_usage_in_maximum_task=row_dict['TRESUsageInMaxTask'],
                tres_usage_in_minimum=row_dict['TRESUsageInMin'],
                tres_usage_in_minimum_node=row_dict['TRESUsageInMinNode'],
                tres_usage_in_minimum_task=row_dict['TRESUsageInMinTask'],
                tres_usage_in_total=row_dict['TRESUsageInTot'],
                tres_usage_out_maximum=row_dict['TRESUsageOutMax'],
                tres_usage_out_maximum_node=row_dict['TRESUsageOutMaxNode'],
                tres_usage_out_maximum_task=row_dict['TRESUsageOutMaxTask'],
                tres_usage_out_average=row_dict['TRESUsageOutAve'],
                tres_usage_out_total=row_dict['TRESUsageOutTot'],
            )

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

        module_logger.debug("Received EOF on 'STDOUT' pipe.")

        # Commit changes to the database and explicitly disconnect so that other threads have access.
        _process_slurm_adaptor.commit()
        _process_slurm_adaptor.disconnect()

        return

    # Open or create a database.
    database_connection = DatabaseConnection(
        file_path=os.path.join(stage.working_directory, database_file_name))
    process_slurm_adaptor = ProcessSLURMAdaptor(database_connection=database_connection)

    process_slurm_list: list[ProcessSLURM] = list()

    # Get all Processes from the database that have no state set.
    process_slurm_list.extend(process_slurm_adaptor.select_all_by_state(state=None))
    # Get all processes from the database that have a state of 'PENDING' or 'RUNNING' as those need checking.
    process_slurm_list.extend(process_slurm_adaptor.select_all_by_states(state_list=['PENDING', 'RUNNING']))

    # Explicitly disconnect here so that the Thread can access the database.

    process_slurm_adaptor.disconnect()

    # Since the SLURM sacct --jobs option just requires JobIDs, remove ProcessSLURM objects from the list,
    # if their job_id represents a job step. JobID.JobStep, i.e., JobID.batch, JobID.external

    process_slurm_list = [x for x in process_slurm_list if x.job_id.find('.') == -1]

    # Return if no ProcessSLURM objects need updating.
    if not process_slurm_list:
        return

    module_logger.debug('Number of ProcessSLURM objects to check: %d', len(process_slurm_list))

    executable_drms = Executable(
        name='sacct',
        program='sacct',
        stdout=StandardOutputStream(
            thread_callable=check_state_stdout,
            thread_kwargs={'_process_slurm_adaptor': process_slurm_adaptor}))
    executable_drms.add_option_long(key='jobs', value=','.join(map(lambda x: x.job_id, process_slurm_list)))
    executable_drms.add_switch_long(key='long')
    executable_drms.add_switch_long(key='parsable')

    exception_str_list = executable_drms.run()

    if exception_str_list:
        exception_str_list.append('Command list representation: ' + repr(executable_drms.command_list()))
        raise Exception('\n'.join(exception_str_list))

    return
