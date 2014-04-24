"""Bio.BSF.Database

A package that centralises (SQLite) Database access.
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

import sqlite3
import string


class DatabaseConnection(object):
    """The BSF Database Connection class wraps the sqlite3.Connection class.

    :ivar file_path: File path
    :type file_path: str, unicode
    # :ivar connection: Connection to an sqlite3 database
    # :type connection: sqlite3.Connection
    """
    def __init__(self, file_path):
        """Initialise a Database Connection object and connect to the sqlite3 database behind.

        :param self: Database Connection
        :type self: DatabaseConnection
        :param file_path: File path
        :type file_path: str, unicode
        :return: Nothing
        :rtype: None
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        self.connection = sqlite3.connect(database=self.file_path)

    def create_schema(self):
        """Create the database schema.

        :param self: BSF Database Connection
        :type self: DatabaseConnection
        :return: Nothing
        :rtype: None
        """

        executable_adaptor = JobSubmissionAdaptor(database_connection=self)
        executable_adaptor.create_table()

        process_slurm_adaptor = ProcessSLURMAdaptor(database_connection=self)
        process_slurm_adaptor.create_table()

        process_sge_adaptor = ProcessSGEAdaptor(database_connection=self)
        process_sge_adaptor.create_table()

        self.connection.commit()


class DatabaseAdaptor(object):
    # Define columns as Python list of Python list objects with
    # SQL column name strings and SQL column constraint strings.
    # This list should be overridden in sub-classes.

    def __init__(self, database_connection):
        """Initialise a Database Adaptor object.

        :param self: Database Adaptor
        :type self: DatabaseAdaptor
        :param database_connection: Database Connection
        :type database_connection: DatabaseConnection
        :return: Nothing
        :rtype: None
        """

        self.table_name = str()

        self.column_definition = list()

        self.table_constraint = list()

        if database_connection:
            self.database_connection = database_connection
        else:
            self.database_connection = DatabaseConnection(file_path=None)

        self.connection = sqlite3.connect(database=self.database_connection.file_path)

    def expression_column_result(self):
        """Build an SQL expression for column results typically used in SELECT statements.

        :param self: DatabaseAdaptor
        :type self: DatabaseAdaptor
        :return: Column result columns expression string
        :rtype: str
        """

        return string.join(map(lambda x: x[0], self.column_definition), sep=", ")

    def expression_column_definition(self):
        """Build an SQL expression for column definitions typically used in CREATE TABLE statements.

        :param self: DatabaseAdaptor
        :type self: DatabaseAdaptor
        :return: Column definition expression string
        :rtype: str
        """

        return string.join(map(lambda x: '{} {}'.format(x[0], x[1]), self.column_definition), sep=", ")

    def expression_column_insert(self):
        """Build an SQL expression for columns and values typically used in INSERT statements.

        :param self: DatabaseAdaptor
        :type self: DatabaseAdaptor
        :return: Column definition expression string
        :rtype: str
        """

        return string.join(map(lambda x: '?', self.column_definition), sep=", ")

    def statement_create_table(self):

        return "CREATE TABLE {!r} ({})".format(self.table_name, self.expression_column_definition())

    def statement_select(self, where_clause=None, group_clause=None, having_clause=None):

        statement = str()
        statement += "SELECT {} FROM {!r}".format(self.expression_column_result(), self.table_name)

        if where_clause:
            statement += " WHERE "
            statement += where_clause

        if group_clause:
            statement += " GROUP BY "
            statement += group_clause

        if having_clause:
            statement += " HAVING "
            statement += having_clause

        return statement

    def create_table(self):

        self.database_connection.connection.cursor().execute(self.statement_create_table())
        self.database_connection.connection.commit()

    def _objects_from_statement(self, statement, parameters=None):
        """BSF Object-specific function to turn results of a SELECT statement into BSF Objects

        This needs to be overridden in the corresponding sub-class.
        :param statement: Complete SQL SELECT statement
        :type statement: str
        :param parameters: Python list of Python str (parameter) objects or None
        :type parameters: list
        :return: Python list of BSF Objects
        :rtype: list
        """

        return []

    def select_all(self):

        statement = self.statement_select()

        return self._objects_from_statement(statement=statement)

    def store(self, data_object):

        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        value_list = list()
        for name in map(lambda x: x[0], self.column_definition):
            value_list.append(data_object.__getattribute__(name))

        cursor = self.database_connection.connection.cursor()
        cursor.execute("INSERT INTO {!r} ({}) VALUES ({})".format(
            self.table_name,
            self.expression_column_result(),
            self.expression_column_insert()),
            value_list)


class JobSubmission(object):
    """The BSF Job Submission class models one job submitted into the DRMS.

    This class is equivalent to the Executable and Command classes, but much less complex.
    Command lines are stored as submitted and not broken down into sub-commands, options and arguments.

    """

    def __init__(self, executable_id=0, name=None, command=None):
        """Initialise a BSF Job Submission object.

        :param self: BSF Job Submission object
        :type self: JobSubmission
        :param executable_id: Primary key
        :type executable_id: int
        :param name: Executable name
        :type name: str
        :param command: Command line
        :type command: str
        :return: Nothing
        :rtype: None
        """
        self.executable_id = executable_id
        self.name = name
        self.command = command


class JobSubmissionAdaptor(DatabaseAdaptor):
    # Define columns as Python list of Python list objects with
    # SQL column name strings and SQL column constraint strings.

    def __init__(self, database_connection):
        """Initialise a BSF Job Submission Adaptor object.

        :param self: BSF Job Submission Adaptor object.
        :type self: JobSubmissionAdaptor
        :param database_connection: BSF Database Connection
        :type database_connection: DatabaseConnection
        :return: Nothing
        :rtype: None
        """

        super(JobSubmissionAdaptor, self).__init__(database_connection=database_connection)

        self.table_name = 'executable'

        self.column_definition = [
            # Primary key
            ['executable_id', 'INTEGER PRIMARY KEY AUTOINCREMENT'],
            # Name
            ['name', 'TEXT UNIQUE'],
            # Command as submitted into the DRMS
            ['command', 'TEXT']
        ]

        self.table_constraint = []

    def _objects_from_statement(self, statement, parameters=None):
        """BSF JobSubmissionAdaptor-specific function to turn results of a SELECT statement into
        BSF JobSubmission objects.

        :param statement: Complete SQL SELECT statement
        :type statement: str
        :param parameters: Python list of Python str (parameter) objects or None
        :type parameters: list
        :return: Python list of BSF Objects
        :rtype: list
        """

        object_list = list()

        cursor = self.database_connection.connection.cursor()

        if parameters:
            cursor.execute(statement, parameters)
        else:
            cursor.execute(statement)

        for row in cursor.fetchall():
            object_instance = JobSubmission()
            object_list.append(object_instance)
            i = 0
            for name in map(lambda x: x[0], self.column_definition):
                object_instance.__setattr__(name, row[i])
                i += 1

        return object_list

    def select_all_by_name(self, name):
        """Select all BSF Job Submission objects by name.

        :param self: BSF Job Submission Adaptor
        :type self: JobSubmissionAdaptor
        :param name: Name
        :type name: str
        :return: Python list of BSF JobSubmission objects
        :rtype: list
        """
        parameters = list()

        statement = self.statement_select(where_clause='name = ?')
        parameters.append(name)

        return self._objects_from_statement(statement=statement, parameters=parameters)


class ProcessSLURM(object):
    def __init__(self, process_slurm_id=None, job_id=None, job_name=None, partition=None,
                 max_vm_size=None, max_vm_size_node=None, max_vm_size_task=None, average_vm_size=None,
                 max_rss=None, max_rss_node=None, max_rss_task=None, average_rss=None,
                 max_pages=None, max_pages_node=None, max_pages_task=None, average_pages=None,
                 min_cpu=None, min_cpu_node=None, min_cpu_task=None, average_cpu=None,
                 number_tasks=None, allocated_cpus=None, elapsed=None, state=None, exit_code=None,
                 average_cpu_frequency=None, requested_cpu_frequency=None, requested_memory=None,
                 consumed_energy=None,
                 max_disk_read=None, max_disk_read_node=None, max_disk_read_task=None, average_disk_read=None,
                 max_disk_write=None, max_disk_write_node=None, max_disk_write_task=None, average_disk_write=None):
        """Initialise a BSF ProcessSLURM object.

        :param process_slurm_id:
        :param job_id:
        :param job_name:
        :param partition:
        :param max_vm_size:
        :param max_vm_size_node:
        :param max_vm_size_task:
        :param average_vm_size:
        :param max_rss:
        :param max_rss_node:
        :param max_rss_task:
        :param average_rss:
        :param max_pages:
        :param max_pages_node:
        :param max_pages_task:
        :param average_pages:
        :param min_cpu:
        :param min_cpu_node:
        :param min_cpu_task:
        :param average_cpu:
        :param number_tasks:
        :param allocated_cpus:
        :param elapsed:
        :param state:
        :param exit_code:
        :param average_cpu_frequency:
        :param requested_cpu_frequency:
        :param requested_memory:
        :param consumed_energy:
        :param max_disk_read:
        :param max_disk_read_node:
        :param max_disk_read_task:
        :param average_disk_read:
        :param max_disk_write:
        :param max_disk_write_node:
        :param max_disk_write_task:
        :param average_disk_write:
        """
        self.process_slurm_id = process_slurm_id
        self.job_id = job_id
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
        self.requested_cpu_frequency = requested_cpu_frequency
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


class ProcessSLURMAdaptor(DatabaseAdaptor):
    # Define columns as Python list of Python list objects with
    # SQL column name strings and SQL column constraint strings.
    # Fields from sacct --user mschuster --parsable --long

    def __init__(self, database_connection):

        super(ProcessSLURMAdaptor, self).__init__(database_connection=database_connection)

        self.table_name = 'process_slurm'

        self.column_definition = [
            # Primary key
            ['process_slurm_id', 'INTEGER PRIMARY KEY AUTOINCREMENT'],
            # JobID
            # The number of the job or job step.  It is in the form: job.jobstep.
            ['job_id', 'TEXT'],
            # JobName
            # The name of the job or job step.
            ['job_name', 'TEXT'],
            # Partition
            # Identifies the partition on which the job ran.
            ['partition', 'TEXT'],
            # MaxVMSize
            # Maximum Virtual Memory size of all tasks in job.
            ['max_vm_size', 'TEXT'],
            # MaxVMSizeNode
            # The node on which the maxvmsize occurred.
            ['max_vm_size_node', 'TEXT'],
            # MaxVMSizeTask
            # The task ID where the maxvmsize occurred.
            ['max_vm_size_task', 'TEXT'],
            # AveVMSize
            # Average Virtual Memory size of all tasks in job.
            ['average_vm_size', 'TEXT'],
            # MaxRSS
            # Maximum resident set size of all tasks in job.
            ['max_rss', 'TEXT'],
            # MaxRSSNode
            # The node on which the maxrss occurred.
            ['max_rss_node', 'TEXT'],
            # MaxRSSTask
            # The task ID where the maxrss occurred.
            ['max_rss_task', 'TEXT'],
            # AveRSS
            # Average resident set size of all tasks in job.
            ['average_rss', 'TEXT'],
            # MaxPages
            # Maximum number of page faults of all tasks in job.
            ['max_pages', 'TEXT'],
            # MaxPagesNode
            # The node on which the maxpages occurred.
            ['max_pages_node', 'TEXT'],
            # MaxPagesTask
            # The task ID where the maxpages occurred.
            ['max_pages_task', 'TEXT'],
            # AvePages
            # Average number of page faults of all tasks in job.
            ['average_pages', 'TEXT'],
            # MinCPU
            # Minimum (system + user) CPU time of all tasks in job.
            ['min_cpu', 'TEXT'],
            # MinCPUNode
            # The node on which the mincpu occurred.
            ['min_cpu_node', 'TEXT'],
            # MinCPUTask
            # The task ID where the mincpu occurred.
            ['min_cpu_task', 'TEXT'],
            # AveCPU
            # Average (system + user) CPU time of all tasks in job.
            ['average_cpu', 'TEXT'],
            # NTasks
            # Total number of tasks in a job or step.
            ['number_tasks', 'TEXT'],
            # AllocCPUS
            # Count of allocated CPUs
            ['allocated_cpus', 'TEXT'],
            # Elapsed
            # The jobs elapsed time.
            ['elapsed', 'TEXT'],
            # State
            # Displays the job status, or state.
            # Output can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or NODE_FAIL.
            ['state', 'TEXT'],
            # ExitCode
            # The exit code returned by the job script or salloc, typically as set by the exit() function.
            # Following the colon is the signal that caused the process to  terminate if it was terminated by a signal.
            ['exit_code', 'TEXT'],
            # AveCPUFreq
            # Average weighted CPU frequency of all tasks in job, in kHz.
            ['average_cpu_frequency', 'TEXT'],
            # ReqCPUFreq
            # Requested CPU frequency for the step, in kHz.
            ['requested_cpu_frequency', 'TEXT'],
            # ReqMem
            # Minimum required memory for the job, in MB.
            ['requested_memory', 'TEXT'],
            # ConsumedEnergy
            # Total energy consumed by all tasks in job, in joules.
            ['consumed_energy', 'TEXT'],
            # MaxDiskRead
            # Maximum number of bytes read by all tasks in job.
            ['max_disk_read', 'TEXT'],
            # MaxDiskReadNode
            # The node on which the maxdiskread occurred.
            ['max_disk_read_node', 'TEXT'],
            # MaxDiskReadTask
            # The task ID where the maxdiskread occurred.
            ['max_disk_read_task', 'TEXT'],
            # AveDiskRead
            # Average number of bytes read by all tasks in job.
            ['average_disk_read', 'TEXT'],
            # MaxDiskWrite
            # Maximum number of bytes written by all tasks in job.
            ['max_disk_write', 'TEXT'],
            # MaxDiskWriteNode
            # The node on which the maxdiskwrite occurred.
            ['max_disk_write_node', 'TEXT'],
            # MaxDiskWriteTask
            # The task ID where the maxdiskwrite occurred.
            ['max_disk_write_task', 'TEXT'],
            # AveDiskWrite
            # Average number of bytes written by all tasks in job.
            ['average_disk_write', 'TEXT']
        ]

        self.table_constraint = []

    def _objects_from_statement(self, statement, parameters=None):
        """BSF ProcessSLURMAdaptor-specific function to turn results of a SELECT statement into
        BSF ProcessSLURM objects.

        :param statement: Complete SQL SELECT statement
        :type statement: str
        :param parameters: Python list of Python str (parameter) objects or None
        :type parameters: list
        :return: Python list of BSF Objects
        :rtype: list
        """

        object_list = list()

        cursor = self.database_connection.connection.cursor()

        if parameters:
            cursor.execute(statement, parameters)
        else:
            cursor.execute(statement)

        for row in cursor.fetchall():
            object_instance = ProcessSLURM()
            object_list.append(object_instance)
            i = 0
            for name in map(lambda x: x[0], self.column_definition):
                object_instance.__setattr__(name, row[i])
                i += 1

        return object_list

    def select_all_by_name(self, name):
        """Select all BSF Job Submission objects by name.

        :param self: BSF Job Submission Adaptor
        :type self: JobSubmissionAdaptor
        :param name: Name
        :type name: str
        :return: Python list of BSF JobSubmission objects
        :rtype: list
        """
        parameters = list()

        statement = self.statement_select(where_clause='job_name = ?')
        parameters.append(name)

        return self._objects_from_statement(statement=statement, parameters=parameters)


class ProcessSGE(object):
    def __init__(self, process_sge_id=None, qname=None, hostname=None, sge_group=None, owner=None, job_name=None,
                 job_number=None, account=None, priority=None, submission_time=None, start_date=None, end_time=None,
                 failed=None, exit_status=None, ru_wallclock=None, project=None, department=None, granted_pe=None,
                 slots=None, task_number=None, cpu=None, mem=None, io=None, category=None, iow=None, pe_taskid=None,
                 maxvmem=None, arid=None):
        """Initialise a BSF ProcessSGE object.

        :param process_sge_id:
        :param qname:
        :param hostname:
        :param sge_group:
        :param owner:
        :param job_name:
        :param job_number:
        :param account:
        :param priority:
        :param submission_time:
        :param start_date:
        :param end_time:
        :param failed:
        :param exit_status:
        :param ru_wallclock:
        :param project:
        :param department:
        :param granted_pe:
        :param slots:
        :param task_number:
        :param cpu:
        :param mem:
        :param io:
        :param category:
        :param iow:
        :param pe_taskid:
        :param maxvmem:
        :param arid:
        """
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


class ProcessSGEAdaptor(DatabaseAdaptor):
    # For SGE accounting.

    # Define columns as Python list of Python list objects with
    # SQL column name strings and SQL column constraint strings.

    def __init__(self, database_connection):

        super(ProcessSGEAdaptor, self).__init__(database_connection=database_connection)

        self.table_name = 'process_sge'

        self.column_definition = [
            # Primary key
            ['process_sge_id', 'INTEGER PRIMARY KEY AUTOINCREMENT'],
            # qname
            # Name of the cluster queue in which the job has run.
            ['qname', 'TEXT'],
            # hostname
            # Name of the execution host.
            ['hostname', 'TEXT'],
            # group
            # The effective group id of the job owner when executing the job.
            # Since group is a reserved word in SQL this had to be renamed to sge_group.
            ['sge_group', 'TEXT'],
            # owner
            # Owner of the Grid Engine job.
            ['owner', 'TEXT'],
            # job_name
            # Job name.
            ['job_name', 'TEXT'],
            # job_number
            # Job identifier (job number).
            ['job_number', 'TEXT'],
            # account
            # An account string as specified by the qsub(1) or qalter(1) -A option.
            ['account', 'TEXT'],
            # priority
            # Priority value assigned to the job, corresponding to the priority parameter in the queue configuration
            # (see queue_conf(5)).
            ['priority', 'TEXT'],
            # submission_time
            # Submission time.
            ['submission_time', 'TEXT'],
            # start_time
            # Start time.
            ['start_date', 'TEXT'],
            # end_time
            # End time.
            ['end_time', 'TEXT'],
            # failed
            # Indicates  the  problem  which  occurred in case a job could not be started on the execution host
            # (e.g. because the owner of the job did not have a valid account on that machine).
            # If Grid Engine tries to start a job multiple times, this may lead to multiple entries in the
            # reporting file corresponding to the same job ID.
            ['failed', 'TEXT'],
            # exit_status
            # Exit status of the job script (or Grid Engine-specific status in case of certain error conditions).
            # The exit status is determined by following the normal shell conventions. If the command terminates
            # normally, the value of the command is its exit status. However, in the case that the command exits
            # abnormally, a value of 0200 (octal), 128 (decimal) is added to the value of the command to make up
            # the exit status.
            #
            # For example: If a job dies through signal 9 (SIGKILL) then the exit status becomes 128 + 9 = 137.
            ['exit_status', 'TEXT'],
            # ru_wallclock
            # Difference between end_time and start_time (see above), except that if the job fails, it is zero.
            ['ru_wallclock', 'TEXT'],
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
            # Depending on the operating system where the job was executed, some of the fields may be 0.
            #
            # project
            # The project which was assigned to the job.
            ['project', 'TEXT'],
            # department
            # The department which was assigned to the job.
            ['department', 'TEXT'],
            # granted_pe
            # The parallel environment which was selected for the job.
            ['granted_pe', 'TEXT'],
            # slots
            # The number of slots which were dispatched to the job by the scheduler.
            ['slots', 'TEXT'],
            # task_number
            # Array job task index number.
            ['task_number', 'TEXT'],
            # cpu
            # The CPU time usage in seconds.
            # The value may be affected by the ACCT_RESERVED_USAGE execd parameter (see sge_conf(5)).
            ['cpu', 'TEXT'],
            # mem
            # The integral memory usage in Gbytes seconds.
            # The value may be affected by the ACCT_RESERVED_USAGE execd parameter (see sge_conf(5)).
            ['mem', 'TEXT'],
            # io
            # The amount of data transferred in input/output operations (if available, otherwise 0).
            ['io', 'TEXT'],
            # category
            # A string specifying the job category.
            # This contains a space-separated pseudo options list for the sub, with components as follows:
            #
            #   -U user_list
            #       An owner/group ACL list composed from host_conf(5), sge_pe(5),
            #       And queue_conf(5) user_lists/xuser_lists entries.
            #       Entries from sge_conf(5) are not considered since they can
            #       only cause a job to be accepted/rejected at submit time.
            #       Omitted if there are no such configuration entries.
            #
            #   -P project_list
            #       Like -U, but for project/xproject entries.
            #
            #   -u owner
            #       The owner's user name, if it was referenced in any RQS (see sge_resource_quota(5)).
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
            ['category', 'TEXT'],
            # iow
            # The input/output wait time in seconds (if available, otherwise 0).
            ['iow', 'TEXT'],
            # pe_taskid
            # If this identifier is set, the task was part of a parallel job, and was passed to Grid Engine
            # via the qrsh -inherit interface.
            ['pe_taskid', 'TEXT'],
            # maxvmem
            # The maximum vmem size in bytes.
            # The value may be affected by the ACCT_RESERVED_USAGE execd parameter (see sge_conf(5)).
            ['maxvmem', 'TEXT'],
            # arid
            # Advance reservation identifier. If the job used the resources of an advance reservation,
            # then this field contains a positive integer identifier; otherwise the value is "0" .
            ['arid', 'TEXT']
        ]

        self.table_constraint = []

    def _objects_from_statement(self, statement, parameters=None):
        """BSF ProcessSLURMAdaptor-specific function to turn results of a SELECT statement into
        BSF ProcessSLURM objects.

        :param statement: Complete SQL SELECT statement
        :type statement: str
        :param parameters: Python list of Python str (parameter) objects or None
        :type parameters: list
        :return: Python list of BSF Objects
        :rtype: list
        """

        object_list = list()

        cursor = self.database_connection.connection.cursor()

        if parameters:
            cursor.execute(statement, parameters)
        else:
            cursor.execute(statement)

        for row in cursor.fetchall():
            object_instance = ProcessSGE()
            object_list.append(object_instance)
            i = 0
            for name in map(lambda x: x[0], self.column_definition):
                object_instance.__setattr__(name, row[i])
                i += 1

        return object_list
