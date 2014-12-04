"""bsf.database

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
    """The C{DatabaseConnection} class encapsulates the C{sqlite3.Connection} class.

    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar connection: Connection to an sqlite3 database
    @type connection: sqlite3.Connection
    """

    def __init__(self, file_path):
        """Initialise a C{DatabaseConnection} object and connect to the sqlite3 database behind.

        @param file_path: File path
        @type file_path: str | unicode
        """

        self.file_path = file_path

        self.connection = sqlite3.connect(database=self.file_path)

    def create_schema(self):
        """Create and commit the database schema.
        """

        job_submission_adaptor = JobSubmissionAdaptor(database_connection=self)
        job_submission_adaptor.create_table()

        process_slurm_adaptor = ProcessSLURMAdaptor(database_connection=self)
        process_slurm_adaptor.create_table()

        process_sge_adaptor = ProcessSGEAdaptor(database_connection=self)
        process_sge_adaptor.create_table()

        self.connection.commit()


class DatabaseAdaptor(object):
    """The C{DatabaseAdaptor} class represents as a super-class for object-specific table adaptors.

    Instance variables should be overridden in sub-classes.
    @ivar table_name: SQL database table name
    @type table_name: str
    @ivar column_definition: Python C{list} of Python C{list} objects with
        Python C{str} (SQL column name) and Python C{str} (SQL column constraint) objects.
    @type column_definition: list
    @ivar table_constraint: SQL table constraint expression
    @type table_constraint: list
    @ivar database_connection: database Connection
    @type database_connection: DatabaseConnection
    @ivar connection: SQLite connection
    @type connection: sqlite3.Connection
    """

    def __init__(self, database_connection, table_name=None, column_definition=None, table_constraint=None,
                 connection=None):
        """Initialise a C{DatabaseAdaptor} object.

        @param database_connection: C{DatabaseConnection}
        @type database_connection: DatabaseConnection
        @param table_name: SQL database table name
        @type table_name: str
        @param column_definition: Python C{list} of Python C{list} objects with
            Python C{str} (SQL column name) and Pyton C{str} (SQL column constraint) objects.
        @type column_definition: list
        @param table_constraint: SQL table constraint expression
        @type table_constraint: list
        @param connection: SQLite connection
        @type connection: sqlite3.Connection
        """

        self.database_connection = database_connection

        if table_name:
            self.table_name = table_name
        else:
            self.table_name = str()

        if column_definition:
            self.column_definition = column_definition
        else:
            self.column_definition = list()

        if table_constraint:
            self.table_constraint = table_constraint
        else:
            self.table_constraint = list()

        if connection:
            self.connection = connection
        else:
            self.connection = sqlite3.connect(database=self.database_connection.file_path)

    def _get_column_name_list_with_primary(self):
        """Build a Python C{list} of SQL column names including the primary key.

        @return: Python C{list} of Python C{str} (SQL column name) objects
        @rtype: list
        """

        return map(lambda x: x[0], self.column_definition)

    def _get_column_name_list_without_primary(self):
        """Build a Python C{list} of SQL column names excluding the primary key.

        This method excludes C{PRIMARY KEY} columns with definition C{AUTOINCREMENT},
        which must not be assigned a value in C{INSERT} or C{UPDATE} statements.
        @return: Python C{list} of Python C{str} (SQL column name) objects
        @rtype: list
        """

        return map(lambda x: x[0], filter(lambda x: 'AUTOINCREMENT' not in x[1], self.column_definition))

    def _get_column_name_for_primary(self):
        """Get the SQL column name for the primary key.
        @return: Column name for primary key
        @rtype: str
        """

        name_list = map(lambda x: x[0], filter(lambda x: 'AUTOINCREMENT' in x[1], self.column_definition))
        return name_list[0]

    def _build_column_result_expression(self):
        """Build an SQL expression of column names typically used in C{SELECT} statements.

        This method simply lists all column names of the column definition.
        @return: Column result expression string
        @rtype: str
        """

        return string.join(words=self._get_column_name_list_with_primary(), sep=', ')

    def _build_column_definition_expression(self):
        """Build an SQL expression of column definitions typically used in C{CREATE TABLE} statements.

        @return: Column definition expression string
        @rtype: str
        """

        return string.join(
            words=map(lambda x: string.join(words=(x[0], x[1]), sep=' '), self.column_definition),
            sep=', ')

    def _build_column_insert_expression(self):
        """Build an SQL expression of column names typically used in C{INSERT} statements.

        This method excludes C{PRIMARY KEY} columns with definition C{AUTOINCREMENT},
        which must not be assigned a value.
        @return: Column definition expression string
        @rtype: str
        """

        return string.join(words=self._get_column_name_list_without_primary(), sep=', ')

    def _build_value_insert_expression(self):
        """Build an SQL expression of value placeholders (?) typically used in C{INSERT} statements.

        @return: Column value expression string
        @rtype: str
        """

        return string.join(
            words=map(lambda x: '?', self._get_column_name_list_without_primary()),
            sep=', ')

    def _build_column_update_expression(self):
        """Build an SQL expression of column name and value placeholder pairs
        typically used in SQL C{UPDATE} statements.

        As in C{INSERT} expressions, leave out the C{PRIMARY KEY} columns with definition C{AUTOINCREMENT}.
        @return: SQL column name and value placeholder pair expression string
        @rtype: str
        """

        return string.join(
            words=map(lambda x: '{} = ?'.format(x), self._get_column_name_list_without_primary()),
            sep=', ')

    def statement_create_table(self):
        """Build an SQL C{CREATE TABLE} statement.

        @return: SQL C{CREATE TABLE} statement
        @rtype: str
        """

        return "CREATE TABLE {!r} ({})".format(self.table_name, self._build_column_definition_expression())

    def statement_select(self, where_clause=None, group_clause=None, having_clause=None):
        """Build an SQL C{SELECT} statement.

        @param where_clause: SQL C{WHERE} clause
        @type where_clause: str
        @param group_clause: SQL C{GROUP BY} clause
        @type group_clause: str
        @param having_clause: SQL C{HAVING} clause
        @type having_clause: str
        @return: SQL SELECT statement
        @rtype: str
        """

        statement = str()
        statement += "SELECT {} FROM {!r}".format(self._build_column_result_expression(), self.table_name)

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
        """Create the canonical table for the C{DatabaseAdaptor} sub-class.

        Before attempting to execute the SQL C{CREATE TABLE} statement, this method checks in 'sqlite_master',
        whether the table already exists in the SQLite database.
        """

        cursor = self.connection.cursor()

        statement = "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ?"
        parameters = list()
        parameters.append(self.table_name)

        cursor.execute(statement, parameters)
        rows_list = cursor.fetchmany()

        # If the list of rows is empty, create the table.

        if not len(rows_list):
            self.database_connection.connection.cursor().execute(self.statement_create_table())
            # The commit method is called in the DatabaseAdaptor.create_schema method.
            # self.database_connection.connection.commit()

    def _objects_from_statement(self, statement, parameters=None):
        """Object-specific function to turn results of a C{SELECT} statement into objects.

        This method needs overriding in the corresponding sub-class.
        @param statement: Complete SQL C{SELECT} statement
        @type statement: str
        @param parameters: Python C{list} of Python C{str} (parameter) objects or None
        @type parameters: list
        @return: Python C{list} of objects
        @rtype: list
        """

        return []

    def select_all(self):
        """Select all canonical objects corresponding to the C{DatabaseAdaptor} sub-class.

        @return: Python C{list} of objects
        @rtype: list
        """

        statement = self.statement_select()

        return self._objects_from_statement(statement=statement)

    def insert(self, data_object):
        """Insert a canonical object corresponding to the C{DatabaseAdaptor} sub-class.

        @param data_object: Object
        @type data_object: object
        """

        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        value_list = map(lambda x: data_object.__getattribute__(x), self._get_column_name_list_without_primary())

        cursor = self.database_connection.connection.cursor()
        try:
            cursor.execute(
                "INSERT INTO {!r} ({}) VALUES ({})".format(
                    self.table_name,
                    self._build_column_insert_expression(),
                    self._build_value_insert_expression()),
                value_list)
        except sqlite3.IntegrityError:
            print "Encountered sqlite3.IntegrityError for table name '{}' on the following SQL fields:". \
                format(self.table_name)
            print "Fields: {}".format(self._build_column_insert_expression())
            print "Values: {!r}".format(value_list)
        except sqlite3.OperationalError:
            print "Encountered sqlite3.OperationalError for table name '{}' on the following SQL fields:". \
                format(self.table_name)
            print "Fields: {}".format(self._build_column_insert_expression())
            print "Values: {!r}".format(value_list)

        # Update the canonical attribute containing the primary key with the last row identifier.
        last_row_identifier = cursor.lastrowid
        column_name = self.table_name + '_id'
        if last_row_identifier and hasattr(data_object, column_name):
            data_object.__setattr__(column_name, last_row_identifier)

    def update(self, data_object):
        """Update a canonical object corresponding to the C{DatabaseAdaptor} sub-class.

        @param data_object: Object
        @type data_object: object
        """
        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        value_list = map(lambda x: data_object.__getattribute__(x), self._get_column_name_list_without_primary())

        primary_name = self._get_column_name_for_primary()
        if primary_name:
            value_list.append(data_object.__getattribute__(primary_name))
        else:
            raise Exception("Cannot update table {!r} without primary key.".format(self.table_name))

        statement = "UPDATE {!r} SET {} WHERE {} = ?".format(
            self.table_name,
            self._build_column_update_expression(),
            primary_name)

        cursor = self.database_connection.connection.cursor()
        try:
            cursor.execute(statement, value_list)
        except sqlite3.IntegrityError:
            print "Encountered SQLite3 integrity error\n" \
                  "  SQL statement: {!r}\n" \
                  "  Values: {!r}".format(statement, value_list)


class JobSubmission(object):
    """The C{JobSubmission} class representing one process submitted into the
    Distributed Resource Management System (DRMS).

    This class is equivalent to the C{Executable} and C{Command} classes, but much less complex.
    Command lines are stored as submitted and not broken down into sub-commands, options and arguments.
    @ivar executable_id: Primary key
    @type executable_id: int
    @ivar name: C{Executable.name}
    @type name: str
    @ivar command: Command line
    @type command: str
    """

    def __init__(self, executable_id=0, name=None, command=None):
        """Initialise a C{JobSubmission} object.

        @param executable_id: Primary key
        @type executable_id: int
        @param name: C{Executable.name}
        @type name: str
        @param command: Command line
        @type command: str
        """
        self.executable_id = executable_id
        self.name = name
        self.command = command


class JobSubmissionAdaptor(DatabaseAdaptor):
    """The C{JobSubmissionAdaptor} class provides database access for the C{JobSubmission} class.
    """

    def __init__(self, database_connection):
        """Initialise a C{JobSubmissionAdaptor} object.

        @param database_connection: C{DatabaseConnection}
        @type database_connection: DatabaseConnection
        """

        super(JobSubmissionAdaptor, self).__init__(
            database_connection=database_connection,
            table_name='executable',
            column_definition=[
                # Primary key
                ['executable_id', 'INTEGER PRIMARY KEY ASC AUTOINCREMENT'],
                # Name
                ['name', 'TEXT UNIQUE'],
                # Command as submitted into the DRMS
                ['command', 'TEXT']
            ])

    def _objects_from_statement(self, statement, parameters=None):
        """C{JobSubmissionAdaptor}-specific function to turn results of a SQL C{SELECT} statement into
        C{JobSubmission} objects.

        @param statement: Complete SQL C{SELECT} statement
        @type statement: str
        @param parameters: Python C{list} of Python C{str} (parameter) objects or C{None}
        @type parameters: list
        @return: Python C{list} of objects
        @rtype: list
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

    def select_by_name(self, name):
        """Select one C{JobSubmission} object by name.

        @param name: Name
        @type name: str
        @return: C{JobSubmission} or C{None}
        @rtype: JobSubmission | None
        """
        parameters = list()

        statement = self.statement_select(where_clause='name = ?')
        parameters.append(name)

        object_list = self._objects_from_statement(statement=statement, parameters=parameters)
        object_length = len(object_list)

        if object_length > 1:
            raise Exception("SQL database returned more than one row for unique field 'name'.")
        elif object_length == 1:
            return object_list[0]
        else:
            return


class ProcessSLURM(object):
    """The C{ProcessSLURM} class models one process in the Simple Linux Utility for Resource Management (SLURM)
    Distributed Resource Management System (DRMS).

    The instance variable names result from the SLURM command C{sacct --parsable --long}
    @ivar process_slurm_id: Primary key
    @type process_slurm_id: int
    @ivar job_id: The number of the job or job step. It is in the form: job.jobstep
    @type job_id: str
    @ivar job_name: The name of the job or job step
    @type job_name: str
    @ivar partition: Identifies the partition on which the job ran
    @type partition: str
    @ivar max_vm_size: Maximum virtual memory size of all tasks in job
    @type max_vm_size: str
    @ivar max_vm_size_node: The node on which the maximum virtual memory size occurred
    @type max_vm_size_node: str
    @ivar max_vm_size_task: The task identifier where the maximum virtual memory size occurred
    @type max_vm_size_task: str
    @ivar average_vm_size: Average virtual memory size of all tasks in job
    @type average_vm_size: str
    @ivar max_rss: Maximum resident set size of all tasks in job
    @type max_rss: str
    @ivar max_rss_node: The node on which the maximum resident set size occurred
    @type max_rss_node: str
    @ivar max_rss_task: The task identifier where the maximum resident set size occurred
    @type max_rss_task: str
    @ivar average_rss: Average resident set size of all tasks in job
    @type average_rss: str
    @ivar max_pages: Maximum number of page faults of all tasks in job
    @type max_pages: str
    @ivar max_pages_node: The node on which the maximum number of page faults occurred
    @type max_pages_node: str
    @ivar max_pages_task: The task identifier where the maximum number of page faults occurred
    @type max_pages_task: str
    @ivar average_pages: Average number of page faults of all tasks in job
    @type average_pages: str
    @ivar min_cpu: Minimum (system + user) CPU time of all tasks in job
    @type min_cpu: str
    @ivar min_cpu_node: The node on which the minimum CPU time occurred
    @type min_cpu_node: str
    @ivar min_cpu_task: The task identifier where the minimum CPU time occurred
    @type min_cpu_task: str
    @ivar average_cpu: Average (system + user) CPU time of all tasks in job
    @type average_cpu: str
    @ivar number_tasks: Total number of tasks in a job or step
    @type number_tasks: str
    @ivar allocated_cpus: Count of allocated CPUs
    @type allocated_cpus: str
    @ivar elapsed: The jobs elapsed time
    @type elapsed: str
    @ivar state: Displays the job status, or state
        Value can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or NODE_FAIL
    @type state: str
    @ivar exit_code: The exit code returned by the job script or salloc, typically as set by the exit() function.
        Following the colon is the signal that caused the process to terminate if it was terminated by a signal.
    @type exit_code: str
    @ivar average_cpu_frequency: Average weighted CPU frequency of all tasks in job, in kHz
    @type average_cpu_frequency: str
    @ivar requested_cpu_frequency: Requested CPU frequency for the step, in kHz
    @type requested_cpu_frequency: str
    @ivar requested_memory: Minimum required memory for the job, in MB
    @type requested_memory: str
    @ivar consumed_energy: Total energy consumed by all tasks in job, in joules
    @type consumed_energy: str
    @ivar max_disk_read: Maximum number of bytes read by all tasks in job
    @type max_disk_read: str
    @ivar max_disk_read_node: The node on which the maximum number of bytes read occurred
    @type max_disk_read_node: str
    @ivar max_disk_read_task: The task identifier where the maximum number of bytes read occurred
    @type max_disk_read_task: str
    @ivar average_disk_read: Average number of bytes read by all tasks in job
    @type average_disk_read: str
    @ivar max_disk_write: Maximum number of bytes written by all tasks in job
    @type max_disk_write: str
    @ivar max_disk_write_node: The node on which the maximum number of bytes written occurred
    @type max_disk_write_node: str
    @ivar max_disk_write_task: The task identifier where the maximum number of bytes written occurred
    @type max_disk_write_task: str
    @ivar average_disk_write: Average number of bytes written by all tasks in job
    @type average_disk_write: str
    """

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
        """Initialise a ProcessSLURM object.

        @param process_slurm_id:
        @type process_slurm_id: int
        @param job_id: The number of the job or job step. It is in the form: job.jobstep
        @type job_id: str
        @param job_name: The name of the job or job step
        @type job_name: str
        @param partition: Identifies the partition on which the job ran
        @type partition: str
        @param max_vm_size: Maximum virtual memory size of all tasks in job
        @type max_vm_size: str
        @param max_vm_size_node: The node on which the maximum virtual memory size occurred
        @type max_vm_size_node: str
        @param max_vm_size_task: The task identifier where the maximum virtual memory size occurred
        @type max_vm_size_task: str
        @param average_vm_size: Average virtual memory size of all tasks in job
        @type average_vm_size: str
        @param max_rss: Maximum resident set size of all tasks in job
        @type max_rss: str
        @param max_rss_node: The node on which the maximum resident set size occurred
        @type max_rss_node: str
        @param max_rss_task: The task identifier where the maximum resident set size occurred
        @type max_rss_task: str
        @param average_rss: Average resident set size of all tasks in job
        @type average_rss: str
        @param max_pages: Maximum number of page faults of all tasks in job
        @type max_pages: str
        @param max_pages_node: The node on which the maximum number of page faults occurred
        @type max_pages_node: str
        @param max_pages_task: The task identifier where the maximum number of page faults occurred
        @type max_pages_task: str
        @param average_pages: Average number of page faults of all tasks in job
        @type average_pages: str
        @param min_cpu: Minimum (system + user) CPU time of all tasks in job
        @type min_cpu: str
        @param min_cpu_node: The node on which the minimum CPU time occurred
        @type min_cpu_node: str
        @param min_cpu_task: The task identifier where the minimum CPU time occurred
        @type min_cpu_task: str
        @param average_cpu: Average (system + user) CPU time of all tasks in job
        @type average_cpu: str
        @param number_tasks: Total number of tasks in a job or step
        @type number_tasks: str
        @param allocated_cpus: Count of allocated CPUs
        @type allocated_cpus: str
        @param elapsed: The jobs elapsed time
        @type elapsed: str
        @param state: Displays the job status, or state.
            Value can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or NODE_FAIL
        @type state: str
        @param exit_code: The exit code returned by the job script or salloc, typically as set by the exit() function.
            Following the colon is the signal that caused the process to  terminate if it was terminated by a signal.
        @type exit_code: str
        @param average_cpu_frequency: Average weighted CPU frequency of all tasks in job, in kHz
        @type average_cpu_frequency: str
        @param requested_cpu_frequency: Requested CPU frequency for the step, in kHz
        @type requested_cpu_frequency: str
        @param requested_memory: Minimum required memory for the job, in MB
        @type requested_memory: str
        @param consumed_energy: Total energy consumed by all tasks in job, in joules
        @type consumed_energy: str
        @param max_disk_read: Maximum number of bytes read by all tasks in job
        @type max_disk_read: str
        @param max_disk_read_node: The node on which the maximum number of bytes read occurred
        @type max_disk_read_node: str
        @param max_disk_read_task: The task identifier where the maximum number of bytes read occurred
        @type max_disk_read_task: str
        @param average_disk_read: Average number of bytes read by all tasks in job
        @type average_disk_read: str
        @param max_disk_write: Maximum number of bytes written by all tasks in job
        @type max_disk_write: str
        @param max_disk_write_node: The node on which the maximum number of bytes written occurred
        @type max_disk_write_node: str
        @param max_disk_write_task: The task identifier where the maximum number of bytes written occurred
        @type max_disk_write_task: str
        @param average_disk_write: Average number of bytes written by all tasks in job
        @type average_disk_write: str
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
    """ProcessSLURMAdaptor class providing database access for the ProcessSLURM class.

    The SQL column names result from SLURM command sacct --parsable --long
    """

    def __init__(self, database_connection):
        """Initialise a ProcessSLURMAdaptor object.

        @param database_connection: C{DatabaseConnection}
        @type database_connection: DatabaseConnection
        """

        super(ProcessSLURMAdaptor, self).__init__(
            database_connection=database_connection,
            table_name='process_slurm',
            column_definition=[
                # Primary key
                ['process_slurm_id', 'INTEGER PRIMARY KEY ASC AUTOINCREMENT'],
                # JobID
                # The number of the job or job step. It is in the form: job.jobstep.
                ['job_id', 'TEXT UNIQUE'],
                # JobName
                # The name of the job or job step.
                ['job_name', 'TEXT'],
                # Partition
                # Identifies the partition on which the job ran.
                ['partition', 'TEXT'],
                # MaxVMSize
                # Maximum virtual memory size of all tasks in job.
                ['max_vm_size', 'TEXT'],
                # MaxVMSizeNode
                # The node on which the maximum virtual memory size occurred.
                ['max_vm_size_node', 'TEXT'],
                # MaxVMSizeTask
                # The task identifier where the maximum virtual memory size occurred.
                ['max_vm_size_task', 'TEXT'],
                # AveVMSize
                # Average virtual memory size of all tasks in job.
                ['average_vm_size', 'TEXT'],
                # MaxRSS
                # Maximum resident set size of all tasks in job.
                ['max_rss', 'TEXT'],
                # MaxRSSNode
                # The node on which the maximum resident set size occurred.
                ['max_rss_node', 'TEXT'],
                # MaxRSSTask
                # The task identifier where the maximum resident set size occurred.
                ['max_rss_task', 'TEXT'],
                # AveRSS
                # Average resident set size of all tasks in job.
                ['average_rss', 'TEXT'],
                # MaxPages
                # Maximum number of page faults of all tasks in job.
                ['max_pages', 'TEXT'],
                # MaxPagesNode
                # The node on which the maximum number of page faults occurred.
                ['max_pages_node', 'TEXT'],
                # MaxPagesTask
                # The task identifier where the maximum number of page faults occurred.
                ['max_pages_task', 'TEXT'],
                # AvePages
                # Average number of page faults of all tasks in job.
                ['average_pages', 'TEXT'],
                # MinCPU
                # Minimum (system + user) CPU time of all tasks in job.
                ['min_cpu', 'TEXT'],
                # MinCPUNode
                # The node on which the minimum CPU time occurred.
                ['min_cpu_node', 'TEXT'],
                # MinCPUTask
                # The task identifier where the minimum CPU time occurred.
                ['min_cpu_task', 'TEXT'],
                # AveCPU
                # Average (system + user) CPU time of all tasks in job.
                ['average_cpu', 'TEXT'],
                # NTasks
                # Total number of tasks in a job or step.
                ['number_tasks', 'TEXT'],
                # AllocCPUS
                # Count of allocated CPUs.
                ['allocated_cpus', 'TEXT'],
                # Elapsed
                # The jobs elapsed time.
                ['elapsed', 'TEXT'],
                # State
                # Displays the job status, or state.
                # Value can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, PREEMPTED or
                # NODE_FAIL.
                ['state', 'TEXT'],
                # ExitCode
                # The exit code returned by the job script or salloc, typically as set by the exit() function.
                # Following the colon is the signal that caused the process to  terminate if it was terminated by
                # a signal.
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
                # The node on which the maximum number of bytes read occurred.
                ['max_disk_read_node', 'TEXT'],
                # MaxDiskReadTask
                # The task identifier where the maximum number of bytes read occurred.
                ['max_disk_read_task', 'TEXT'],
                # AveDiskRead
                # Average number of bytes read by all tasks in job.
                ['average_disk_read', 'TEXT'],
                # MaxDiskWrite
                # Maximum number of bytes written by all tasks in job.
                ['max_disk_write', 'TEXT'],
                # MaxDiskWriteNode
                # The node on which the maximum number of bytes written occurred.
                ['max_disk_write_node', 'TEXT'],
                # MaxDiskWriteTask
                # The task identifier where the maximum number of bytes written occurred.
                ['max_disk_write_task', 'TEXT'],
                # AveDiskWrite
                # Average number of bytes written by all tasks in job.
                ['average_disk_write', 'TEXT']
            ])

    def _objects_from_statement(self, statement, parameters=None):
        """ProcessSLURMAdaptor-specific function to turn results of a SELECT statement into
        ProcessSLURM objects.

        @param statement: Complete SQL SELECT statement
        @type statement: str
        @param parameters: Python list of Python str (parameter) objects or None
        @type parameters: list
        @return: Python list of objects
        @rtype: list
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

    def select_all_by_job_name(self, name):
        """Select all JobSubmission objects by job_name.

        The same Executable can be submitted more than once into the DRMS.
        @param name: Job name
        @type name: str
        @return: Python list of JobSubmission objects
        @rtype: list
        """
        statement = self.statement_select(where_clause='job_name = ?')
        parameters = list()
        parameters.append(name)

        return self._objects_from_statement(statement=statement, parameters=parameters)

    def select_by_job_id(self, job_id):
        """Select one JobSubmission object by job_id.

        @param job_id: Job identifier
        @type job_id: str
        @return: Python list of JobSubmission objects
        @rtype: list
        """
        statement = self.statement_select(where_clause='job_id = ?')
        parameters = list()
        parameters.append(job_id)

        object_list = self._objects_from_statement(statement=statement, parameters=parameters)
        object_length = len(object_list)

        if object_length > 1:
            raise Exception("SQL database returned more than one row for unique field 'job_id'.")
        elif object_length == 1:
            return object_list[0]
        else:
            return


class ProcessSGE(object):
    """The ProcessSLURM class models one process in the Son of Grid Engine (SGE)
    Distributed Resource Management System.

    The instance variable names result from the SGE accounting file. See man 5 accounting.
    @ivar process_sge_id: Primary key
    @type process_sge_id: int
    @ivar qname: Name of the cluster queue in which the job has run
    @type qname: str
    @ivar hostname: Name of the execution host
    @type hostname: str
    @ivar sge_group: The effective group id of the job owner when executing the job
    @type sge_group: str
    @ivar owner: Owner of the Grid Engine job
    @type owner: str
    @ivar job_name: Job name
    @type job_name: str
    @ivar job_number: Job identifier (job number)
    @type job_number: str
    @ivar account: An account string as specified by the qsub(1) or qalter(1) -A option
    @type account: str
    @ivar priority: Priority value assigned to the job, corresponding to the priority parameter in the
        queue configuration (see queue_conf(5))
    @type priority: str
    @ivar submission_time: Submission time
    @type submission_time: str
    @ivar start_date: Start time
    @type start_date: str
    @ivar end_time: End time
    @type end_time: str
    @ivar failed: Indicates the problem which occurred in case a job could not be started on the execution host
    @type failed: str
    @ivar exit_status: Exit status of the job script (or Grid Engine-specific status in case of certain error
        conditions). The exit status is determined by following the normal shell conventions. If the command
        terminates normally, the value of the command is its exit status. However, in the case that the command
        exits abnormally, a value of 0200 (octal), 128 (decimal) is added to the value of the command to make up
        the exit status.
    @type exit_status: str
    @ivar ru_wallclock: Difference between end_time and start_time (see above), except that if the job fails,
        it is zero
    @type ru_wallclock: str
    @ivar project: The department which was assigned to the job
    @type project: str
    @ivar department: The parallel environment which was selected for the job
    @type department: str
    @ivar granted_pe: The number of slots which were dispatched to the job by the scheduler
    @type granted_pe: str
    @ivar slots: The number of slots which were dispatched to the job by the scheduler
    @type slots: str
    @ivar task_number: Array job task index number
    @type task_number: str
    @ivar cpu: The CPU time usage in seconds
    @type cpu: str
    @ivar mem: The integral memory usage in Gbytes seconds
    @type mem: str
    @ivar io: The amount of data transferred in input/output operations in GB (if available, otherwise 0)
    @type io: str
    @ivar category: A string specifying the job category
    @type category: str
    @ivar iow: The input/output wait time in seconds (if available, otherwise 0)
    @type iow: str
    @ivar pe_taskid: If this identifier is set, the task was part of a parallel job, and was passed to Grid Engine
        via the qrsh -inherit interface.
    @type pe_taskid: str
    @ivar maxvmem: The maximum vmem size in bytes
    @type maxvmem: str
    @ivar arid: Advance reservation identifier
    @type arid: str
    """

    def __init__(self, process_sge_id=None, qname=None, hostname=None, sge_group=None, owner=None, job_name=None,
                 job_number=None, account=None, priority=None, submission_time=None, start_date=None, end_time=None,
                 failed=None, exit_status=None, ru_wallclock=None, project=None, department=None, granted_pe=None,
                 slots=None, task_number=None, cpu=None, mem=None, io=None, category=None, iow=None, pe_taskid=None,
                 maxvmem=None, arid=None):
        """Initialise a ProcessSGE object.

        @param process_sge_id: Primary key
        @type process_sge_id: int
        @param qname: Name of the cluster queue in which the job has run
        @type qname: str
        @param hostname: Name of the execution host
        @type hostname: str
        @param sge_group: The effective group id of the job owner when executing the job
        @type sge_group: str
        @param owner: Owner of the Grid Engine job
        @type owner: str
        @param job_name: Job name
        @type job_name: str
        @param job_number: Job identifier (job number)
        @type job_number: str
        @param account: An account string as specified by the qsub(1) or qalter(1) -A option
        @type account: str
        @param priority: Priority value assigned to the job, corresponding to the priority parameter in the
            queue configuration (see queue_conf(5))
        @type priority: str
        @param submission_time: Submission time
        @type submission_time: str
        @param start_date: Start time
        @type start_date: str
        @param end_time: End time
        @type end_time: str
        @param failed: Indicates the problem which occurred in case a job could not be started on the execution host
        @type failed: str
        @param exit_status: Exit status of the job script (or Grid Engine-specific status in case of certain error
            conditions). The exit status is determined by following the normal shell conventions. If the command
            terminates normally, the value of the command is its exit status. However, in the case that the command
            exits abnormally, a value of 0200 (octal), 128 (decimal) is added to the value of the command to make up
            the exit status.
        @type exit_status: str
        @param ru_wallclock: Difference between end_time and start_time (see above), except that if the job fails,
            it is zero
        @type ru_wallclock: str
        @param project: The department which was assigned to the job
        @type project: str
        @param department: The parallel environment which was selected for the job
        @type department: str
        @param granted_pe: The number of slots which were dispatched to the job by the scheduler
        @type granted_pe: str
        @param slots: The number of slots which were dispatched to the job by the scheduler
        @type slots: str
        @param task_number: Array job task index number
        @type task_number: str
        @param cpu: The CPU time usage in seconds
        @type cpu: str
        @param mem: The integral memory usage in Gbytes seconds
        @type mem: str
        @param io: The amount of data transferred in input/output operations in GB (if available, otherwise 0)
        @type io: str
        @param category: A string specifying the job category
        @type category: str
        @param iow: The input/output wait time in seconds (if available, otherwise 0)
        @type iow: str
        @param pe_taskid: If this identifier is set, the task was part of a parallel job, and was passed to Grid Engine
            via the qrsh -inherit interface.
        @type pe_taskid: str
        @param maxvmem: The maximum vmem size in bytes
        @type maxvmem: str
        @param arid: Advance reservation identifier
        @type arid: str
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
    """ProcessSLURMAdaptor class providing database access for the ProcessSLURM class.

    The SQL column names result from the SGE accounting file. See man 5 accounting.
    """

    def __init__(self, database_connection):
        """Initialise a ProcessSGEAdaptor object.

        @param database_connection: C{DatabaseConnection}
        @type database_connection: DatabaseConnection
        """

        super(ProcessSGEAdaptor, self).__init__(
            database_connection=database_connection,
            table_name='process_sge',
            column_definition=[
                # Primary key
                ['process_sge_id', 'INTEGER PRIMARY KEY ASC AUTOINCREMENT'],
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
                # Indicates the problem which occurred in case a job could not be started on the execution host
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
            ])

    def _objects_from_statement(self, statement, parameters=None):
        """ProcessSLURMAdaptor-specific function to turn results of a SELECT statement into
        ProcessSLURM objects.

        @param statement: Complete SQL SELECT statement
        @type statement: str
        @param parameters: Python list of Python str (parameter) objects or None
        @type parameters: list
        @return: Python list of objects
        @rtype: list
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
