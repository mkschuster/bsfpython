"""bsf.database

A package that centralises (SQLite) Database access.
"""

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

import sqlite3


class DatabaseConnection(object):
    """The C{DatabaseConnection} class encapsulates the C{sqlite3.Connection} class.

    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar _connection: C{sqlite3.Connection}
    @type _connection: sqlite3.Connection
    """

    def __init__(
            self,
            file_path=None):
        """Initialise a C{DatabaseConnection} object.

        The underlying C{sqlite3.Connection} object is instantiated only upon calling C{connect}.
        @param file_path: File path
        @type file_path: str | unicode
        """

        super(DatabaseConnection, self).__init__()

        if file_path is None:
            self.file_path = ':memory:'
        else:
            self.file_path = file_path

        self._connection = None

        return

    def __del__(self):
        """Delete a C{DatabaseConnection} object, which implies closing its C{sqlite3.Connection} first.

        @return:
        @rtype:
        """

        if self._connection is None:
            return
        else:
            return self._connection.close()

    def connect(self):
        """Connect a C{DatabaseConnection} by instantiating the underlying C{sqlite3.Connection} object.

        Just returns if the C{sqlite3.Connection} object already exists.
        @return:
        @rtype:
        """

        if self._connection is None:
            self._connection = sqlite3.connect(database=self.file_path)

        return

    def disconnect(self):
        """Disconnect a C{DatabaseConnection} by closing the underlying C{sqlite3.Connection} object and
        replacing it by C{None}.

        Just returns if the C{Connection} object does not exist.
        @return:
        @rtype:
        """

        if self._connection is not None:
            self._connection.close()
            self._connection = None

        return

    def commit(self):
        """Commit changes to the underlying C{sqlite3.Connection} object.

        @return:
        @rtype:
        """

        if self._connection is not None:
            return self._connection.commit()
        else:
            return

    def get_cursor(self):
        """Get a C{sqlite3.Cursor} object of the underlying C{sqlite3.Connection} object.

        Creates (connects) the underlying C{sqlite3.Connection} object if it does not exist already.
        @return: C{sqlite3.Cursor}
        @rtype: sqlite3.Cursor
        """

        self.connect()

        return self._connection.cursor()


class DatabaseAdaptor(object):
    """The C{DatabaseAdaptor} class represents as a super-class for object-specific table adaptors.

    Instance variables should be overridden in sub-classes.
    @ivar database_connection: C{DatabaseConnection}
    @type database_connection: DatabaseConnection
    @ivar table_name: SQL database table name
    @type table_name: str
    @ivar column_definition: Python C{list} of Python C{list} objects with
        Python C{str} (SQL column name) and Python C{str} (SQL column constraint) objects.
    @type column_definition: list[list[str]]
    @ivar table_constraint: SQL table constraint expression
    @type table_constraint: list[str]
    """

    def __init__(
            self,
            database_connection,
            object_type,
            table_name=None,
            column_definition=None,
            table_constraint=None):
        """Initialise a C{DatabaseAdaptor} object.

        @param database_connection: C{DatabaseConnection}
        @type database_connection: DatabaseConnection
        @param object_type: Object type
        @type object_type: type
        @param table_name: SQL database table name
        @type table_name: str
        @param column_definition: Python C{list} of Python C{list} objects with
            Python C{str} (SQL column name) and Python C{str} (SQL column constraint) objects.
        @type column_definition: list[list[str]]
        @param table_constraint: SQL table constraint expression
        @type table_constraint: list[str]
        """

        assert isinstance(database_connection, DatabaseConnection)
        assert isinstance(object_type, type)

        super(DatabaseAdaptor, self).__init__()

        self.database_connection = database_connection  # Can be None.
        self.object_type = object_type  # Can be None.

        if table_name is None:
            self.table_name = str()
        else:
            self.table_name = table_name

        if column_definition is None:
            self.column_definition = list()
        else:
            self.column_definition = column_definition

        if table_constraint is None:
            self.table_constraint = list()
        else:
            self.table_constraint = table_constraint

        # FIXME: Experimentally create the table in the __init__() method.
        # TODO: Should this be moved into the module that uses the Adaptor?
        # Try creating the table once the DatabaseAdaptor gets instantiated.
        # It then needs committing at some stage.
        self.create_table()

        return

    def _get_column_name_list_with_primary(self):
        """Build a Python C{list} of SQL column names including the primary key.

        @return: Python C{list} of Python C{str} (SQL column name) objects
        @rtype: list[str]
        """

        return map(lambda x: x[0], self.column_definition)

    def _get_column_name_list_without_primary(self):
        """Build a Python C{list} of SQL column names excluding the primary key.

        This method excludes I{PRIMARY KEY} columns with definition I{AUTOINCREMENT},
        which must not be assigned a value in I{INSERT} or I{UPDATE} statements.
        @return: Python C{list} of Python C{str} (SQL column name) objects
        @rtype: list[str]
        """

        return map(lambda x: x[0], filter(lambda x: 'AUTOINCREMENT' not in x[1], self.column_definition))

    def _get_column_name_for_primary(self):
        """Get the SQL column name for the primary key.

        This method returns the I{PRIMARY KEY} column with definition I{AUTOINCREMENT}.
        @return: Column name for primary key
        @rtype: str
        """

        name_list = map(lambda x: x[0], filter(lambda x: 'AUTOINCREMENT' in x[1], self.column_definition))

        list_length = len(name_list)
        if list_length < 1:
            return str()
        elif list_length == 1:
            return name_list[0]
        else:
            raise Exception("SQL column definition with more than one 'AUTOINCREMENT' attribute.")

    def _build_column_result_expression(self):
        """Build an SQL expression of column names typically used in I{SELECT} statements.

        This method simply lists all column names of the column definition.
        @return: Column result expression string
        @rtype: str
        """

        return ', '.join(self._get_column_name_list_with_primary())

    def _build_column_definition_expression(self):
        """Build an SQL expression of column definitions typically used in I{CREATE TABLE} statements.

        @return: Column definition expression string
        @rtype: str
        """

        return ', '.join(map(lambda x: ' '.join((x[0], x[1])), self.column_definition))

    def _build_column_insert_expression(self):
        """Build an SQL expression of column names typically used in I{INSERT} statements.

        This method excludes C{PRIMARY KEY} columns with definition I{AUTOINCREMENT},
        which must not be assigned a value.
        @return: Column definition expression string
        @rtype: str
        """

        return ', '.join(self._get_column_name_list_without_primary())

    def _build_value_insert_expression(self):
        """Build an SQL expression of value placeholders (?) typically used in I{INSERT} statements.

        @return: Column value expression string
        @rtype: str
        """

        return ', '.join(map(lambda x: '?', self._get_column_name_list_without_primary()))

    def _build_column_update_expression(self):
        """Build an SQL expression of column name and value placeholder pairs
        typically used in SQL I{UPDATE} statements.

        As in I{INSERT} expressions, leave out the I{PRIMARY KEY} columns with definition I{AUTOINCREMENT}.
        @return: SQL column name and value placeholder pair expression string
        @rtype: str
        """

        return ', '.join(map(lambda x: '{} = ?'.format(x), self._get_column_name_list_without_primary()))

    def connect(self):
        """Connect a C{DatabaseAdaptor} via its C{DatabaseConnection} by instantiating the underlying
        C{sqlite3.Connection} object.

        Just returns if the C{sqlite3.Connection} object already exists.
        @return:
        @rtype:
        """

        return self.database_connection.connect()

    def disconnect(self):
        """Explicitly disconnect the underlying C{sqlite3.Connection} via the C{DatabaseConnection}.

        @return:
        @rtype:
        """

        return self.database_connection.disconnect()

    def commit(self):
        """Commit changes to the underlying C{sqlite3.Connection} via the C{DatabaseConnection}.

        @return:
        @rtype:
        """

        return self.database_connection.commit()

    def get_cursor(self):
        """Get the C{sqlite3.Cursor} object of the underlying C{DatabaseConnection} and C{sqlite3.Connection} objects.

        @return: C{sqlite3.Cursor}
        @rtype: sqlite3.Cursor
        """

        return self.database_connection.get_cursor()

    def statement_create_table(self):
        """Build an SQL I{CREATE TABLE} statement.

        @return: SQL I{CREATE TABLE} statement
        @rtype: str
        """

        return "CREATE TABLE {!r} ({})".format(self.table_name, self._build_column_definition_expression())

    def statement_select(self, where_clause=None, group_clause=None, having_clause=None):
        """Build an SQL I{SELECT} statement.

        @param where_clause: SQL I{WHERE} clause
        @type where_clause: str
        @param group_clause: SQL I{GROUP BY} clause
        @type group_clause: str
        @param having_clause: SQL I{HAVING} clause
        @type having_clause: str
        @return: SQL I{SELECT} statement
        @rtype: str
        """

        statement = str()
        statement += "SELECT {} FROM {!r}".format(self._build_column_result_expression(), self.table_name)

        if where_clause is not None and len(where_clause):
            statement += " WHERE "
            statement += where_clause

        if group_clause is not None and len(group_clause):
            statement += " GROUP BY "
            statement += group_clause

        if having_clause is not None and len(having_clause):
            statement += " HAVING "
            statement += having_clause

        return statement

    def create_table(self):
        """Execute a SQL I{CREATE TABLE} statement for the canonical table of the C{DatabaseAdaptor} sub-class.

        Before attempting to execute the SQL I{CREATE TABLE} statement, this method checks in 'sqlite_master',
        whether the table already exists in the SQLite database.
        After calling, the commit() method has to be called at some stage.
        """

        cursor = self.get_cursor()

        statement = "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ?"
        parameters = list()
        parameters.append(self.table_name)

        cursor.execute(statement, parameters)
        rows_list = cursor.fetchmany()

        # If the list of rows is empty, create the table.

        if not len(rows_list):
            self.get_cursor().execute(self.statement_create_table())

        return

    def select(self, statement, parameters=None):
        """Execute a SQL I{SELECT} statement and return canonical Python C{object} instances of the
        C{DatabaseAdaptor} sub-class.

        @param statement: Complete SQL I{SELECT} statement
        @type statement: str
        @param parameters: Python C{list} of Python C{str} (parameter) objects or C{None}
        @type parameters: list[str]
        @return: Python C{list} of Python C{object} objects
        @rtype: list[object]
        """

        object_list = list()

        cursor = self.get_cursor()

        if parameters:
            cursor.execute(statement, parameters)
        else:
            cursor.execute(statement)

        # FIXME: Setting attributes directly only works for type str only!
        # TODO: Investigate automatic type mapping.
        for row in cursor.fetchall():
            object_instance = self.object_type()
            object_list.append(object_instance)
            i = 0
            for name in map(lambda x: x[0], self.column_definition):
                object_instance.__setattr__(name, row[i])
                i += 1

        return object_list

    def select_all(self):
        """Select all canonical Python C{object} instances corresponding to the C{DatabaseAdaptor} sub-class.

        @return: Python C{list} of Python C{object} objects
        @rtype: list[object]
        """

        statement = self.statement_select()

        return self.select(statement=statement)

    def select_by_identifier(self, identifier=0):
        """Select one canonical Python C{object} instance corresponding to the C{DatabaseAdaptor} sub-class
        by its primary key identifier.

        @param identifier: Primary key identifier
        @type identifier: int
        @return: Python C{object} instance
        @rtype: object
        """
        # Check if the table has a primary identifier.

        primary_key = self._get_column_name_for_primary()
        if not primary_key:
            return

        parameters = list()

        # statement = self.statement_select(where_clause='{}_id = ?'.format(self.table_name))
        statement = self.statement_select(where_clause='{} = ?'.format(primary_key))
        parameters.append(identifier)

        object_list = self.select(statement=statement, parameters=parameters)
        object_length = len(object_list)

        if object_length > 1:
            raise Exception("SQL database returned more than one row for unique field '{}_id'.".format(
                self.table_name))
        elif object_length == 1:
            return object_list[0]
        else:
            return

    def insert(self, object_instance):
        """Execute a SQL I{INSERT} statement for a canonical Python C{object} instance corresponding to the
        C{DatabaseAdaptor} sub-class.

        @param object_instance: Python C{object} object
        @type object_instance: object
        """

        assert isinstance(object_instance, self.object_type)

        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        value_list = map(lambda x: object_instance.__getattribute__(x), self._get_column_name_list_without_primary())

        cursor = self.get_cursor()
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
        if last_row_identifier and hasattr(object_instance, column_name):
            object_instance.__setattr__(column_name, last_row_identifier)

        return

    def update(self, object_instance):
        """Execute a SQL I{UPDATE} statement for a canonical Python C{object} instance corresponding to the
        C{DatabaseAdaptor} sub-class.

        @param object_instance: Python C{object} instance
        @type object_instance: object
        """

        assert isinstance(object_instance, self.object_type)
        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        value_list = map(lambda x: object_instance.__getattribute__(x), self._get_column_name_list_without_primary())

        primary_name = self._get_column_name_for_primary()
        if primary_name:
            value_list.append(object_instance.__getattribute__(primary_name))
        else:
            raise Exception("Cannot update table {!r} without primary key.".format(self.table_name))

        statement = "UPDATE {!r} SET {} WHERE {} = ?".format(
            self.table_name,
            self._build_column_update_expression(),
            primary_name)

        cursor = self.get_cursor()
        try:
            cursor.execute(statement, value_list)
        except sqlite3.IntegrityError:
            print "Encountered SQLite3 integrity error\n" \
                  "  SQL statement: {!r}\n" \
                  "  Values: {!r}".format(statement, value_list)

        return


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

    def __init__(
            self,
            executable_id=0,
            name=None,
            command=None):
        """Initialise a C{JobSubmission} object.

        @param executable_id: Primary key
        @type executable_id: int
        @param name: C{Executable.name}
        @type name: str
        @param command: Command line
        @type command: str
        """

        super(JobSubmission, self).__init__()

        self.executable_id = executable_id  # Can be 0.
        self.name = name  # Can be None.
        self.command = command  # Can be None.

        return


class JobSubmissionAdaptor(DatabaseAdaptor):
    """The C{JobSubmissionAdaptor} class provides database access for the C{JobSubmission} class.
    """

    def __init__(
            self,
            database_connection):
        """Initialise a C{JobSubmissionAdaptor} object.

        @param database_connection: C{DatabaseConnection}
        @type database_connection: DatabaseConnection
        """

        super(JobSubmissionAdaptor, self).__init__(
            database_connection=database_connection,
            object_type=JobSubmission,
            table_name='executable',
            column_definition=[
                # Primary key
                ['executable_id', 'INTEGER PRIMARY KEY ASC AUTOINCREMENT'],
                # Name
                ['name', 'TEXT UNIQUE'],
                # Command as submitted into the DRMS
                ['command', 'TEXT'],
            ])

        return

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

        object_list = self.select(statement=statement, parameters=parameters)
        object_length = len(object_list)

        if object_length > 1:
            raise Exception("SQL database returned more than one row for unique field 'name'.")
        elif object_length == 1:
            return object_list[0]
        else:
            return
