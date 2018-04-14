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


from __future__ import print_function

import sqlite3


class DatabaseConnection(object):
    """C{bsf.database.DatabaseConnection} class encapsulating the C{sqlite3.Connection} class.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    """

    def __init__(
            self,
            file_path=None):
        """Initialise a C{bsf.database.DatabaseConnection}.

        The underlying C{sqlite3.Connection} is instantiated only upon calling
        C{bsf.database.DatabaseConnection.connect}.

        @param file_path: File path
        @type file_path: str | unicode
        @return:
        @rtype:
        """
        super(DatabaseConnection, self).__init__()

        if file_path is None:
            self.file_path = ':memory:'
        else:
            self.file_path = file_path

        self._connection = None
        """ @type _connection: sqlite3.Connection | None """

        return

    def __del__(self):
        """Delete a C{bsf.database.DatabaseConnection} object, which implies closing its C{sqlite3.Connection} first.

        @return:
        @rtype:
        """
        if self._connection is None:
            return
        else:
            return self._connection.close()

    def connect(self):
        """Connect a C{bsf.database.DatabaseConnection} by instantiating the underlying C{sqlite3.Connection}.

        Just returns if the C{sqlite3.Connection} already exists.
        @return:
        @rtype:
        """
        if self._connection is None:
            self._connection = sqlite3.connect(database=self.file_path)

        return

    def disconnect(self):
        """Disconnect a C{bsf.database.DatabaseConnection} by closing the underlying C{sqlite3.Connection}.

        Just returns if the C{sqlite3.Connection} object does not exist, replaces it with C{None} otherwise.

        @return:
        @rtype:
        """
        if self._connection is not None:
            self._connection.close()
            self._connection = None

        return

    def commit(self):
        """Commit changes to the underlying C{sqlite3.Connection}.

        @return:
        @rtype:
        """
        if self._connection is not None:
            return self._connection.commit()
        else:
            return

    def get_cursor(self):
        """Get a C{sqlite3.Cursor} of the underlying C{sqlite3.Connection}.

        Creates (connects) the underlying C{sqlite3.Connection} if it does not exist already.

        @return: C{sqlite3.Cursor}
        @rtype: sqlite3.Cursor
        """
        self.connect()

        return self._connection.cursor()


class SQLiteMaster(object):
    """The C{bsf.database.SQLiteMaster} models a row of the SQLite-specific 'sqlite_master' system table.

    Attributes:
    @ivar sql_object_type: SQLite object type
    @type sql_object_type: str
    @ivar sql_object_name: SQLite object name
    @type sql_object_name: str
    @ivar sql_table_name: SQLite table name
    @type sql_table_name: str
    @ivar root_page: B-Tree root page
    @type root_page: int
    @ivar sql_statement: SQLite CREATE TABLE statement
    @type sql_statement: str
    """

    def __init__(self, sql_object_type, sql_object_name, sql_table_name, root_page, sql_statement):
        """Initialise a C{bsf.database.SQLiteMaster} object.

        @param sql_object_type: SQLite object type
        @type sql_object_type: str
        @param sql_object_name: SQLite object name
        @type sql_object_name: str
        @param sql_table_name: SQLite table name
        @type sql_table_name: str
        @param root_page: B-Tree root page
        @type root_page: int
        @param sql_statement: SQLite CREATE TABLE statement
        @type sql_statement: str
        @return:
        @rtype:
        """
        self.sql_object_type = sql_object_type
        self.sql_object_name = sql_object_name
        self.sql_table_name = sql_table_name
        self.root_page = root_page
        self.sql_statement = sql_statement

        return


class SQLiteMasterAdaptor(object):
    """The C{bsf.database.SQLiteMasterAdaptor} models access to the SQLite-specific 'sqlite_master' system table.

    The C{bsf.database.SQLiteMasterAdaptor} does not depend on the C{bsf.database.DatabaseAdaptor} class, since
    the table schema is intrinsic to SQLite so that most methods do not apply.
    Attributes:
    @ivar database_connection: C{bsf.database.DatabaseConnection}
    @type database_connection: bsf.database.DatabaseConnection
    """

    def __init__(self, database_connection):
        """Initialise a C{bsf.database.SQLiteMasterAdaptor} object.

        @param database_connection: C{bsf.database.DatabaseConnection}
        @type database_connection: bsf.database.DatabaseConnection
        @return:
        @rtype:
        """
        assert isinstance(database_connection, DatabaseConnection)

        self.database_connection = database_connection  # Can be None.
        self.table_name = 'sqlite_master'

        return

    @staticmethod
    def _build_column_result_expression():
        """Build a SQL expression of column names typically used in I{SELECT} statements.

        @return: Column result expression string
        @rtype: str
        """
        return ', '.join(('type', 'name', 'tbl_name', 'rootpage', 'sql'))

    def statement_select(self, where_clause=None):
        """Build a SQL I{SELECT} statement.

        @param where_clause: SQL I{WHERE} clause
        @type where_clause: str
        @return: SQL I{SELECT} statement
        @rtype: str
        """
        statement_list = list()
        """ @type statement_list: list[str] """

        statement_list.append('SELECT')
        statement_list.append(self._build_column_result_expression())
        statement_list.append('FROM')
        statement_list.append(self.table_name)

        if where_clause is not None and where_clause:
            statement_list.append('WHERE')
            statement_list.append(where_clause)

        return ' '.join(statement_list)

    def select(self, statement, parameters=None):
        """Execute a SQL I{SELECT} statement and return canonical Python C{object} instances.

        @param statement: Complete SQL I{SELECT} statement
        @type statement: str
        @param parameters: Python C{list} of Python C{str} (parameter) objects or C{None}
        @type parameters: list[str]
        @return: Python C{list} of C{bsf.database.SQLiteMaster} objects
        @rtype: list[bsf.database.SQLiteMaster]
        """
        object_list = list()
        """ @type object_list: list[bsf.database.SQLiteMaster] """

        cursor = self.database_connection.get_cursor()

        if parameters:
            cursor.execute(statement, parameters)
        else:
            cursor.execute(statement)

        for row_tuple in cursor.fetchall():
            object_list.append(SQLiteMaster(*row_tuple))

        return object_list

    def select_all_by_type(self, sql_object_type):
        """Select all C{bsf.database.SQLiteMaster} objects by type.

        @param sql_object_type: SQL object type
        @type sql_object_type: str
        @return: Python C{list} of C{bsf.database.SQLiteMaster} objects
        @rtype: list[bsf.database.SQLiteMaster]
        """
        return self.select(
            statement=self.statement_select(where_clause='type = ?'),
            parameters=[sql_object_type])

    def select_by_type_and_name(self, sql_object_type, sql_object_name):
        """Select a C{bsf.database.SQLiteMaster} object by SQL type and SQL name.

        @param sql_object_type: SQL object type
        @type sql_object_type: str
        @param sql_object_name: SQL object name
        @type sql_object_name: str
        @return: C{bsf.database.SQLiteMaster}
        @rtype: bsf.database.SQLiteMaster | None
        """
        object_list = self.select(
            statement=self.statement_select(where_clause='type = ? and name = ?'),
            parameters=[sql_object_type, sql_object_name])

        if len(object_list) == 1:
            return object_list[0]
        else:
            return


class SQLiteTableInfo(object):
    """SQLite-specific PRAGMA table_info() row object

    Attributes:
    @ivar column_identifier: SQLite column identifier
    @type column_identifier: int
    @ivar column_name: SQLite column name
    @type column_name: str
    @ivar column_type: SQLite column type
    @type column_type: str
    @ivar column_not_null: Column not NULL
    @type column_not_null: int
    @ivar default_value: Default value
    @type default_value: str
    @ivar primary_key: Primary key
    @type primary_key: str
    """

    def __init__(
            self,
            column_identifier=0,
            column_name=None,
            column_type=None,
            column_not_null=0,
            default_value=None,
            primary_key=None):
        """Initialise a C{bsf.database.SQLiteTableInfo} object.

        @param column_identifier: SQLite column identifier
        @type column_identifier: int
        @param column_name: SQLite column name
        @type column_name: str
        @param column_type: SQLite column type
        @type column_type: str
        @param column_not_null: Column not NULL
        @type column_not_null: int
        @param default_value: Default value
        @type default_value: str
        @param primary_key: Primary key
        @type primary_key: str
        """
        self.column_identifier = column_identifier
        self.column_name = column_name  # Can be None.
        self.column_type = column_type  # Can be None.
        self.not_null = column_not_null
        self.default_value = default_value
        self.primary_key = primary_key

        return


class SQLiteTableInfoAdaptor(object):
    """SQLite-specific PRAGMA table_info() adaptor.

    Attributes:
    @ivar database_connection: C{bsf.database.DatabaseConnection}
    @type database_connection: bsf.database.DatabaseConnection
    """

    def __init__(self, database_connection):
        """Initialise a C{bsf.database.SQLiteTableInfoAdaptor}.

        @param database_connection: C{bsf.database.DatabaseConnection}
        @type database_connection: bsf.database.DatabaseConnection
        @return:
        @rtype:
        """
        assert isinstance(database_connection, DatabaseConnection)

        self.database_connection = database_connection  # Can be None.

        return

    @staticmethod
    def statement_pragma_table_info(table_name):
        """Build a SQLite I{PRAGMA table_info} statement.

        @param table_name: SQL table name
        @type table_name: str
        @return: SQLite I{PRAGMA table_info} statement
        @rtype: str
        """
        return "PRAGMA table_info('" + table_name + "')"

    def select_all_by_table_name(self, table_name):
        """Select all C{bsf.database.SQLiteTableInfo} objects by SQLite table name.

        @param table_name: SQLite table name
        @type table_name: str
        @return: Python C{list} of C{bsf.database.SQLiteTableInfo} objects
        @rtype: list[bsf.database.SQLiteTableInfo]
        """
        object_list = list()
        """ @type object_list: list[bsf.database.SQLiteTableInfo] """

        cursor = self.database_connection.get_cursor()

        cursor.execute(self.statement_pragma_table_info(table_name=table_name))

        for row_tuple in cursor.fetchall():
            object_list.append(SQLiteTableInfo(*row_tuple))

        return object_list


class DatabaseAdaptor(object):
    """C{bsf.database.DatabaseAdaptor} class representing as a super-class of object-specific table adaptors.

    Instance variables should be overridden in sub-classes.

    Attributes:
    @ivar database_connection: C{bsf.database.DatabaseConnection}
    @type database_connection: bsf.database.DatabaseConnection
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
        """Initialise a C{bsf.database.DatabaseAdaptor}.

        @param database_connection: C{bsf.database.DatabaseConnection}
        @type database_connection: bsf.database.DatabaseConnection
        @param object_type: Object type
        @type object_type: type
        @param table_name: SQL database table name
        @type table_name: str
        @param column_definition: Python C{list} of Python C{list} objects with
            Python C{str} (SQL column name) and Python C{str} (SQL column constraint) objects
        @type column_definition: list[list[str]]
        @param table_constraint: SQL table constraint expression
        @type table_constraint: list[str]
        @return:
        @rtype:
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
        """Build a SQL expression of column names typically used in I{SELECT} statements.

        This method simply lists all column names of the column definition.

        @return: Column result expression string
        @rtype: str
        """
        return ', '.join(self._get_column_name_list_with_primary())

    def _build_column_definition_expression(self):
        """Build a SQL expression of column definitions typically used in I{CREATE TABLE} statements.

        @return: Column definition expression string
        @rtype: str
        """
        return ', '.join(map(lambda x: ' '.join((x[0], x[1])), self.column_definition))

    def _build_column_insert_expression(self):
        """Build a SQL expression of column names typically used in I{INSERT} statements.

        This method excludes C{PRIMARY KEY} columns with definition I{AUTOINCREMENT},
        which must not be assigned a value.

        @return: Column definition expression string
        @rtype: str
        """
        return ', '.join(self._get_column_name_list_without_primary())

    def _build_value_insert_expression(self):
        """Build a SQL expression of value placeholders (?) typically used in I{INSERT} statements.

        @return: Column value expression string
        @rtype: str
        """
        return ', '.join(map(lambda x: '?', self._get_column_name_list_without_primary()))

    def _build_column_update_expression(self):
        """Build a SQL expression of column name and value placeholder pairs typically used in SQL I{UPDATE} statements.

        As in I{INSERT} expressions, leave out the I{PRIMARY KEY} columns with definition I{AUTOINCREMENT}.

        @return: SQL column name and value placeholder pair expression string
        @rtype: str
        """
        return ', '.join(map(lambda x: ' '.join((x, '=', '?')), self._get_column_name_list_without_primary()))

    def connect(self):
        """Convenience method to connect a C{bsf.database.DatabaseAdaptor}.

        This method instantiates the underlying C{sqlite3.Connection} via its C{bsf.database.DatabaseConnection}.
        Just returns if the C{sqlite3.Connection} object already exists.

        @return:
        @rtype:
        """
        return self.database_connection.connect()

    def disconnect(self):
        """Convenience method to explicitly disconnect a C{bsf.database.DatabaseAdaptor}.

        This method disconnects the underlying C{sqlite3.Connection} via the C{bsf.database.DatabaseConnection}.

        @return:
        @rtype:
        """
        return self.database_connection.disconnect()

    def commit(self):
        """Convenience method to commit changes to a C{bsf.database.DatabaseAdaptor}.

        This method commits to the underlying C{sqlite3.Connection} via the C{bsf.database.DatabaseConnection}.

        @return:
        @rtype:
        """
        return self.database_connection.commit()

    def get_cursor(self):
        """Convenience method to get a C{sqlite3.Cursor} of a C{bsf.database.DatabaseAdaptor}.

        This method gets the C{sqlite3.Cursor} from the C{sqlite3.Connection} of the underlying
        C{bsf.database.DatabaseConnection}.

        @return: C{sqlite3.Cursor}
        @rtype: sqlite3.Cursor
        """
        return self.database_connection.get_cursor()

    def statement_alter_table_rename(self, table_name_old=None, table_name_new=None):
        """Build a SQL I{ALTER TABLE} statement.

        @param table_name_old: Old table name, defaults to table_name.
        @type table_name_old: str
        @param table_name_new: New table name, defaults to table_name_altered.
        @type table_name_new: str
        @return: SQL I{ALTER TABLE table RENAME TO table_name} statement
        @rtype: str
        """
        if not table_name_old:
            table_name_old = self.table_name

        if not table_name_new:
            table_name_new = '_'.join((self.table_name, 'altered'))

        statement_list = list()
        """ @type statement_list: list[str] """

        statement_list.append('ALTER')
        statement_list.append('TABLE')
        statement_list.append("'" + table_name_old + "'")
        statement_list.append('RENAME')
        statement_list.append('TO')
        statement_list.append(table_name_new)

        return ' '.join(statement_list)

    def statement_create_table(self):
        """Build a SQL I{CREATE TABLE} statement.

        @return: SQL I{CREATE TABLE} statement
        @rtype: str
        """
        statement_list = list()
        """ @type statement_list: list[str] """

        statement_list.append('CREATE')
        statement_list.append('TABLE')
        statement_list.append("'" + self.table_name + "'")
        statement_list.append('(' + self._build_column_definition_expression() + ')')

        return ' '.join(statement_list)

    @staticmethod
    def statement_drop_table(table_name):
        """Build a SQL I{DROP TABLE} statement.

        @return: SQL I{DROP TABLE} statement
        @rtype: str
        """
        statement_list = list()
        """ @type statement_list: list[str] """

        statement_list.append('DROP')
        statement_list.append('TABLE')
        statement_list.append('IF')
        statement_list.append('EXISTS')
        statement_list.append("'" + table_name + "'")

        return ' '.join(statement_list)

    def statement_insert(self):
        """Build a SQL I{INSERT INTO} statement.

        @return: SQL I{INSERT INTO} statement
        @rtype: str
        """
        statement_list = list()
        """ @type statement_list: list[str] """

        statement_list.append('INSERT')
        statement_list.append('INTO')
        statement_list.append("'" + self.table_name + "'")
        statement_list.append('(' + self._build_column_insert_expression() + ')')
        statement_list.append('VALUES')
        statement_list.append('(' + self._build_value_insert_expression() + ')')

        return ' '.join(statement_list)

    def statement_select(self, where_clause=None, group_clause=None, having_clause=None):
        """Build a SQL I{SELECT} statement.

        @param where_clause: SQL I{WHERE} clause
        @type where_clause: str
        @param group_clause: SQL I{GROUP BY} clause
        @type group_clause: str
        @param having_clause: SQL I{HAVING} clause
        @type having_clause: str
        @return: SQL I{SELECT} statement
        @rtype: str
        """
        statement_list = list()
        """ @type statement_list: list[str] """

        statement_list.append('SELECT')
        statement_list.append(self._build_column_result_expression())
        statement_list.append('FROM')
        statement_list.append("'" + self.table_name + "'")

        if where_clause is not None and len(where_clause):
            statement_list.append('WHERE')
            statement_list.append(where_clause)

        if group_clause is not None and len(group_clause):
            statement_list.append('GROUP BY')
            statement_list.append(group_clause)

        if having_clause is not None and len(having_clause):
            statement_list.append('HAVING')
            statement_list.append(having_clause)

        return ' '.join(statement_list)

    def create_table(self):
        """Execute a SQL I{CREATE TABLE} statement for the canonical C{bsf.database.DatabaseAdaptor} table.

        Before attempting to execute the SQL I{CREATE TABLE} statement, this method checks in 'sqlite_master',
        whether the table already exists in the SQLite database.
        After calling, the C{bsf.database.DatabaseAdaptor.commit} method has to be called at some stage.

        @return:
        @rtype:
        """

        sqlite_master_adaptor = SQLiteMasterAdaptor(database_connection=self.database_connection)

        sqlite_master = sqlite_master_adaptor.select_by_type_and_name(
            sql_object_type='table',
            sql_object_name=self.table_name)

        if sqlite_master is None:
            self.get_cursor().execute(self.statement_create_table())

        return

    def select(self, statement, parameters=None):
        """Execute a SQL I{SELECT} statement and return canonical Python C{object} instances.

        @param statement: Complete SQL I{SELECT} statement
        @type statement: str
        @param parameters: Python C{list} of Python C{str} (parameter) objects or C{None}
        @type parameters: list[None | int | float | str | unicode]
        @return: Python C{list} of Python C{object} objects
        @rtype: list[object]
        """
        object_list = list()

        cursor = self.get_cursor()

        if parameters:
            cursor.execute(statement, parameters)
        else:
            cursor.execute(statement)

        for row_tuple in cursor.fetchall():
            object_instance = self.object_type()
            object_list.append(object_instance)
            i = 0
            for name in map(lambda x: x[0], self.column_definition):
                object_instance.__setattr__(name, row_tuple[i])
                i += 1

        return object_list

    def select_all(self):
        """Select all canonical Python C{object} instances.

        @return: Python C{list} of Python C{object} objects
        @rtype: list[object]
        """
        return self.select(statement=self.statement_select())

    def select_by_identifier(self, identifier=0):
        """Select one canonical Python C{object} instance corresponding to the primary key identifier.

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

        statement = self.statement_select(where_clause=' '.join((primary_key, '=', '?')))
        parameters.append(identifier)

        object_list = self.select(statement=statement, parameters=parameters)
        object_length = len(object_list)

        if object_length > 1:
            raise Exception('SQL database returned more than one row for primary key ' + repr(primary_key) + '.')
        elif object_length == 1:
            return object_list[0]
        else:
            return

    def insert(self, object_instance):
        """Execute a SQL I{INSERT} statement for a canonical Python C{object} instance.

        @param object_instance: Python C{object} object
        @type object_instance: object
        @return:
        @rtype:
        """
        assert isinstance(object_instance, self.object_type)

        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        value_list = map(lambda x: object_instance.__getattribute__(x), self._get_column_name_list_without_primary())

        try:
            self.get_cursor().execute(self.statement_insert(), value_list)
        except sqlite3.IntegrityError:
            print('Encountered sqlite3.IntegrityError for table name ' +
                  self.table_name + 'on the following SQL fields:\n' +
                  'Fields: ' + self._build_column_insert_expression() + '\n' +
                  'Values: ' + repr(value_list))
            raise
        except sqlite3.OperationalError:
            print('Encountered sqlite3.OperationalError for table name ' +
                  self.table_name + ' on the following SQL fields:\n' +
                  'Fields: ' + self._build_column_insert_expression() + '\n' +
                  'Values: ' + repr(value_list))
            raise

        # Update the canonical attribute containing the primary key with the last row identifier.
        last_row_identifier = self.get_cursor().lastrowid
        column_name = self.table_name + '_id'
        if last_row_identifier and hasattr(object_instance, column_name):
            object_instance.__setattr__(column_name, last_row_identifier)

        return

    def update(self, object_instance):
        """Execute a SQL I{UPDATE} statement for a canonical Python C{object} instance.

        @param object_instance: Python C{object} instance
        @type object_instance: object
        @return:
        @rtype:
        """
        assert isinstance(object_instance, self.object_type)
        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        value_list = map(lambda x: object_instance.__getattribute__(x), self._get_column_name_list_without_primary())

        primary_name = self._get_column_name_for_primary()
        if primary_name:
            value_list.append(object_instance.__getattribute__(primary_name))
        else:
            raise Exception('Cannot update table ' + repr(self.table_name) + ' without primary key.')

        statement_list = list()
        """ @type statement_list: list[str] """

        statement_list.append('UPDATE')
        statement_list.append("'" + self.table_name + "'")
        statement_list.append('SET')
        statement_list.append(self._build_column_update_expression())
        statement_list.append('WHERE')
        statement_list.append("'" + primary_name + "'")
        statement_list.append('=')
        statement_list.append('?')

        try:
            self.get_cursor().execute(' '.join(statement_list), value_list)
        except sqlite3.IntegrityError:
            print('Encountered SQLite3 integrity error.\n',
                  '  SQL statement:', ' '.join(statement_list), '\n',
                  '  Values:', value_list)
            raise

        return

    def compare_table_definitions(self):
        """Compare the current table definition to the SQLite PRAGMA table_info().

        @return: Python C{tuple} of Python C{dict} objects of
            Python C{str} (column name) key and
            Python C{None} value for the table_info() and the column definition.
        @rtype: (dict[str, None], dict[str, None])
        """
        pragma_table_info_adaptor = SQLiteTableInfoAdaptor(
            database_connection=self.database_connection)

        pragma_table_info_list = pragma_table_info_adaptor.select_all_by_table_name(
            table_name=self.table_name)

        column_dict_old = dict(map(lambda x: (x, None), map(lambda x: x.column_name, pragma_table_info_list)))
        """ @type column_dict_old: dict[str, None] """
        column_dict_new = dict(map(lambda x: (x, None), map(lambda x: x[0], self.column_definition)))
        """ @type column_dict_new: dict[str, None] """

        for key in column_dict_new.keys():
            if key in column_dict_old:
                del column_dict_new[key]
                del column_dict_old[key]

        return column_dict_old, column_dict_new


class JobSubmission(object):
    """C{bsf.database.JobSubmission} class representing a Distributed Resource Management System (DRMS) process.

    This class is equivalent to the C{bsf.process.Executable} and C{bsf.process.Command} classes, but much less complex.
    Command lines are stored as submitted and not broken down into sub-commands, options and arguments.

    Attributes:
    @ivar executable_id: Primary key
    @type executable_id: int
    @ivar name: C{bsf.process.Executable.name}
    @type name: str
    @ivar command: Command line
    @type command: str
    """

    def __init__(
            self,
            executable_id=0,
            name=None,
            command=None):
        """Initialise a C{bsf.database.JobSubmission}.

        @param executable_id: Primary key
        @type executable_id: int
        @param name: C{bsf.process.Executable.name}
        @type name: str
        @param command: Command line
        @type command: str
        @return:
        @rtype:
        """
        super(JobSubmission, self).__init__()

        self.executable_id = executable_id  # Can be 0.
        self.name = name  # Can be None.
        self.command = command  # Can be None.

        return


class JobSubmissionAdaptor(DatabaseAdaptor):
    """C{bsf.database.JobSubmissionAdaptor} class providing database access for C{bsf.database.JobSubmission}.

    Attributes:
    """

    def __init__(
            self,
            database_connection):
        """Initialise a C{bsf.database.JobSubmissionAdaptor}.

        @param database_connection: C{bsf.database.DatabaseConnection}
        @type database_connection: bsf.database.DatabaseConnection
        @return:
        @rtype:
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
                # Command as submitted into the Stage
                ['command', 'TEXT'],
            ])

        return

    def select_by_name(self, name):
        """Select one C{bsf.database.JobSubmission} object by name.

        @param name: Name
        @type name: str
        @return: C{bsf.database.JobSubmission} or C{None}
        @rtype: bsf.database.JobSubmission | None
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
