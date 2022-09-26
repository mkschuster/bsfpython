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
"""The :py:mod:`bsf.database` module provides classes and methods to centralises :literal:`SQLite` Database access.
"""
from sqlite3 import Connection, Cursor, IntegrityError, OperationalError, connect
from typing import Dict, List


class DatabaseConnection(object):
    """The :py:class:`bsf.database.DatabaseConnection` class encapsulates the :py:class:`sqlite3.Connection` class.

    :ivar file_path: A file path.
    :type file_path: str
    :ivar _connection: A :py:class:`sqlite3.Connection` object.
    :type _connection: Connection | None
    """

    def __init__(
            self,
            file_path=None):
        """Initialise a :py:class:`bsf.database.DatabaseConnection` object.

        The underlying :py:class:`sqlite3.Connection` object is instantiated only upon calling
        :py:meth:`bsf.database.DatabaseConnection.connect`.

        :param file_path: A file path.
        :type file_path: str | None
        """
        super(DatabaseConnection, self).__init__()

        if file_path is None:
            self.file_path = ':memory:'
        else:
            self.file_path = file_path

        self._connection = None

        return

    def __del__(self):
        """Delete a :py:class:`bsf.database.DatabaseConnection` object,
        which implies closing its :py:class:`sqlite3.Connection` object first.
        """
        if self._connection is None:
            return
        else:
            return self._connection.close()

    def connect(self):
        """Connect a :py:class:`bsf.database.DatabaseConnection` object by
        instantiating the underlying :py:class:`sqlite3.Connection` object.

        Just returns if the :py:class:`sqlite3.Connection` object already exists.
        """
        if self._connection is None:
            self._connection = connect(database=self.file_path)

        return

    def disconnect(self):
        """Disconnect a :py:class:`bsf.database.DatabaseConnection` object by
        closing the underlying :py:class:`sqlite3.Connection` object.

        Just returns if the :py:class:`sqlite3.Connection` object does not exist,
        replaces it with :py:const:`None` otherwise.
        """
        if self._connection is not None:
            self._connection.close()
            self._connection = None

        return

    def commit(self):
        """Commit changes to the underlying :py:class:`sqlite3.Connection` object.
        """
        if self._connection is not None:
            return self._connection.commit()
        else:
            return

    def get_cursor(self):
        """Get a :py:class:`sqlite3.Cursor` object of the underlying :py:class:`sqlite3.Connection` object.

        Creates (connects) the underlying :py:class:`sqlite3.Connection` object if it does not exist already.

        :return: A :py:class:`sqlite3.Cursor` object.
        :rtype: Cursor
        """
        self.connect()

        return self._connection.cursor()


class SQLiteMaster(object):
    """The :py:class:`bsf.database.SQLiteMaster` object models a row of the
    :literal:`SQLite`-specific :literal:`sqlite_master` system table.

    :ivar sql_object_type: A :literal:`SQLite` object type.
    :type sql_object_type: str | None
    :ivar sql_object_name: A :literal:`SQLite` object name.
    :type sql_object_name: str | None
    :ivar sql_table_name: A :literal:`SQLite` table name.
    :type sql_table_name: str | None
    :ivar root_page: The B-Tree root page.
    :type root_page: int | None
    :ivar sql_statement: A :literal:`SQLite` :literal:`CREATE TABLE` statement.
    :type sql_statement: str | None
    """

    def __init__(
            self,
            sql_object_type=None,
            sql_object_name=None,
            sql_table_name=None,
            root_page=None,
            sql_statement=None):
        """Initialise a :py:class:`bsf.database.SQLiteMaster` object.

        :param sql_object_type: A :literal:`SQLite` object type.
        :type sql_object_type: str | None
        :param sql_object_name: A :literal:`SQLite` object name.
        :type sql_object_name: str | None
        :param sql_table_name: A :literal:`SQLite` table name.
        :type sql_table_name: str | None
        :param root_page: The B-Tree root page.
        :type root_page: int | None
        :param sql_statement: A :literal:`SQLite` :literal:`CREATE TABLE` statement.
        :type sql_statement: str | None
        """
        super(SQLiteMaster, self).__init__()

        self.sql_object_type = sql_object_type
        self.sql_object_name = sql_object_name
        self.sql_table_name = sql_table_name
        self.root_page = root_page
        self.sql_statement = sql_statement

        return


class SQLiteMasterAdaptor(object):
    """The :py:class:`bsf.database.SQLiteMasterAdaptor` class models access to the
    :literal:`SQLite`-specific :literal:`sqlite_master` system table.

    The :py:class:`bsf.database.SQLiteMasterAdaptor` class does not depend on the
    :py:class:`bsf.database.DatabaseAdaptor` class, since
    the table schema is intrinsic to :literal:`SQLite` so that most methods do not apply.

    :ivar database_connection: A :py:class:`bsf.database.DatabaseConnection` object.
    :type database_connection: DatabaseConnection
    """

    def __init__(self, database_connection):
        """Initialise a :py:class:`bsf.database.SQLiteMasterAdaptor` object.

        :param database_connection: A :py:class:`bsf.database.DatabaseConnection` object.
        :type database_connection: DatabaseConnection
        """
        super(SQLiteMasterAdaptor, self).__init__()

        self.database_connection = database_connection
        self.table_name = 'sqlite_master'

        return

    @staticmethod
    def _build_column_result_expression():
        """Build a SQL expression of column names typically used in :literal:`SELECT` statements.

        :return: A column result expression string.
        :rtype: str
        """
        return ', '.join(('type', 'name', 'tbl_name', 'rootpage', 'sql'))

    def statement_select(self, where_clause=None):
        """Build a SQL :literal:`SELECT` statement.

        :param where_clause: A SQL :literal:`WHERE` clause.
        :type where_clause: str | None
        :return: A SQL :literal:`SELECT` statement.
        :rtype: str
        """
        statement_list: List[str] = list()

        statement_list.append('SELECT')
        statement_list.append(self._build_column_result_expression())
        statement_list.append('FROM')
        statement_list.append(self.table_name)

        if where_clause:
            statement_list.append('WHERE')
            statement_list.append(where_clause)

        return ' '.join(statement_list)

    def select(self, statement, parameters=None):
        """Execute a SQL :literal:`SELECT` statement and return canonical Python :py:class:`object` instances.

        :param statement: A complete SQL :literal:`SELECT` statement.
        :type statement: str
        :param parameters: A Python :py:class:`list` object of
            Python :py:class:`str` (parameter) objects or :py:const:`None`.
        :type parameters: list[int | float | str | None] | None
        :return: A Python :py:class:`list` object of :py:class:`bsf.database.SQLiteMaster` objects.
        :rtype: list[SQLiteMaster]
        """
        object_list: List[SQLiteMaster] = list()

        cursor = self.database_connection.get_cursor()

        if parameters:
            cursor.execute(statement, parameters)
        else:
            cursor.execute(statement)

        for row_tuple in cursor.fetchall():
            object_list.append(SQLiteMaster(*row_tuple))

        return object_list

    def select_all_by_type(self, sql_object_type):
        """Select all :py:class:`bsf.database.SQLiteMaster` objects by type.

        :param sql_object_type: A SQL object type.
        :type sql_object_type: str
        :return: A Python :py:class:`list` object of :py:class:`bsf.database.SQLiteMaster` objects.
        :rtype: list[SQLiteMaster]
        """
        return self.select(
            statement=self.statement_select(where_clause='type = ?'),
            parameters=[sql_object_type])

    def select_by_type_and_name(self, sql_object_type, sql_object_name):
        """Select a :py:class:`bsf.database.SQLiteMaster` object by SQL type and SQL name.

        :param sql_object_type: A SQL object type.
        :type sql_object_type: str
        :param sql_object_name: A SQL object name.
        :type sql_object_name: str
        :return: A :py:class:`bsf.database.SQLiteMaster` object.
        :rtype: SQLiteMaster | None
        """
        object_list = self.select(
            statement=self.statement_select(where_clause='type = ? and name = ?'),
            parameters=[sql_object_type, sql_object_name])

        if len(object_list) == 1:
            return object_list[0]
        else:
            return


class SQLiteTableInfo(object):
    """The :py:class:`bsf.database.SQLiteTableInfo` class models a
    :literal:`SQLite`-specific :literal:`PRAGMA table_info()` row object.

    :ivar column_identifier: A :literal:`SQLite` column identifier.
    :type column_identifier: int | None
    :ivar column_name: A :literal:`SQLite` column name.
    :type column_name: str | None
    :ivar column_type: A :literal:`SQLite` column type.
    :type column_type: str | None
    :ivar column_not_null: A column not NULL.
    :type column_not_null: int | None
    :ivar default_value: A default value.
    :type default_value: str | None
    :ivar primary_key: A primary key.
    :type primary_key: str | None
    """

    def __init__(
            self,
            column_identifier=None,
            column_name=None,
            column_type=None,
            column_not_null=None,
            default_value=None,
            primary_key=None):
        """Initialise a :py:class:`bsf.database.SQLiteTableInfo` object.

        :param column_identifier: A :literal:`SQLite` column identifier.
        :type column_identifier: int | None
        :param column_name: A :literal:`SQLite` column name.
        :type column_name: str | None
        :param column_type: A :literal:`SQLite` column type.
        :type column_type: str | None
        :param column_not_null: A column not NULL.
        :type column_not_null: int | None
        :param default_value: A default value.
        :type default_value: str | None
        :param primary_key: A primary key.
        :type primary_key: str | None
        """
        super(SQLiteTableInfo, self).__init__()

        self.column_identifier = column_identifier
        self.column_name = column_name
        self.column_type = column_type
        self.column_not_null = column_not_null
        self.default_value = default_value
        self.primary_key = primary_key

        return


class SQLiteTableInfoAdaptor(object):
    """The :py:class:`bsf.database.SQLiteTableInfoAdaptor` class models a
    :literal:`SQLite`-specific :literal:`PRAGMA table_info()` adaptor.

    :ivar database_connection: A :py:class:`bsf.database.DatabaseConnection` object.
    :type database_connection: DatabaseConnection
    """

    def __init__(self, database_connection):
        """Initialise a :py:class:`bsf.database.SQLiteTableInfoAdaptor` object.

        :param database_connection: A :py:class:`bsf.database.DatabaseConnection` object.
        :type database_connection: DatabaseConnection
        """
        assert isinstance(database_connection, DatabaseConnection)

        super(SQLiteTableInfoAdaptor, self).__init__()

        self.database_connection = database_connection

        return

    @staticmethod
    def statement_pragma_table_info(table_name):
        """Build a :literal:`SQLite` :literal:`PRAGMA table_info` statement.

        :param table_name: A SQL table name.
        :type table_name: str
        :return: A :literal:`SQLite` :literal:`PRAGMA table_info` statement.
        :rtype: str
        """
        return "PRAGMA table_info('" + table_name + "')"

    def select_all_by_table_name(self, table_name):
        """Select all :py:class:`bsf.database.SQLiteTableInfo` objects by SQLite table name.

        :param table_name: A SQLite table name.
        :type table_name: str
        :return: A Python :py:class:`list` object of :py:class:`bsf.database.SQLiteTableInfo` objects.
        :rtype: list[SQLiteTableInfo]
        """
        object_list: List[SQLiteTableInfo] = list()

        cursor = self.database_connection.get_cursor()

        cursor.execute(self.statement_pragma_table_info(table_name=table_name))

        for row_tuple in cursor.fetchall():
            object_list.append(SQLiteTableInfo(*row_tuple))

        return object_list


class DatabaseAdaptor(object):
    """The :py:class:`bsf.database.DatabaseAdaptor` class represents a super-class of object-specific table adaptors.

    Instance variables should be overridden in subclasses.

    :ivar database_connection: A :py:class:`bsf.database.DatabaseConnection` object.
    :type database_connection: DatabaseConnection
    :ivar table_name: A SQL database table name.
    :type table_name: str
    :ivar column_definition: A Python :py:class:`list` object of
        Python :py:class:`tuple` objects of
        Python :py:class:`str` (SQL column name) and
        Python :py:class:`str` (SQL column constraint) objects.
    :type column_definition: list[(str, str)]
    :ivar table_constraint: A SQL table constraint expression.
    :type table_constraint: list[str]
    """

    def __init__(
            self,
            database_connection,
            object_type,
            table_name,
            column_definition,
            table_constraint=None):
        """Initialise a :py:class:`bsf.database.DatabaseAdaptor` object.

        :param database_connection: A :py:class:`bsf.database.DatabaseConnection` object.
        :type database_connection: DatabaseConnection
        :param object_type: An object type.
        :type object_type: type
        :param table_name: A SQL database table name.
        :type table_name: str
        :param column_definition: A Python :py:class:`list` object of
            Python :py:class:`tuple` objects of
            Python :py:class:`str` (SQL column name) and
            Python :py:class:`str` (SQL column constraint) objects.
        :type column_definition: list[(str, str)]
        :param table_constraint: A SQL table constraint expression.
        :type table_constraint: list[str] | None
        """
        assert isinstance(database_connection, DatabaseConnection)
        assert isinstance(object_type, type)

        super(DatabaseAdaptor, self).__init__()

        self.database_connection = database_connection
        self.object_type = object_type
        self.table_name = table_name
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
        """Build a Python :py:class:`list` object of SQL column names including the primary key.

        :return: A Python :py:class:`list` object of Python :py:class:`str` (SQL column name) objects.
        :rtype: list[str]
        """
        return map(lambda x: x[0], self.column_definition)

    def _get_column_name_list_without_primary(self):
        """Build a Python :py:class:`list` object of SQL column names excluding the primary key.

        This method excludes :literal:`PRIMARY KEY` columns with definition :py:class:`AUTOINCREMENT`,
        which must not be assigned a value in :py:class:`INSERT` or :py:class:`UPDATE` statements.

        :return: A Python :py:class:`list` object of Python :py:class:`str` (SQL column name) objects.
        :rtype: list[str]
        """
        return map(lambda x: x[0], filter(lambda x: 'AUTOINCREMENT' not in x[1], self.column_definition))

    def _get_column_name_for_primary(self):
        """Get the SQL column name for the primary key.

        This method returns the :literal:`PRIMARY KEY` column with definition :literal:`AUTOINCREMENT`.

        :return: A column name for primary key.
        :rtype: str | None
        """
        primary_key = None

        for item in self.column_definition:
            if 'AUTOINCREMENT' in item[1]:
                if primary_key:
                    raise Exception("SQL column definition with more than one 'AUTOINCREMENT' attribute.")
                else:
                    primary_key = item[0]

        return primary_key

    def _build_column_result_expression(self):
        """Build a SQL expression of column names typically used in :literal:`SELECT` statements.

        This method simply lists all column names of the column definition.

        :return: A column result expression string.
        :rtype: str
        """
        return ', '.join(self._get_column_name_list_with_primary())

    def _build_column_definition_expression(self):
        """Build a SQL expression of column definitions typically used in :literal:`CREATE TABLE` statements.

        :return: A column definition expression string.
        :rtype: str
        """
        return ', '.join(map(lambda x: ' '.join((x[0], x[1])), self.column_definition))

    def _build_column_insert_expression(self):
        """Build a SQL expression of column names typically used in :literal:`INSERT` statements.

        This method excludes :literal:`PRIMARY KEY` columns with definition :literal:`AUTOINCREMENT`,
        which must not be assigned a value.

        :return: A column definition expression string.
        :rtype: str
        """
        return ', '.join(self._get_column_name_list_without_primary())

    def _build_value_insert_expression(self):
        """Build a SQL expression of value placeholders (:literal:`?`) typically used in :literal:`INSERT` statements.

        :return: A column value expression string.
        :rtype: str
        """
        return ', '.join(map(lambda x: '?', self._get_column_name_list_without_primary()))

    def _build_column_update_expression(self):
        """Build a SQL expression of column name and value placeholder pairs typically used in
        SQL :literal:`UPDATE` statements.

        As in :literal:`INSERT` expressions, leave out the :literal:`PRIMARY KEY` columns with definition
        :literal:`AUTOINCREMENT`.

        :return: A SQL column name and value placeholder pair expression string.
        :rtype: str
        """
        return ', '.join(map(lambda x: ' '.join((x, '=', '?')), self._get_column_name_list_without_primary()))

    def connect(self):
        """Convenience method to connect a :py:class:`bsf.database.DatabaseAdaptor` object.

        This method instantiates the underlying :py:class:`sqlite3.Connection` object via its
        :py:class:`bsf.database.DatabaseConnection` object.
        Just returns if the :py:class:`sqlite3.Connection` object already exists.
        """
        return self.database_connection.connect()

    def disconnect(self):
        """Convenience method to explicitly disconnect a :py:class:`bsf.database.DatabaseAdaptor` object.

        This method disconnects the underlying :py:class:`sqlite3.Connection` object via the
        :py:class:`bsf.database.DatabaseConnection` object.
        """
        return self.database_connection.disconnect()

    def commit(self):
        """Convenience method to commit changes to a :py:class:`bsf.database.DatabaseAdaptor` object.

        This method commits to the underlying :py:class:`sqlite3.Connection` object via the
        :py:class:`bsf.database.DatabaseConnection` object.
        """
        return self.database_connection.commit()

    def get_cursor(self):
        """Convenience method to get a :py:class:`sqlite3.Cursor` object of a
        :py:class:`bsf.database.DatabaseAdaptor` object.

        This method gets the :py:class:`sqlite3.Cursor` object from the
        :py:class:`sqlite3.Connection` object of the underlying
        :py:class:`bsf.database.DatabaseConnection` object.

        :return: A :py:class:`sqlite3.Cursor` object.
        :rtype: Cursor
        """
        return self.database_connection.get_cursor()

    def statement_alter_table_rename(self, table_name_old=None, table_name_new=None):
        """Build a SQL :literal:`ALTER TABLE` statement.

        :param table_name_old: An old table name, defaults to table_name.
        :type table_name_old: str
        :param table_name_new: A new table name, defaults to table_name_altered.
        :type table_name_new: str
        :return: A SQL :literal:`ALTER TABLE <table> RENAME TO <table_name>` statement.
        :rtype: str
        """
        if not table_name_old:
            table_name_old = self.table_name

        if not table_name_new:
            table_name_new = '_'.join((self.table_name, 'altered'))

        statement_list: List[str] = list()

        statement_list.append('ALTER')
        statement_list.append('TABLE')
        statement_list.append(table_name_old)
        statement_list.append('RENAME')
        statement_list.append('TO')
        statement_list.append(table_name_new)

        return ' '.join(statement_list)

    def statement_create_table(self):
        """Build a SQL :literal:`CREATE TABLE` statement.

        :return: A SQL :literal:`CREATE TABLE` statement.
        :rtype: str
        """
        statement_list: List[str] = list()

        statement_list.append('CREATE')
        statement_list.append('TABLE')
        statement_list.append(self.table_name)
        statement_list.append('(' + self._build_column_definition_expression() + ')')

        return ' '.join(statement_list)

    @staticmethod
    def statement_drop_table(table_name):
        """Build a SQL :literal:`DROP TABLE` statement.

        :return: A SQL :literal:`DROP TABLE` statement.
        :rtype: str
        """
        statement_list: List[str] = list()

        statement_list.append('DROP')
        statement_list.append('TABLE')
        statement_list.append('IF')
        statement_list.append('EXISTS')
        statement_list.append(table_name)

        return ' '.join(statement_list)

    def statement_insert(self):
        """Build a SQL :literal:`INSERT INTO` statement.

        :return: A SQL :literal:`INSERT INTO` statement.
        :rtype: str
        """
        statement_list: List[str] = list()

        statement_list.append('INSERT')
        statement_list.append('INTO')
        statement_list.append(self.table_name)
        statement_list.append('(' + self._build_column_insert_expression() + ')')
        statement_list.append('VALUES')
        statement_list.append('(' + self._build_value_insert_expression() + ')')

        return ' '.join(statement_list)

    def statement_select(self, where_clause=None, group_clause=None, having_clause=None):
        """Build a SQL :literal:`SELECT` statement.

        :param where_clause: A SQL :literal:`WHERE` clause.
        :type where_clause: str | None
        :param group_clause: A SQL :literal:`GROUP BY` clause.
        :type group_clause: str | None
        :param having_clause: A SQL :literal:`HAVING` clause.
        :type having_clause: str | None
        :return: A SQL :literal:`SELECT` statement.
        :rtype: str
        """
        statement_list: List[str] = list()

        statement_list.append('SELECT')
        statement_list.append(self._build_column_result_expression())
        statement_list.append('FROM')
        statement_list.append(self.table_name)

        if where_clause:
            statement_list.append('WHERE')
            statement_list.append(where_clause)

        if group_clause:
            statement_list.append('GROUP BY')
            statement_list.append(group_clause)

        if having_clause:
            statement_list.append('HAVING')
            statement_list.append(having_clause)

        return ' '.join(statement_list)

    def statement_update(self):
        """Build a SQL :literal:`UPDATE` statement.

        :return: A SQL :literal:`UPDATE` statement.
        :rtype: str
        """
        primary_name = self._get_column_name_for_primary()

        if not primary_name:
            raise Exception(
                'Cannot create a SQL UPDATE statement ' + repr(self.table_name) + ' without a primary key.')

        statement_list: List[str] = list()

        statement_list.append('UPDATE')
        statement_list.append(self.table_name)
        statement_list.append('SET')
        statement_list.append(self._build_column_update_expression())
        statement_list.append('WHERE')
        statement_list.append(primary_name)
        statement_list.append('=')
        statement_list.append('?')

        return ' '.join(statement_list)

    def create_table(self):
        """Execute a SQL :literal:`CREATE TABLE` statement for the canonical
        :py:class:`bsf.database.DatabaseAdaptor` object table.

        Before attempting to execute the SQL :literal:`CREATE TABLE` statement,
        this method checks in :literal:`sqlite_master`,
        whether the table already exists in the :literal:`SQLite` database.
        After calling, the :py:meth:`bsf.database.DatabaseAdaptor.commit` method has to be called at some stage.
        """
        sqlite_master_adaptor = SQLiteMasterAdaptor(database_connection=self.database_connection)

        sqlite_master = sqlite_master_adaptor.select_by_type_and_name(
            sql_object_type='table',
            sql_object_name=self.table_name)

        if sqlite_master is None:
            self.get_cursor().execute(self.statement_create_table())

        return

    def select(self, statement, parameters=None):
        """Execute a SQL :literal:`SELECT` statement and return canonical Python :py:class:`object` instances.

        :param statement: A complete SQL :literal:`SELECT` statement.
        :type statement: str
        :param parameters: A Python :py:class:`list` object of
            Python :py:class:`str` (parameter) objects or :py:const:`None`.
        :type parameters: list[int | float | str | None] | None
        :return: A Python :py:class:`list` object of Python :py:class:`object` objects.
        :rtype: list[object]
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
        """Select all canonical Python :py:class:`object` instances.

        :return: A Python :py:class:`list` object of Python :py:class:`object` objects.
        :rtype: list[object]
        """
        return self.select(statement=self.statement_select())

    def select_by_identifier(self, identifier=0):
        """Select one canonical Python :py:class:`object` instance corresponding to the primary key identifier.

        :param identifier: A primary key identifier.
        :type identifier: int
        :return: A Python :py:class:`object` instance.
        :rtype: object
        """
        # Check if the table has a primary identifier.
        primary_key = self._get_column_name_for_primary()
        if not primary_key:
            return

        object_list = self.select(
            statement=self.statement_select(where_clause=' '.join((primary_key, '=', '?'))),
            parameters=[identifier])
        object_length = len(object_list)

        if object_length > 1:
            raise Exception('SQL database returned more than one row for primary key ' + repr(primary_key) + '.')
        elif object_length == 1:
            return object_list[0]
        else:
            return

    def insert(self, object_instance):
        """Execute a SQL :literal:`INSERT` statement for a canonical Python :py:class:`object` instance.

        :param object_instance: A Python :py:class:`object` object.
        :type object_instance: object
        """
        assert isinstance(object_instance, self.object_type)

        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        value_list = [object_instance.__getattribute__(x) for x in self._get_column_name_list_without_primary()]

        try:
            self.get_cursor().execute(self.statement_insert(), value_list)
        except IntegrityError:
            print('Encountered sqlite3.IntegrityError for table name ' +
                  self.table_name + 'on the following SQL fields:\n' +
                  'Fields: ' + self._build_column_insert_expression() + '\n' +
                  'Values: ' + repr(value_list))
            raise
        except OperationalError:
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
        """Execute a SQL :literal:`UPDATE` statement for a canonical Python :py:class:`object` instance.

        :param object_instance: A Python :py:class:`object` instance.
        :type object_instance: object
        """
        assert isinstance(object_instance, self.object_type)
        # Get the list of values by using the column definition and reading attributes of the same name
        # from the Python object.

        primary_name = self._get_column_name_for_primary()

        if not primary_name:
            raise Exception('Cannot update table ' + repr(self.table_name) + ' without a primary key.')

        value_list = [object_instance.__getattribute__(item) for item in self._get_column_name_list_without_primary()]
        value_list.append(object_instance.__getattribute__(primary_name))

        try:
            self.get_cursor().execute(self.statement_update(), value_list)
        except IntegrityError:
            print('Encountered SQLite3 integrity error.\n',
                  '  SQL statement:', self.statement_update(), '\n',
                  '  Values:', value_list)
            raise

        return

    def compare_table_definitions(self):
        """Compare the current table definition to the :literal:`SQLite` :literal:`PRAGMA table_info()`.

        :return: A Python :py:class:`tuple` object of
            Python :py:class:`dict` objects of
            Python :py:class:`str` (column name) key and
            Python :py:const:`None` value for the :literal:`PRAGMA table_info()` and the column definition.
        :rtype: (dict[str, None], dict[str, None])
        """
        pragma_table_info_adaptor = SQLiteTableInfoAdaptor(
            database_connection=self.database_connection)

        pragma_table_info_list = pragma_table_info_adaptor.select_all_by_table_name(
            table_name=self.table_name)

        column_dict_old: Dict[str, None] = dict(
            map(lambda x: (x, None), map(lambda x: x.column_name, pragma_table_info_list)))

        column_dict_new: Dict[str, None] = dict(
            map(lambda x: (x, None), map(lambda x: x[0], self.column_definition)))

        # Use a list comprehension to create a list of key objects since the dict gets modified in the loop.
        for key in [key for key in column_dict_new]:
            if key in column_dict_old:
                del column_dict_new[key]
                del column_dict_old[key]

        return column_dict_old, column_dict_new


class JobSubmission(object):
    """The :py:class:`bsf.database.JobSubmission` class represents a
    :literal:`Distributed Resource Management System` (DRMS) process.

    This class is equivalent to the :py:class:`bsf.process.Executable` and :py:class:`bsf.process.Command` classes,
    but much less complex. Command lines are stored as submitted and not broken down into sub-commands,
    options and arguments.

    :ivar executable_id: A primary key.
    :type executable_id: int | None
    :ivar name: A :py:attr:`bsf.process.Executable.name` attribute.
    :type name: str | None
    :ivar command: A command line.
    :type command: str | None
    """

    def __init__(
            self,
            executable_id=None,
            name=None,
            command=None):
        """Initialise a :py:class:`bsf.database.JobSubmission` object.

        :param executable_id: A primary key.
        :type executable_id: int | None
        :param name: A :py:attr:`bsf.process.Executable.name` attribute.
        :type name: str | None
        :param command: A command line.
        :type command: str | None
        """
        super(JobSubmission, self).__init__()

        self.executable_id = executable_id
        self.name = name
        self.command = command

        return


class JobSubmissionAdaptor(DatabaseAdaptor):
    """The :py:class:`bsf.database.JobSubmissionAdaptor` class provides database access for
    the :py:class:`bsf.database.JobSubmission` class.
    """

    def __init__(
            self,
            database_connection):
        """Initialise a :py:class:`bsf.database.JobSubmissionAdaptor` object.

        :param database_connection: A :py:class:`bsf.database.DatabaseConnection` object.
        :type database_connection: DatabaseConnection
        """
        super(JobSubmissionAdaptor, self).__init__(
            database_connection=database_connection,
            object_type=JobSubmission,
            table_name='executable',
            column_definition=[
                # Primary key
                ('executable_id', 'INTEGER PRIMARY KEY ASC AUTOINCREMENT'),
                # Name
                ('name', 'TEXT UNIQUE'),
                # Command as submitted into the Stage
                ('command', 'TEXT'),
            ])

        return

    def select_by_name(self, name):
        """Select one :py:class:`bsf.database.JobSubmission` object by name.

        :param name: A name.
        :type name: str
        :return: A :py:class:`bsf.database.JobSubmission` or :py:const:`None`.
        :rtype: JobSubmission | None
        """
        object_list = self.select(statement=self.statement_select(where_clause='name = ?'), parameters=[name])
        object_length = len(object_list)

        if object_length > 1:
            raise Exception("SQL database returned more than one row for unique field 'name'.")
        elif object_length == 1:
            return object_list[0]
        else:
            return
