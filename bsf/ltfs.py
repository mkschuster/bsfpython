# -*- coding: utf-8 -*-
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
"""The :py:mod:`bsf.ltfs` module provides classes supporting the
`Linear Tape File System (LTFS) <https://en.wikipedia.org/wiki/Linear_Tape_File_System>`_.

This implementation is based on `IBM Spectrum® Archive Single Drive Edition 2.4.5.1
<https://www.ibm.com/docs/en/spectrum-archive-sde/2.4.5>`_ and the
`IBM Spectrum® Archive Library Edition 2.4.5.1 <https://www.ibm.com/docs/en/spectrum-archive-le/2.4.5>`_,
which are based on
`LTFS version 2.4.0 <https://www.snia.org/tech_activities/standards/curr_standards/ltfs>`_ standardised by the
`Storage Networking Industry Association (SNIA) <https://www.snia.org>`_.
"""
import errno
import json
import logging
import os
import platform
import re
import shutil
import stat
import sys
from argparse import ArgumentParser
from datetime import datetime, timezone
from tempfile import NamedTemporaryFile
from typing import Optional, TextIO, TypeVar
from xml.etree.ElementTree import ElementTree, Element

from bsf.connector import ConnectorTextIO, StandardOutputRedirection
from bsf.md5sum import MD5Sum, MD5SumArchive
from bsf.process import Executable
from bsf.standards import LinearTapeFileSystem

LTFSObjectType = TypeVar(name='LTFSObjectType', bound='LTFSObject')
LTFSCollectionType = TypeVar(name='LTFSCollectionType', bound='LTFSCollection')
LTFSArchiveType = TypeVar(name='LTFSArchiveType', bound='LTFSArchive')

module_logger = logging.getLogger(name=__name__)

# https://www.ibm.com/docs/en/spectrum-archive-sde/2.4.5?topic=attributes-virtual-extended

_attribute_dict: dict[str, dict[str, bool]] = {
    # LTFS Format 2.4.0 Specification
    'software': {
        # C.1 Software Metadata
        # Available for files and directories, in read only mode.
        'ltfs.softwareProduct': True,
        'ltfs.softwareVendor': True,  # Read only
        'ltfs.softwareVersion': True,
        'ltfs.softwareFormatSpec': True,
    },
    'drive': {
        # C.2 Drive Metadata
        # Available only at the root directory.
        'ltfs.driveEncryptionState': True,
        'ltfs.driveEncryptionMethod': True,
        'ltfs.driveCaptureDump': False,
    },
    'object': {
        # C.3 Object Metadata
        'ltfs.accessTime': True,
        'ltfs.backupTime': True,
        'ltfs.changeTime': True,
        'ltfs.createTime': True,
        'ltfs.fileUID': True,  # Read Only
        'ltfs.modifyTime': True,
        'ltfs.partition': True,  # Read Only
        'ltfs.startblock': True,  # Read Only
        # F.1 Spanning Files across Multiple Tape Volumes in LTFS
        'ltfs.spannedFileOffset': False,  # For files spanning more than one tape.
        # LTFS Format 2.4.0 Specification
        # F.2 File Permissions in LTFS
        'ltfs.permissions.unix': False,  # UNIX only
        'ltfs.permissions.posixacl': False,  # POSIX ACL
        'ltfs.permissions.ntfsacl': False,  # NTFS ACL
        'ltfs.permissions.nfsv4acl': False,  # NFS v4 ACL
        # F.3 Storing File Hash Values in LTFS
        'ltfs.hash.crc32sum': False,
        'ltfs.hash.md5sum': True,
        'ltfs.hash.sha1sum': False,
        'ltfs.hash.sha256sum': False,
        'ltfs.hash.sha512sum': False,
    },
    'volume': {
        # C.4 Volume Metadata
        'ltfs.commitMessage': True,
        # 'ltfs.indexVersion': True,
        'ltfs.indexCreator': True,  # Read Only
        'ltfs.indexGeneration': True,  # Read Only
        'ltfs.indexLocation': True,  # Read Only
        'ltfs.indexPrevious': True,  # Read Only
        'ltfs.indexTime': True,  # Read Only
        # 'ltfs.labelVersion': True,
        'ltfs.labelCreator': True,  # Read Only
        'ltfs.mamApplicationFormatVersion': True,  # Read Only
        'ltfs.mamApplicationVendor': True,  # Read Only
        'ltfs.mamApplicationVersion': True,  # Read Only
        'ltfs.mamBarcode': True,
        'ltfs.partitionMap': True,  # Read Only
        'ltfs.policyAllowUpdate': True,  # Read Only
        'ltfs.policyExists': True,  # Read Only
        'ltfs.policyMaxFileSize': True,  # Read Only
        'ltfs.sync': False,  # Write Only
        'ltfs.volumeBlocksize': True,  # Read Only
        'ltfs.volumeCompression': True,  # Read Only
        'ltfs.volumeFormatTime': True,  # Read Only
        'ltfs.volumeLockState': True,
        'ltfs.volumeName': True,
        'ltfs.volumeSerial': True,  # Read Only
        'ltfs.volumeUUID': True,  # Read Only
    },
    'media': {
        # C.5 Media Metadata
        # Available only at the root directory.
        'ltfs.mediaBeginningMediumPasses': True,
        'ltfs.mediaDataPartitionAvailableSpace': True,
        'ltfs.mediaDataPartitionTotalCapacity': True,
        'ltfs.mediaDatasetsRead': True,
        'ltfs.mediaDatasetsWritten': True,
        'ltfs.mediaEfficiency': True,
        'ltfs.mediaIndexPartitionAvailableSpace': True,
        'ltfs.mediaIndexPartitionTotalCapacity': True,
        'ltfs.mediaEncrypted': True,
        'ltfs.mediaLoads': True,
        'ltfs.mediaMBRead': True,
        'ltfs.mediaMBWritten': True,
        'ltfs.mediaMiddleMediumPasses': True,
        'ltfs.mediaPermanentReadErrors': True,
        'ltfs.mediaPermanentWriteErrors': True,
        'ltfs.mediaPool.additionalInfo': False,
        'ltfs.mediaPool.name': False,
        'ltfs.mediaPreviousPermanentReadErrors': True,
        'ltfs.mediaPreviousPermanentWriteErrors': True,
        'ltfs.mediaRecoveredReadErrors': True,
        'ltfs.mediaRecoveredWriteErrors': True,
        'ltfs.mediaStorageAlert': True,
    },
}

_cartridge_dict: dict[str, int] = {
    # IBM LTO Ultrium Cartridge Label Specification (Revision 6)
    # Part Number 19P0034
    # EC - M10321
    # http://www-01.ibm.com/support/docview.wss?uid=ssg1S7000429
    # https://www.ibm.com/docs/en/ts4300-tape-library?topic=overview-supported-tape-cartridges
    # 'L1': # Generation 1 Type A (100 GB)
    # 'LA': # Generation 1 Type B (50 GB)
    # 'LB': # Generation 1 Type C (30 GB)
    # 'LC': # Generation 1 Type D (10 GB)
    # 'L2': # Generation 2 Type A (200 GB)
    # 'L3': # Generation 3 Type A (400 GB)
    # 'LT': # Generation 3 WORM   (400 GB)
    # 'L4': # Generation 4 Type A (600 GB)
    # 'LU': # Generation 4 WORM   (600 GB)
    'L5': 1358985 * 1024 * 1024,  # Generation 5 Type A (1500 GB)
    'LV': 1358985 * 1024 * 1024,  # Generation 5 WORM   (1500 GB)
    'L6': 2296532 * 1024 * 1024,  # Generation 6 Type A (2500 GB)
    'LW': 2296532 * 1024 * 1024,  # Generation 6 WORM   (2500 GB)
    # 'L7': # Generation 7 Type A (6 TB)
    # 'LX': # Generation 7 WORM   (6 TB)
    # 'M8': # Generation M8       (9 TB)
    'L8': 11168993 * 1024 * 1024,  # Generation 8 Type A (12 TB)
    'LY': 11168993 * 1024 * 1024,  # Generation 8 WORM   (12 TB)
    # 'L9': # Generation 9 Type A (18 TB)
    # 'LZ': # Generation 9 WORM   (18 TB)

    # LTO-5
    # Virtual Extended Attribute "ltfs.mediaDataPartitionTotalCapacity": "1358985"
    # Virtual Extended Attribute "ltfs.mediaIndexPartitionTotalCapacity": "36239"
    # LTO-6
    # Virtual Extended Attribute "ltfs.mediaDataPartitionTotalCapacity": "2296532"
    # Virtual Extended Attribute "ltfs.mediaIndexPartitionTotalCapacity": "35060"
    # LTO-8
    # Virtual Extended Attribute "ltfs.mediaDataPartitionTotalCapacity": "11168993"
    # Virtual Extended Attribute "ltfs.mediaIndexPartitionTotalCapacity": "110038"
    #
    #  1,358,985 * 1,024 =  1,391,600,640 (LTO-5)
    #  2,296,532 * 1,024 =  2,351,648,768 (LTO-6)
    # 11,168,993 * 1,024 = 11,437,048,832 (LTO-8)
}

_platform_system = platform.system()

_platform_prefix: Optional[str] = None
if _platform_system == 'Darwin':
    _platform_prefix = None
elif _platform_system == 'Linux':
    _platform_prefix = 'user'
else:
    raise Exception('Unsupported platform: ' + _platform_system)

_ltfs_prefix = 'ltfs'


def get_ltfs_vea_key(key: str) -> str:
    """Get a :emphasis:`Virtual Extended Attribute` (VEA) key.

    The key is split by '.' characters and the resulting list is checked for the platform and LTFS prefixes.

    :param key: A VEA key.
    :type key: str
    :return: A VEA key.
    :rtype: str
    """
    key_index = 0
    key_list = key.split(sep='.')

    if _platform_prefix:
        if _platform_prefix in key_list:
            if key_list.index(_platform_prefix) != key_index:
                raise Exception(
                    f'The platform prefix {_platform_prefix!r} needs to be at index {key_index!r} of the key {key!r}.')
        else:
            key_list.insert(key_index, _platform_prefix)

        key_index += 1

    if _ltfs_prefix in key_list:
        if key_list.index(_ltfs_prefix) != key_index:
            raise Exception(
                f'The LTFS prefix {_ltfs_prefix!r} needs to be at index {key_index!r} of the key {key!r}.')
    else:
        key_list.insert(key_index, _ltfs_prefix)

    return '.'.join(key_list)


def get_ltfs_vea_value(file_path: str, key: str) -> str:
    """Get a :emphasis:`Virtual Extended Attribute` (VEA) value for a
    :emphasis:`Linear Tape File System` (LTFS) file path and a VEA key.

    :param file_path: A file path.
    :type file_path: str
    :param key: A VEA key.
    :type key: str
    :return: A VEA value.
    :rtype: str
    """
    key = get_ltfs_vea_key(key=key)
    value_bytes = None

    try:
        value_bytes = os.getxattr(path=file_path, attribute=key, follow_symlinks=False)
    except AttributeError as exception:
        module_logger.error('Getting extended attributes is only supported on Linux systems.')
        module_logger.error(exception)
        raise exception
    except OSError as exception:
        module_logger.debug('Getting extended attribute %r for file path %r raises OSError errno: %d',
                            key, file_path, exception.errno)
        if exception.errno == errno.ENODATA:
            # 61 errno.ENODATA: The attribute is not defined, but the namespace is otherwise correct.
            module_logger.debug('Getting extended attribute %r for file path %r is not available.', key, file_path)
        elif exception.errno == errno.ENOTSUP:
            # 95 errno.ENOTSUP: The namespace is not defined. (The Lustre file system has also 'lustre' defined.)
            module_logger.error('Getting extended attribute %r for file path %r is not supported.', key, file_path)
        else:
            module_logger.error('Getting extended attribute %r for file path %r raises OSError errno: %d error: %s',
                                key, file_path, exception.errno, exception)
            raise exception
    except Exception as exception:
        module_logger.error('Check destination: %s', str(exception))
        raise exception

    if value_bytes:
        return value_bytes.decode(encoding='utf-8')
    else:
        return ''


def set_ltfs_vea_value(file_path: str, key: str, value: str) -> None:
    """Set a :emphasis:`Virtual Extended Attribute` (VEA) value for a
    :emphasis:`Linear Tape File System` (LTFS) file path and a VEA key.

    :param file_path: A file path.
    :type file_path: str
    :param key: A VEA key.
    :type key: str
    :param value: A VEA value.
    :type value: str
    """
    key = get_ltfs_vea_key(key=key)

    # The file needs to be writable by the user, otherwise an OSError [Errno 13] "Permission denied" gets raised.

    path_stat_result = os.stat(path=file_path, follow_symlinks=True)

    if not path_stat_result.st_mode & stat.S_IWUSR:
        # If the file is not writable by the user, change file permissions prior to setting the extended attribute.
        os.chmod(
            path=file_path,
            mode=path_stat_result.st_mode | stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH,
            follow_symlinks=True)

    try:
        os.setxattr(path=file_path, attribute=key, value=value.encode(encoding='utf-8'))
    except AttributeError as exception:
        module_logger.error('Setting extended attributes is only supported on Linux systems.')
        module_logger.error(exception)
    except OSError as exception:
        module_logger.debug('Setting extended attribute %r for file path %r raises OSError errno: %d',
                            key, file_path, exception.errno)
        if exception.errno == errno.ENODATA:
            # 61 errno.ENODATA: The attribute is not defined, but the namespace is otherwise correct.
            module_logger.error('Setting extended attribute %r for file path %r is not available.', key, file_path)
        elif exception.errno == errno.ENOTSUP:
            # 95 errno.ENOTSUP: The namespace is not defined. (The Lustre file system has also 'lustre' defined.)
            module_logger.error('Setting extended attribute %r for file path %r is not supported.', key, file_path)
        else:
            module_logger.error('Setting extended attribute %r for file path %r raises OSError errno: %d error: %s',
                                key, file_path, exception.errno, exception)
            raise exception
    except Exception as exception:
        module_logger.error('Check destination: %s', str(exception))
        raise exception

    if not path_stat_result.st_mode & stat.S_IWUSR:
        # If the file was not writable by the user, reset the file permissions to tne original ones.
        os.chmod(path=file_path, mode=path_stat_result.st_mode, follow_symlinks=True)

    return


def get_ltfs_domain_attributes(
        ltfs_path: str,
        vea_dict: Optional[dict[str, dict[str, str]]] = None) -> dict[str, dict[str, str]]:
    """Get all :emphasis:`Virtual Extended Attributes` for LTFS domains :literal:`software`, :literal:`drive`,
    :literal:`media` and :literal:`volume`.

    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) directory path.
    :type ltfs_path: str
    :param vea_dict: A Python :py:class:`dict` (VEA) object.
    :type vea_dict: dict[str, dict[str, str]] | None
    :return: A Python :py:class:`dict` (VEA) object.
    :rtype: dict[str, dict[str, str]]
    """
    if vea_dict is None:
        vea_dict = dict()

    for domain_key in ('software', 'drive', 'media', 'volume'):
        if domain_key not in vea_dict:
            vea_dict[domain_key] = dict()

        for attribute_key in _attribute_dict[domain_key]:
            if _attribute_dict[domain_key][attribute_key]:
                vea_dict[domain_key][attribute_key] = get_ltfs_vea_value(file_path=ltfs_path, key=attribute_key)

    return vea_dict


def get_ltfs_object_attributes(
        ltfs_path: str,
        vea_dict: Optional[dict[str, dict[str, str]]] = None) -> dict[str, dict[str, str]]:
    """Get all :emphasis:`Virtual Extended Attributes` for LTFS domain :literal:`object` for all file system objects.

    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) directory path.
    :type ltfs_path: str
    :param vea_dict: A Python :py:class:`dict` (VEA) object.
    :type vea_dict: dict[str, dict[str, str]] | None
    :return: A Python :py:class:`dict` (VEA) object.
    :rtype: dict[str, dict[str, str]]
    """
    if vea_dict is None:
        vea_dict = dict()

    if 'objects' not in vea_dict:
        vea_dict['objects'] = list()

    for directory_path, directory_list, file_list in os.walk(top=ltfs_path):
        for object_name in file_list:
            object_path = os.path.join(directory_path, object_name)
            path_stat_result = os.stat(path=object_path, follow_symlinks=True)

            # Record the relative path from the mount point, so that LTFS subdirectories would be supported.
            # Record the size in bytes via the stat system call.

            object_dict = {
                'file.path': os.path.relpath(object_path, ltfs_path),
                'stat.st_size': path_stat_result.st_size,
                # 'stat.st_atime': path_stat_result.st_atime,
                # 'stat.st_mtime': path_stat_result.st_mtime,
                # 'stat.st_ctime': path_stat_result.st_ctime,
            }

            # Record LTFS VEAs for each object.

            domain_key = 'object'
            for attribute_key in _attribute_dict[domain_key]:
                if _attribute_dict[domain_key][attribute_key]:
                    object_dict[attribute_key] = get_ltfs_vea_value(file_path=object_path, key=attribute_key)

            vea_dict['objects'].append(object_dict)

    return vea_dict


def get_ltfs_attributes(
        ltfs_path: str,
        vea_dict: Optional[dict[str, dict[str, str]]] = None) -> dict[str, dict[str, str]]:
    """Get all :emphasis:`Virtual Extended Attributes` for all LTFS domains for a file path.

    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) directory path.
    :type ltfs_path: str
    :param vea_dict: A Python :py:class:`dict` (VEA) object.
    :type vea_dict: dict[str, dict[str, str]] | None
    :return: A Python :py:class:`dict` (VEA) object.
    :rtype: dict[str, dict[str, str]]
    """
    if vea_dict is None:
        vea_dict = dict()

    vea_dict = get_ltfs_domain_attributes(ltfs_path=ltfs_path, vea_dict=vea_dict)
    vea_dict = get_ltfs_object_attributes(ltfs_path=ltfs_path, vea_dict=vea_dict)

    return vea_dict


def is_ltfs(file_path: str) -> bool:
    """Check that a file path is backed by a :emphasis:`Linear Tape File System` (LTFS)
    via the :literal:`ltfs.softwareProduct` :emphasis:`Virtual Extended Attribute`.

    :param file_path: A file path.
    :type file_path: str
    :return: True if backed by LTFS, False otherwise.
    :rtype: bool
    """
    vea_value = get_ltfs_vea_value(file_path=file_path, key='softwareProduct')

    if vea_value is None:
        return False

    return vea_value.startswith('LTFS')


def do_ltfs_drive_mount(device_selector: str, ltfs_path: str, output_text_io: TextIO) -> None:
    """Mount a :emphasis:`Linear Tape File System` (LTFS) file system in a tape drive.

    :param device_selector: A tape drive :emphasis:`serial number` or :emphasis:`device path`.
    :type device_selector: str
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param output_text_io: A Python :py:class:`io.TextIOWrapper` (or :py:class:`typing.TextIO`) object.
    :type output_text_io: TextIO
    """
    if is_ltfs(file_path=ltfs_path):
        return

    executable = Executable(
        name='ltfs',
        program='ltfs',
        stdout=ConnectorTextIO(text_io=output_text_io),
        stderr=StandardOutputRedirection())

    mount_options = LinearTapeFileSystem.get_drive_ltfs_options()

    if mount_options:
        executable.add_option_short(key='o', value=f'devname={device_selector!s},{mount_options!s}')
    else:
        executable.add_option_short(key='o', value=f'devname={device_selector!s}')

    executable.arguments.append(ltfs_path)

    exception_str_list = executable.run()

    if exception_str_list:
        raise Exception('\n'.join(exception_str_list))

    return


def do_ltfs_library_mount(device_selector: str, ltfs_path: str, output_text_io: TextIO) -> None:
    """Mount a :emphasis:`Linear Tape File System` (LTFS) file system in a tape library.

    :param device_selector: A tape library changer :emphasis:`serial number` or :emphasis:`device path`.
    :type device_selector: str
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param output_text_io: A Python :py:class:`io.TextIOWrapper` (or :py:class:`typing.TextIO`) object.
    :type output_text_io: TextIO
    """
    if is_ltfs(file_path=ltfs_path):
        return

    executable = Executable(
        name='ltfs',
        program='ltfs',
        stdout=ConnectorTextIO(text_io=output_text_io),
        stderr=StandardOutputRedirection())

    mount_options = LinearTapeFileSystem.get_library_ltfs_options()

    if mount_options:
        executable.add_option_short(key='o', value=f'changer_devname={device_selector!s},{mount_options!s}')
    else:
        executable.add_option_short(key='o', value=f'changer_devname={device_selector!s}')

    executable.arguments.append(ltfs_path)

    exception_str_list = executable.run()

    if exception_str_list:
        raise Exception('\n'.join(exception_str_list))

    return


def do_ltfs_unmount(ltfs_path: str, output_text_io: TextIO) -> None:
    """Unmount a :emphasis:`Linear Tape File System` (LTFS) file system via the
    `Filesystem in Userspace (FUSE) <https://en.wikipedia.org/wiki/Filesystem_in_Userspace>`_
    :literal:`fusermount` command.

    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param output_text_io: A Python :py:class:`io.TextIOWrapper` (or :py:class:`typing.TextIO`) object.
    :type output_text_io: TextIO
    """
    executable = Executable(
        name='fusermount',
        program='fusermount',
        stdout=ConnectorTextIO(text_io=output_text_io),
        stderr=StandardOutputRedirection())

    executable.add_switch_short(key='u')

    executable.arguments.append(ltfs_path)

    exception_str_list = executable.run()

    if exception_str_list:
        raise Exception('\n'.join(exception_str_list))

    return


def do_ltfs_drive_format(device_selector: str, cartridge_code: str, output_text_io: TextIO) -> None:
    """Make a :emphasis:`Linear Tape File System` (LTFS) in a tape drive.

    :param device_selector: A tape drive :emphasis:`serial number` or :emphasis:`device path`.
    :type device_selector: str
    :param cartridge_code: A :emphasis:`cartridge barcode` including
        a :emphasis:`volume serial number` and
        a :emphasis:`cartridge type`.
    :type cartridge_code: str
    :param output_text_io: A Python :py:class:`io.TextIOWrapper` (or :py:class:`typing.TextIO`) object.
    :type output_text_io: TextIO
    """
    executable = Executable(
        name='mkltfs',
        program='mkltfs',
        stdout=ConnectorTextIO(text_io=output_text_io),
        stderr=StandardOutputRedirection())

    executable.add_option_pair_long(key='device', value=device_selector)
    # The tape-serial is the volume serial number (i.e., 6 alphanumeric characters, e.g. AB1234).
    executable.add_option_pair_long(key='tape-serial', value=cartridge_code[:6])
    # The volume-name is the volume serial number and the cartridge type (e.g., AB1234L6)
    executable.add_option_pair_long(key='volume-name', value=cartridge_code[:8])
    executable.add_switch_long(key='force')

    exception_str_list = executable.run()

    if exception_str_list:
        raise Exception('\n'.join(exception_str_list))

    return


def get_readable_bytes(integer_bytes: int) -> tuple[float, Optional[str]]:
    """Convert an integer number of bytes into a human-readable number with an SI suffix.

    :param integer_bytes: An integer number of bytes.
    :type integer_bytes: int
    :return: A Python :py:class:`tuple` of
        Python :py:class:`float` (readable bytes) and
        Python :py:class:`str` (SI suffix) object.
    :rtype: (float, str | None)
    """
    readable_bytes = float(integer_bytes)
    si_prefix: Optional[str] = None

    for si_prefix in ('B', 'kiB', 'MiB', 'GiB', 'TiB', 'PiB', 'EiB', 'ZiB', 'YiB',):
        if readable_bytes < 1024.0:
            readable_prefix = si_prefix
            break

        readable_bytes /= 1024.0
    else:
        readable_prefix = si_prefix
        readable_bytes *= 1024.0

    return readable_bytes, readable_prefix


def get_cartridge_total_bytes(cartridge_code: str) -> int:
    """Get the total number of bytes for a :emphasis:`cartridge barcode`.

    :param cartridge_code: A :emphasis:`cartridge barcode` including
        a :emphasis:`volume serial number` and
        a :emphasis:`cartridge type`.
    :type cartridge_code: str
    :return: Total number of bytes.
    :rtype: int
    """
    if cartridge_code[-2:] in _cartridge_dict:
        return _cartridge_dict[cartridge_code[-2:]]
    else:
        raise Exception(f'The cartridge type {cartridge_code[-2:]!r} is not supported.')


class LTFSObject(object):
    """The :py:class:`bsf.ltfs.LTFSObject` represents a file system object (i.e., file or directory).

    :ivar file_path: A file path.
    :type file_path: str | None
    :ivar md5_sum: An MD5 check sum.
    :type md5_sum: str | None
    :ivar ltfs_partition: An LTFS partition.
    :type ltfs_partition: str | None
    :ivar ltfs_start_block: An LTFS start block.
    :type ltfs_start_block: int | None
    :ivar st_size: File status file size.
    :type st_size: int | None
    :ivar st_atime: File status access time.
    :type st_atime: float | None
    :ivar st_mtime: File status modification time.
    :type st_mtime: float | None
    :ivar st_ctime: File status metadata change time.
    :type st_ctime: float | None
    """

    def __init__(
            self,
            file_path: Optional[str] = None,
            md5_sum: Optional[str] = None,
            ltfs_partition: Optional[str] = None,
            ltfs_start_block: Optional[int] = None,
            st_size: Optional[int] = None,
            st_atime: Optional[float] = None,
            st_mtime: Optional[float] = None,
            st_ctime: Optional[float] = None) -> None:
        """Initialise a :py:class:`bsf.ltfs.LTFSObject` object.

        :param file_path: A file path.
        :type file_path: str | None
        :param md5_sum: An MD5 check sum.
        :type md5_sum: str | None
        :param ltfs_partition: An LTFS partition.
        :type ltfs_partition: str | None
        :param ltfs_start_block: An LTFS start block.
        :type ltfs_start_block: int | None
        :param st_size: File status file size.
        :type st_size: int | None
        :param st_atime: File status access time.
        :type st_atime: float | None
        :param st_mtime: File status modification time.
        :type st_mtime: float | None
        :param st_ctime: File status metadata change time.
        :type st_ctime: float | None
        """
        super(LTFSObject, self).__init__()

        self.file_path = file_path
        self.md5_sum = md5_sum
        self.ltfs_partition = ltfs_partition
        self.ltfs_start_block = ltfs_start_block
        self.st_size = st_size
        self.st_atime = st_atime
        self.st_mtime = st_mtime
        self.st_ctime = st_ctime

        return

    def __repr__(self) -> str:
        """Return a printable representation of a :py:class:`bsf.ltfs.LTFSObject` object.

        :return: A printable representation.
        :rtype: str
        """
        return \
            f'{self.__class__.__name__}(' \
            f'file_path={self.file_path!r}, ' \
            f'md5_sum={self.md5_sum!r}, ' \
            f'ltfs_partition={self.ltfs_partition!r}, ' \
            f'ltfs_start_block={self.ltfs_start_block!r}, ' \
            f'st_size={self.st_size!r}, ' \
            f'st_atime={self.st_atime!r}, ' \
            f'st_mtime={self.st_mtime!r}, ' \
            f'st_ctime={self.st_ctime!r})'

    def merge(self, ltfs_object: LTFSObjectType) -> None:
        """Merge an :py:class:`bsf.ltfs.LTFSObject` object into another.

        Only those instance variables that exist in the second object are merged into the first.

        :param ltfs_object: A :py:class:`bsf.ltfs.LTFSObject` object.
        :type ltfs_object: LTFSObject
        """
        if ltfs_object.file_path != self.file_path:
            raise Exception(
                f'Cannot merge LTFSObject objects with non-matching file path instance variables. '
                f'{ltfs_object.file_path!r} versus {self.file_path!r}')

        if ltfs_object.md5_sum:
            self.md5_sum = ltfs_object.md5_sum

        if ltfs_object.ltfs_partition:
            self.ltfs_partition = ltfs_object.ltfs_partition

        if ltfs_object.ltfs_start_block:
            self.ltfs_start_block = ltfs_object.ltfs_start_block

        if ltfs_object.st_size:
            self.st_size = ltfs_object.st_size

        if ltfs_object.st_atime:
            self.st_atime = ltfs_object.st_atime

        if ltfs_object.st_mtime:
            self.st_mtime = ltfs_object.st_mtime

        if ltfs_object.st_ctime:
            self.st_ctime = ltfs_object.st_ctime

        return


class LTFSCollection(object):
    """The :py:class:`bsf.ltfs.LTFSCollection` class supports a collection of :py:class:`bsf.ltfs.LTFSObject` objects
    indexed by file path.

    :ivar ltfs_object_dict: A Python :py:class:`dict` of
        Python :py:class:`str` (file path) key and
        :py:class:`bsf.ltfs.LTFSObject` value data.
    :type ltfs_object_dict: dict[str, LTFSObject]
    :ivar cartridge_code: A :emphasis:`cartridge barcode` including
        a :emphasis:`volume serial number` and
        a :emphasis:`cartridge type`.
    :type cartridge_code: str | None
    :ivar cartridge_path: A cartridge file path.
    :type cartridge_path: str | None
    :ivar directory_path: A directory path.
    :type directory_path: str | None
    """

    def __init__(self, ltfs_object_dict: Optional[dict[str, LTFSObject]] = None) -> None:
        """Initialise a py:class:`bsf.ltfs.LTFSObject` object.

        :param ltfs_object_dict: A Python :py:class:`dict` of
            Python :py:class:`str` (file path) key and
            :py:class:`bsf.ltfs.LTFSObject` value data.
        :type ltfs_object_dict: dict[str, LTFSObject] | None
        """
        super(LTFSCollection, self).__init__()

        if ltfs_object_dict is None:
            self.ltfs_object_dict = dict()
        else:
            self.ltfs_object_dict = ltfs_object_dict

        self.cartridge_code: Optional[str] = None
        self.cartridge_path: Optional[str] = None
        self.directory_path: Optional[str] = None

        return

    def __repr__(self) -> str:
        """Return a printable representation of a :py:class:`bsf.ltfs.LTFSCollection` object.

        :return: A printable representation.
        :rtype: str
        """
        return \
            f'{self.__class__.__name__}(' \
            f'ltfs_object_dict={self.ltfs_object_dict!r}, ' \
            f'cartridge_code={self.cartridge_code!r}, ' \
            f'cartridge_path={self.cartridge_path!r}, ' \
            f'directory_path={self.directory_path!r})'

    def add_ltfs_object(self, ltfs_object: LTFSObject) -> bool:
        """Add an :py:class:`bsf.ltfs.LTFSObject` to a :py:class:`bsf.ltfs.LTFSCollection` object.

        :param ltfs_object: A :py:class:`bsf.ltfs.LTFSObject`.
        :type ltfs_object: LTFSObject
        :return: :py:const:`True` if the :py:class:`bsf.ltfs.LTFSObject` was added, :py:const:`False` otherwise.
        :rtype: bool
        """
        if ltfs_object is None:
            return False

        if not ltfs_object.file_path:
            module_logger.warning('LTFSObject without file_path.')
            return False

        if ltfs_object.file_path in self.ltfs_object_dict:
            # Issue a warning, but add the LTFSObject object regardless.
            module_logger.warning('Replacing LTFSObject.file_path %r in LTFSCollection.', ltfs_object.file_path)

        self.ltfs_object_dict[ltfs_object.file_path] = ltfs_object

        return True

    def update_ltfs_object(self, ltfs_object: LTFSObject) -> bool:
        """Update an :py:class:`bsf.ltfs.LTFSObject`.

        :param ltfs_object: A :py:class:`bsf.ltfs.LTFSObject`
        :type ltfs_object: LTFSObject
        :return: :py:const:`True` upon success, :py:const:`False` otherwise
        :rtype: bool
        """
        if ltfs_object is None:
            return False

        if not ltfs_object.file_path:
            module_logger.warning('LTFSObject without file_path.')
            return False

        if ltfs_object.file_path in self.ltfs_object_dict:
            self.ltfs_object_dict[ltfs_object.file_path].merge(ltfs_object=ltfs_object)
        else:
            self.ltfs_object_dict[ltfs_object.file_path] = ltfs_object

        return True

    def verify(self) -> None:
        """Report any :py:class:`bsf.ltfs.LTFSObject` objects of a :py:class:`bsf.ltfs.LTFSCollection` object
        that do not have a file size and an MD5 sum.
        """
        for file_path in self.ltfs_object_dict:
            ltfs_object = self.ltfs_object_dict[file_path]

            if ltfs_object.st_size is None:
                module_logger.warning('MD5 checksum without file size: %r', file_path)

            if not ltfs_object.md5_sum:
                module_logger.warning('File path without MD5 checksum: %r', file_path)

        return

    def get_total_bytes(self) -> int:
        """Get the total number of bytes the :py:class:`bsf.ltfs.LTFSCollection` object represents.

        :return: Total number of bytes.
        :rtype: int
        """
        total_bytes = 0

        for ltfs_object in self.ltfs_object_dict.values():
            if ltfs_object.st_size is not None:
                total_bytes += ltfs_object.st_size

        return total_bytes

    def copy(
            self,
            target_directory_path: str,
            output_text_io: TextIO,
            file_path_pattern_list: Optional[list[str]] = None) -> None:
        """Copy all :py:class:`bsf.ltfs.LTFSObject` objects of a
        :py:class:`bsf.ltfs.LTFSCollection` object to a target directory path.

        :param target_directory_path: A target directory path.
        :type target_directory_path: str
        :param output_text_io: A Python :py:class:`io.TextIOWrapper` (or :py:class:`typing.TextIO`) object.
        :type output_text_io: TextIO
        :param file_path_pattern_list: An optional Python :py:class:`list` object of
            Python :py:class:`str` (file path pattern) objects.
        :type file_path_pattern_list: list[str] | None
        """

        def get_ltfs_file_path(_ltfs_object: LTFSObject) -> str:
            """Get the file path from a :py:class:`bsf.ltfs.LTFSObject`.

            This is a helper function for the :py:meth:`list.sort` method.

            :param _ltfs_object: A :py:class:`bsf.ltfs.LTFSObject` object.
            :type _ltfs_object: LTFSObject
            :return: A LTFS start block.
            :rtype: int
            """
            return _ltfs_object.file_path

        def get_ltfs_start_block(_ltfs_object: LTFSObject) -> int:
            """Get the LTFS start block from a :py:class:`bsf.ltfs.LTFSObject`.

            This is a helper function for the :py:meth:`list.sort` method.

            :param _ltfs_object: A :py:class:`bsf.ltfs.LTFSObject` object.
            :type _ltfs_object: LTFSObject
            :return: A LTFS start block.
            :rtype: int
            """
            return _ltfs_object.ltfs_start_block

        def copy_ltfs_object(_ltfs_object: LTFSObject) -> None:
            """Copy one :py:class:`bsf.ltfs.LTFSObject` object to the target directory path.

            :param _ltfs_object: A :py:class:`bsf.ltfs.LTFSObject` object.
            :type _ltfs_object: LTFSObject
            """
            # Check that the LTFSObject.file_path exists.

            if not os.path.isfile(_ltfs_object.file_path):
                raise Exception('File path does not exist: ' + _ltfs_object.file_path)

            target_file_path = os.path.join(target_directory_path, os.path.basename(_ltfs_object.file_path))

            # Check if the target file path already exists and has the same size.

            try:
                target_file_path_result = os.stat(target_file_path)
            except FileNotFoundError:
                target_file_path_result = None

            if target_file_path_result:
                if target_file_path_result.st_size == _ltfs_object.st_size:
                    print('Skipping:', repr(_ltfs_object.file_path), file=output_text_io)
                else:
                    raise Exception(
                        f'File {_ltfs_object.file_path!r} source ({_ltfs_object.st_size:%d}) and '
                        f'target {target_file_path_result.st_size:%d} sizes differ.')
            else:
                executable = Executable(
                    name='cp',
                    program='cp',
                    stdout=ConnectorTextIO(text_io=output_text_io),
                    stderr=StandardOutputRedirection())

                executable.add_option_pair_long(key='preserve', value='timestamps')
                executable.add_option_pair_long(key='target-directory', value=target_directory_path)
                executable.add_switch_long(key='verbose')

                # IBM recommends setting the 'sparse' option to 'never' to avoid degrading of performance
                # when reading from an LTFS.
                # https://www.ibm.com/docs/en/spectrum-archive-le/2.4.5?topic=edition-typical-user-scenario-linux-system
                # https://www.ibm.com/docs/en/spectrum-archive-sde/2.4.5?topic=systems-accessing-tape-media
                if is_ltfs(file_path=_ltfs_object.file_path):
                    executable.add_option_pair_long(key='sparse', value='never')

                executable.arguments.append(_ltfs_object.file_path)

                exception_str_list = executable.run()

                if exception_str_list:
                    raise Exception('\n'.join(exception_str_list))

            # Set the MD5 sum as virtual extended attribute or write it to a GNU md5sum file.

            if _ltfs_object.md5_sum:
                module_logger.debug('Setting MD5 checksum %r for object %r.',
                                    _ltfs_object.md5_sum,
                                    os.path.basename(_ltfs_object.file_path))

                if target_is_ltfs:
                    set_ltfs_vea_value(
                        file_path=target_file_path,
                        key='ltfs.hash.md5sum',
                        value=_ltfs_object.md5_sum)
                else:
                    md5sum_archive = MD5SumArchive(file_path=target_file_path + '.md5')
                    md5sum_archive.add_md5sum(
                        md5sum=MD5Sum(
                            file_path=os.path.basename(_ltfs_object.file_path),
                            check_sum=_ltfs_object.md5_sum,
                            check_mode='*'))
                    md5sum_archive.to_file_path()
            else:
                module_logger.debug('Not setting MD5 checksum for object %r.',
                                    os.path.basename(_ltfs_object.file_path))

            return

        target_is_ltfs = is_ltfs(file_path=target_directory_path)

        # Optionally, filter the file path list by a list of file path patterns if one was provided.

        file_path_list: list[str] = []
        if file_path_pattern_list:
            for file_path_pattern_str in file_path_pattern_list:
                re_pattern = re.compile(pattern=file_path_pattern_str)

                for file_path in self.ltfs_object_dict:
                    if re_pattern.search(string=file_path):
                        file_path_list.append(file_path)
        else:
            file_path_list.extend(self.ltfs_object_dict.keys())

        # Order LTFSObject objects by LTFS partition and LTFS start block before copying.

        ltfs_partition_dict: dict[str, list[LTFSObject]] = dict()

        for file_path in file_path_list:
            ltfs_object = self.ltfs_object_dict[file_path]

            # Since None is a valid Python dict key, LTFSObjects resulting either from an LTFS partition or not
            # can be processed along.
            if ltfs_object.ltfs_partition not in ltfs_partition_dict:
                ltfs_partition_dict[ltfs_object.ltfs_partition] = list()

            ltfs_partition_dict[ltfs_object.ltfs_partition].append(ltfs_object)

        for ltfs_partition in ltfs_partition_dict:
            ltfs_object_list = ltfs_partition_dict[ltfs_partition]

            # LTFSObjects from an LTFS partition must be sorted by their start block,
            # the others by their file path.
            if ltfs_partition is None:
                ltfs_object_list.sort(key=get_ltfs_file_path, reverse=False)
            else:
                ltfs_object_list.sort(key=get_ltfs_start_block, reverse=False)

            for ltfs_object in ltfs_object_list:
                if ltfs_partition is None:
                    module_logger.debug('Copying file path: %r', ltfs_object.file_path)
                else:
                    module_logger.debug('Copying partition: %r start block: %d file path: %r',
                                        ltfs_object.ltfs_partition,
                                        ltfs_object.ltfs_start_block,
                                        ltfs_object.file_path)

                copy_ltfs_object(_ltfs_object=ltfs_object)

        return

    @classmethod
    def from_collection_path(cls, collection_path: str) -> LTFSCollectionType:
        """Create a :py:class:`bsf.ltfs.LTFSCollection` object from a collection file path
        listing (absolute) file paths or :emphasis:`Uniform Resource Locator` entries
        obeying the :emphasis:`file` schema.

        :param collection_path: A collection file path obeying a :literal:`<cartridge barcode>.txt` schema.
        :type collection_path: str
        :return: A :py:class:`bsf.ltfs.LTFSCollection` object.
        :rtype: LTFSCollection
        """
        if not collection_path:
            raise Exception(f"The 'collection_path' is required.")

        if not os.path.isfile(collection_path):
            raise Exception(f"The 'collection_path' path {collection_path!r} is not a regular file.")

        # Check that the collection name obeys the cartridge barcode schema.
        re_pattern = re.compile(pattern='[0-9A-Z]{6,6}L[5-9V-Z].txt')
        re_match = re.search(pattern=re_pattern, string=os.path.basename(collection_path))
        if re_match is None:
            raise Exception(f'The collection name {os.path.basename(collection_path)!r} does not obey the '
                            f'cartridge barcode pattern {re_pattern.pattern!r}.')

        ltfs_collection = cls()

        ltfs_collection.cartridge_path = collection_path

        ltfs_collection.cartridge_code = os.path.basename(ltfs_collection.cartridge_path)
        if ltfs_collection.cartridge_code.endswith('.txt'):
            # In case shell completion is used, remove the trailing *.txt suffix.
            ltfs_collection.cartridge_code = ltfs_collection.cartridge_code[:-4]

        with open(file=collection_path, mode='rt') as text_io:
            for file_path in text_io:
                file_path = file_path.strip()

                # Exclude empty lines.
                if not file_path:
                    continue

                if file_path.startswith('#'):
                    continue

                # Convert a file schema URL, which results from copy and paste in the Gnome desktop environment.
                if file_path.startswith('file://'):
                    file_path = file_path[7:]

                module_logger.debug('Processing file path: %r', file_path)

                try:
                    path_stat_result = os.stat(path=file_path, follow_symlinks=True)
                except FileNotFoundError:
                    module_logger.error('File path %r does not exist.', file_path)
                    continue

                if not stat.S_ISREG(path_stat_result.st_mode):
                    module_logger.error('File path %r is not a regular file.', file_path)
                    continue

                if file_path.endswith('.md5'):
                    module_logger.debug('Reading MD5 checksum file path: %r', file_path)

                    # The MD5SumArchive class cannot be used here, because this method has to deal with
                    # GNU md5sum files and also bare Picard MD5 sums.
                    with open(file=file_path, mode='rt') as md5_text_io:
                        for md5sum_str in md5_text_io:
                            md5_sum = MD5Sum.from_line(md5sum_str=md5sum_str.strip())

                    ltfs_collection.update_ltfs_object(
                        ltfs_object=LTFSObject(
                            # Strip the *.md5 extension to register the LTFSObject under the relevant file path.
                            file_path=file_path[:-4],
                            md5_sum=md5_sum.check_sum))
                else:
                    ltfs_collection.update_ltfs_object(
                        ltfs_object=LTFSObject(
                            file_path=file_path,
                            st_size=path_stat_result.st_size,
                            st_atime=path_stat_result.st_atime,
                            st_mtime=path_stat_result.st_mtime,
                            st_ctime=path_stat_result.st_ctime))

        ltfs_collection.verify()

        return ltfs_collection

    @classmethod
    def from_directory_path(cls, directory_path: str) -> LTFSCollectionType:
        """Create a :py:class:`bsf.ltfs.LTFSCollection` object from a directory path.

        :param directory_path: A directory path.
        :type directory_path: str
        :return: A :py:class:`bsf.ltfs.LTFSCollection` object.
        :rtype: LTFSCollection
        """
        if not directory_path:
            raise Exception("The 'directory_path' is required.")

        if not os.path.isdir(directory_path):
            raise Exception(f"The 'directory_path' {directory_path!r} is not a directory.")

        ltfs_collection = cls()

        ltfs_collection.directory_path = directory_path

        for _directory_path, directory_name_list, file_name_list in os.walk(top=directory_path):
            for file_name in file_name_list:
                file_path = os.path.join(_directory_path, file_name)

                ltfs_start_block_str = get_ltfs_vea_value(file_path=file_path, key='ltfs.startblock')
                if ltfs_start_block_str:
                    ltfs_start_block_int = int(ltfs_start_block_str)
                else:
                    ltfs_start_block_int = None

                path_stat_result = os.stat(path=file_path, follow_symlinks=True)

                ltfs_collection.add_ltfs_object(
                    ltfs_object=LTFSObject(
                        file_path=file_path,
                        md5_sum=get_ltfs_vea_value(file_path=file_path, key='ltfs.hash.md5sum'),
                        ltfs_partition=get_ltfs_vea_value(file_path=file_path, key='ltfs.partition'),
                        ltfs_start_block=ltfs_start_block_int,
                        st_size=path_stat_result.st_size,
                        st_ctime=path_stat_result.st_ctime,
                        st_mtime=path_stat_result.st_mtime,
                        st_atime=path_stat_result.st_atime))

        ltfs_collection.verify()

        return ltfs_collection

    def report_summary(
            self,
            ltfs_path: Optional[str] = None,
            output_text_io: TextIO = sys.stdout) -> None:
        """Report a size summary for a :py:class:`bsf.ltfs.LTFSCollection` object.

        Either the remaining or the exceeded space regarding the cartridge barcode or the mounted LTFS path
        gets reported.

        :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) directory path.
        :type ltfs_path: str | None
        :param output_text_io: A Python :py:class:`io.TextIOWrapper` (or :py:class:`typing.TextIO`) object.
        :type output_text_io: TextIO
        """
        if ltfs_path and is_ltfs(file_path=ltfs_path):
            # If a LTFS cartridge is already mounted, calculate sizes on the basis of virtual extended attributes.
            ltfs_block_size = int(get_ltfs_vea_value(file_path=ltfs_path, key='ltfs.volumeBlocksize'))
            ltfs_total_bytes = int(get_ltfs_vea_value(file_path=ltfs_path, key='ltfs.mediaDataPartitionTotalCapacity'))
            ltfs_free_bytes = int(get_ltfs_vea_value(file_path=ltfs_path, key='ltfs.mediaDataPartitionAvailableSpace'))
            ltfs_free_bytes *= ltfs_block_size
            ltfs_total_bytes += ltfs_block_size
        else:
            if self.cartridge_code:
                # The LTFSCollection object was instantiated from a cartridge file path.
                ltfs_total_bytes = ltfs_free_bytes = get_cartridge_total_bytes(cartridge_code=self.cartridge_code)
            else:
                raise Exception('Neither a cartridge barcode nor an LTFS path are available.')

        total_bytes = self.get_total_bytes()

        ltfs_total_readable, ltfs_total_unit = get_readable_bytes(integer_bytes=ltfs_total_bytes)
        ltfs_free_readable, ltfs_free_unit = get_readable_bytes(integer_bytes=ltfs_free_bytes)
        total_readable, total_unit = get_readable_bytes(integer_bytes=total_bytes)

        if ltfs_free_bytes <= total_bytes:
            difference_bytes = total_bytes - ltfs_free_bytes
            difference_readable, difference_unit = get_readable_bytes(integer_bytes=difference_bytes)

            print(
                f'LTFS total: {ltfs_total_bytes:,d} ({ltfs_total_readable:0.2f} {ltfs_total_unit!s})',
                f'LTFS free:  {ltfs_free_bytes:,d} ({ltfs_free_readable:0.2f} {ltfs_free_unit!s})',
                f'Total size: {total_bytes:,d} ({total_readable:0.1f} {total_unit!s})',
                f'Exceeding:  {difference_bytes:,d} ({difference_readable:0.1f} {difference_unit!s})',
                sep='\n',
                file=output_text_io)
        else:
            difference_bytes = ltfs_free_bytes - total_bytes
            difference_readable, difference_unit = get_readable_bytes(integer_bytes=difference_bytes)

            print(
                f'LTFS total: {ltfs_total_bytes:,d} ({ltfs_total_readable:0.2f} {ltfs_total_unit!s})',
                f'LTFS free:  {ltfs_free_bytes:,d} ({ltfs_free_readable:0.2f} {ltfs_free_unit!s})',
                f'Total size: {total_bytes:,d} ({total_readable:0.2f} {total_unit!s})',
                f'Remaining:  {difference_bytes:,d} ({difference_readable:0.2f} {difference_unit!s})',
                sep='\n',
                file=output_text_io)

        return


class LTFSCopyDirectory(object):
    """The :py:class:`bsf.ltfs.LTFSCopyDirectory` class represents a source (and target) directory in the copy process.

    :ivar source_path: A source directory path.
    :type source_path: str | None
    :ivar target_path: A target directory path.
    :type target_path: str | None
    :ivar source_specification: A source specification pattern including wildcard characters.
    :type source_specification: str | None
    :ivar source_file_path_list: A Python :py:class:`list` object of Python :py:class:`str` (source file path) objects.
    :type source_file_path_list: list[str]
    """

    def __init__(
            self,
            source_path: Optional[str] = None,
            target_path: Optional[str] = None,
            source_specification: Optional[str] = None,
            source_file_path_list: Optional[list[str]] = None) -> None:
        """Initialise a :py:class:`bsf.ltfs.LTFSCopyDirectory` object.

        :param source_path: A source directory path.
        :type source_path: str | None
        :param target_path: A target directory path.
        :type target_path: str | None
        :param source_specification: A source specification pattern including wildcard characters.
        :type source_specification: str | None
        :param source_file_path_list: A Python :py:class:`list` object of
            Python :py:class:`str` (source file path) objects.
        :type source_file_path_list: list[str] | None
        """
        super(LTFSCopyDirectory, self).__init__()

        self.source_path = source_path
        self.target_path = target_path
        self.source_specification = source_specification

        if source_file_path_list is None:
            self.source_file_path_list = list()
        else:
            self.source_file_path_list = source_file_path_list

        return

    def add_source_file_path(self, source_file_path: str) -> None:
        """Add a source file path to a :py:class:`bsf.ltfs.LTFSCopyDirectory` object.

        :param source_file_path: A source file path.
        :type source_file_path: str
        """
        if source_file_path is None:
            return

        if source_file_path not in self.source_file_path_list:
            self.source_file_path_list.append(source_file_path)

        return


class LTFSCopy(object):
    """The :py:class:`bsf.ltfs.LTFSCopy` class represents one :emphasis:`Linear Tape File System Copy`
    (:literal:`ltfscp`) tool process.

    :ivar total_buffer_size: A total buffer size.
    :type total_buffer_size: str | None
    :ivar buffer_size: A buffer size per thread.
    :type buffer_size: str | None
    :ivar log_level: A log level.
    :type log_level: str | None
    :ivar recursive: Recursive processing.
    :type recursive: bool | None
    :ivar sparse: Support sparse files.
    :type sparse: bool | None
    :ivar default_target_path: A default target path for all :py:class:`bsf.ltfs.LTFSCopyDirectory` objects.
    :type default_target_path: str | None
    :ivar ltfs_directory_dict: A Python :py:class:`dict` object of Python :py:class:`str` (directory path) key and
        :py:class:`bsf.ltfs.LTFSCopyDirectory` value objects.
    :type ltfs_directory_dict: dict[str, LTFSCopyDirectory]

    Definition of the `XML batch file
    <https://www.ibm.com/docs/en/spectrum-archive-le?topic=tool-use-batch-file>`_ format::

        <?xml version="1.0" encoding="UTF-8"?>
        <ltfscpspec version="1.0">
         <params>
          [parameters]
         </params>
         <data>
          <file>
           <srcpath>path</srcpath>
           <dstpath>path</dstpath>
           <srcspec>expression</srcspec>
           <sf>filename</sf>
           <sf rename="newname">filename</sf>
          </file>
          <file>
           ...
          </file>
         </data>
        </ltfscpspec>
    """

    def __init__(
            self,
            total_buffer_size: Optional[str] = None,
            buffer_size: Optional[str] = None,
            log_level: Optional[str] = None,
            recursive: Optional[bool] = None,
            sparse: Optional[bool] = None,
            default_target_path: Optional[str] = None,
            ltfs_directory_dict: Optional[dict[str, LTFSCopyDirectory]] = None) -> None:
        """Initialise a :py:class:`bsf.ltfs.LTFSCopy` object.

        :param total_buffer_size: A total buffer size.
        :type total_buffer_size: str | None
        :param buffer_size: A buffer size per thread.
        :type buffer_size: str | None
        :param log_level: A log level.
        :type log_level: str | None
        :param recursive: Recursive processing.
        :type recursive: bool | None
        :param sparse: Support sparse files.
        :type sparse: bool | None
        :param default_target_path: A default target path for all :py:class:`bsf.ltfs.LTFSCopyDirectory` objects.
        :type default_target_path: str | None
        :param ltfs_directory_dict: A Python :py:class:`dict` object of Python :py:class:`str` (directory path) key and
            :py:class:`bsf.ltfs.LTFSCopyDirectory` value objects.
        :type ltfs_directory_dict: dict[str, LTFSCopyDirectory] | None
        """
        super(LTFSCopy, self).__init__()

        self.total_buffer_size = total_buffer_size
        self.buffer_size = buffer_size
        self.log_level = log_level
        self.recursive = recursive
        self.sparse = sparse
        self.default_target_path = default_target_path

        if ltfs_directory_dict is None:
            self.ltfs_directory_dict = dict()
        else:
            self.ltfs_directory_dict = ltfs_directory_dict

        return

    def get_or_add_ltfs_directory(
            self, ltfs_directory: LTFSCopyDirectory) -> Optional[LTFSCopyDirectory]:
        """Add a :py:class:`bsf.ltfs.LTFSCopyDirectory` object to the :py:class:`bsf.ltfs.LTFSCopy` object.

        The :py:class:`bsf.ltfs.LTFSCopyDirectory` object is only added, if it does not exist already.

        :param ltfs_directory: A :py:class:`bsf.ltfs.LTFSCopyDirectory` object.
        :type ltfs_directory: LTFSCopyDirectory
        :return: A :py:class:`bsf.ltfs.LTFSCopyDirectory` object.
        :rtype: LTFSCopyDirectory | None
        """
        if ltfs_directory is None:
            return

        if ltfs_directory.source_path in self.ltfs_directory_dict:
            return self.ltfs_directory_dict[ltfs_directory.source_path]
        else:
            self.ltfs_directory_dict[ltfs_directory.source_path] = ltfs_directory
            return ltfs_directory

    def get_or_add_source_directory(
            self,
            source_path: str,
            source_specification: Optional[str] = None) -> LTFSCopyDirectory:
        """Add a source directory to a :py:class:`LTFSCopy` object.

        The corresponding :py:class:`bsf.ltfs.LTFSCopyDirectory` object is returned.

        :param source_path: A source directory path.
        :type source_path: str
        :param source_specification: A source specification.
        :type source_specification: str | None
        :return: A :py:class:`bsf.ltfs.LTFSCopyDirectory` object.
        :rtype: LTFSCopyDirectory
        """
        if source_path in self.ltfs_directory_dict:
            return self.ltfs_directory_dict[source_path]
        else:
            ltfs_directory = LTFSCopyDirectory(
                source_path=source_path,
                source_specification=source_specification)
            self.ltfs_directory_dict[ltfs_directory.source_path] = ltfs_directory
            return ltfs_directory

    def add_source_file_path(self, source_path: str) -> None:
        """Add a source file path to a :py:class:`bsf.ltfs.LTFSCopy` object.

        :param source_path: A source file path.
        :type source_path: str
        """
        source_path = os.path.normpath(source_path)
        source_directory = os.path.dirname(source_path)
        source_file_path = os.path.basename(source_path)

        ltfs_directory = self.get_or_add_source_directory(source_path=source_directory)
        ltfs_directory.add_source_file_path(source_file_path=source_file_path)

        return

    def get_element_tree(self) -> ElementTree:
        """Get the :py:class:`xml.etree.ElementTree.ElementTree` object representation of the
        :py:class:`bsf.ltfs.LTFSCopy` object.

        :return: A :py:class:`xml.etree.ElementTree.ElementTree` object.
        :rtype: ElementTree
        """
        ltfs_specification = Element(tag='ltfscpspec', attrib={'version': '1.0'})

        ltfs_element_tree = ElementTree(element=ltfs_specification)

        ltfs_parameters = Element(tag='params')
        ltfs_specification.append(ltfs_parameters)

        if self.total_buffer_size:
            param_total_buffer_size = Element(tag='total-buffer-size')
            param_total_buffer_size.text = self.total_buffer_size
            ltfs_parameters.append(param_total_buffer_size)

        if self.buffer_size:
            param_buffer_size = Element(tag='buffer-size')
            param_buffer_size.text = self.buffer_size
            ltfs_parameters.append(param_buffer_size)

        if self.log_level:
            param_log_level = Element(tag='loglevel')
            param_log_level.text = self.log_level
            ltfs_parameters.append(param_log_level)

        if self.recursive:
            param_recursive = Element(tag='recursive')
            param_recursive.text = 'enable'
            ltfs_parameters.append(param_recursive)

        if self.sparse:
            param_sparse = Element(tag='sparse')
            param_sparse.text = 'enable'
            ltfs_parameters.append(param_sparse)

        ltfs_data = Element(tag='data')
        ltfs_specification.append(ltfs_data)

        for source_path in sorted(self.ltfs_directory_dict):
            ltfs_directory = self.ltfs_directory_dict[source_path]

            ltfs_file = Element(tag='file')
            ltfs_data.append(ltfs_file)

            ltfs_file_source_path = Element(tag='srcpath')
            ltfs_file_source_path.text = ltfs_directory.source_path
            ltfs_file.append(ltfs_file_source_path)

            ltfs_file_destination_path = Element(tag='dstpath')
            if ltfs_directory.target_path:
                ltfs_file_destination_path.text = ltfs_directory.target_path
            else:
                ltfs_file_destination_path.text = self.default_target_path
            ltfs_file.append(ltfs_file_destination_path)

            if ltfs_directory.source_specification:
                ltfs_file_source_specification = Element(tag='srcspec')
                ltfs_file_source_specification.text = ltfs_directory.source_specification
                ltfs_file.append(ltfs_file_source_specification)

            ltfs_directory.source_file_path_list.sort()

            for source_file_path in ltfs_directory.source_file_path_list:
                ltfs_source_file = Element(tag='sf')
                ltfs_source_file.text = source_file_path
                ltfs_file.append(ltfs_source_file)

        return ltfs_element_tree

    def write_batch_file(self, file_path: str) -> None:
        """Write a Linear Tape File System Copy tool batch file.

        :param file_path: A file path.
        :type file_path: str
        """
        ltfs_element_tree = self.get_element_tree()
        ltfs_element_tree.write(file_or_filename=file_path)

        return


def get_ltfs_attributes_file_name(cartridge_code: str) -> str:
    """Get a :emphasis:`Virtual Extended Attribute` (VEA) JSON file name including a timestamp.

    :param cartridge_code: A :emphasis:`cartridge barcode` including
        a :emphasis:`volume serial number` and
        a :emphasis:`cartridge type`.
    :type cartridge_code: str
    :return: A :emphasis:`Virtual Extended Attribute` (VEA) JSON file name.
    :rtype: str
    """
    return f'{cartridge_code!s}_{datetime.now(tz=timezone.utc):%Y%m%d_%H%M%S}.json'


def get_ltfs_output_file_name(cartridge_code: str) -> str:
    """Get an output file name including a timestamp.

    :param cartridge_code: A :emphasis:`cartridge barcode` including
        a :emphasis:`volume serial number` and
        a :emphasis:`cartridge type`.
    :type cartridge_code: str
    :return: An output file name.
    :rtype: str
    """
    return f'{cartridge_code!s}_{datetime.now(tz=timezone.utc):%Y%m%d_%H%M%S}.txt'


def console_drive_archive(
        device_selector: str,
        ltfs_path: str,
        collection_path: str,
        json_path: Optional[str] = None,
        output_path: Optional[str] = None,
        format_cartridge: Optional[bool] = None) -> int:
    """Console function to archive files on a :emphasis:`Linear Tape File System` (LTFS) on a tape drive.

    - Check that the :literal:`device_selector` is a :emphasis:`character` device.
    - Check that the :literal:`ltfs_path` is a :emphasis:`directory`.
    - Check that the :literal:`collection_path` is a :emphasis:`regular file`.
    - Read and check all file paths of the :literal:`collection_path` option.
    - If the :literal:`output_path` was not provided, construct it from the :emphasis:`cartridge barcode`.
    - If the :literal:`json_path` was not provided, construct it from the :emphasis:`cartridge barcode`.
    - Check that the :literal:`ltfs_path` is on an already mounted LTFS.
    - If the :literal:`format_cartridge` option is set, the cartridge is formatted setting the :literal:`mkltfs`
      :literal:`--tape-serial=AB1234` and :literal:`--volume-name=AB1234L8` options so that the
      :literal:`ltfs.volumeName` virtual extended attribute (VEA) matches the
      :emphasis:`cartridge barcode` (AB1234L8).
    - Mount the LTFS.
    - Copy all file paths of the :literal:`collection_path` to the :literal:`ltfs_path`,
      but set :emphasis:`GNU md5sum` files as :literal:`ltfs.hash.md5sum` VEAs.
    - Write all virtual extended attributes (VEAs) of all LTFS objects to a JSON file.
    - Automatically unmount the LTFS.

    :param device_selector: A tape drive :emphasis:`serial number` or :emphasis:`device path`.
    :type device_selector: str
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param collection_path: A collection file path obeying a :literal:`<cartridge barcode>.txt` schema.
    :type collection_path: str
    :param json_path: A JSON (data) file path.
    :type json_path: str | None
    :param output_path: An output (log) file path.
    :type output_path: str | None
    :param format_cartridge: Request formatting the cartridge.
    :type format_cartridge: bool | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if not device_selector:
        raise Exception(f"The 'device_selector' option is required.")

    if not ltfs_path:
        raise Exception(f"The 'ltfs_path' option is required.")

    if not os.path.isdir(ltfs_path):
        raise Exception(f"The 'ltfs_path' {ltfs_path!r} is not a directory.")

    if not collection_path:
        raise Exception(f"The 'collection_path' option is required.")

    if not os.path.isfile(collection_path):
        raise Exception(f"The 'collection_path' {collection_path!r} is not a regular file.")

    module_logger.debug('Reading an LTFSCollection.')
    ltfs_collection: LTFSCollection = LTFSCollection.from_collection_path(collection_path=collection_path)

    # Create the output file name, if it was not provided already.

    if not output_path:
        output_path = get_ltfs_output_file_name(cartridge_code=ltfs_collection.cartridge_code)

    output_text_io = open(file=output_path, mode='at')

    # If the LTFS-specific virtual extended attribute 'ltfs.softwareProduct' does not exist,
    # the LTFS needs to be made and mounted before.

    if not is_ltfs(file_path=ltfs_path):
        if format_cartridge:
            # Create the LTFS if requested.

            module_logger.debug('Creating an LTFS on cartridge %r', ltfs_collection.cartridge_code)
            do_ltfs_drive_format(
                device_selector=device_selector,
                cartridge_code=ltfs_collection.cartridge_code,
                output_text_io=output_text_io)

        # Mount the LTFS.

        module_logger.debug('Mounting an LTFS on cartridge %r', ltfs_collection.cartridge_code)
        do_ltfs_drive_mount(
            device_selector=device_selector,
            ltfs_path=ltfs_path,
            output_text_io=output_text_io)

    # Since ltfs returns 0, even if LTFS mounting failed, the LTFS needs testing once more.

    if not is_ltfs(file_path=ltfs_path):
        raise Exception('Failed to mount a Linear Tape File System.')

    # Copy the source file paths onto the target file system only if they do not exist already.

    module_logger.debug('Copying files')
    ltfs_collection.copy(target_directory_path=ltfs_path, output_text_io=output_text_io)

    # Extract LTFS VEAs and store them as a JSON file.
    # Create the JSON file path, if it was not provided already.

    if not json_path:
        json_path = get_ltfs_attributes_file_name(cartridge_code=ltfs_collection.cartridge_code)

    module_logger.debug('Extracting Virtual Extended Attributes.')
    with open(file=json_path, mode='wt') as json_text_io:
        json.dump(fp=json_text_io, obj=get_ltfs_attributes(ltfs_path=ltfs_path), indent=2)

    # Unmount the LTFS.

    module_logger.debug('Unmounting the LTFS.')
    do_ltfs_unmount(ltfs_path=ltfs_path, output_text_io=output_text_io)

    print('All done.', file=output_text_io)
    output_text_io.close()

    module_logger.debug('All done.')

    return 0


def entry_point_drive_archive() -> int:
    """Console entry point to archive files on a :emphasis:`Linear Tape File System` (LTFS) on a tape drive.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Archive files on a Linear Tape File System (LTFS) in a tape drive.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--format-cartridge',
        action='store_true',
        help='format the cartridge in LTFS [False]')

    argument_parser.add_argument(
        '--device-selector',
        default=LinearTapeFileSystem.get_drive_selector(),
        help=f'tape drive serial number or device path [{LinearTapeFileSystem.get_drive_selector()}]')

    argument_parser.add_argument(
        '--json-path',
        help='JSON (data) file path [<cartridge barcode>_<data>_<time>.json]')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_drive_path(),
        help=f'LTFS mount point path [{LinearTapeFileSystem.get_drive_path()}]')

    argument_parser.add_argument(
        '--output-path',
        help='output (log) file path [<cartridge barcode>_<data>_<time>_log.txt]')

    argument_parser.add_argument(
        'collection_path',
        help='collection file path obeying a <cartridge barcode>.txt schema (e.g., AB1234L8.txt)')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_drive_archive(
        device_selector=name_space.device_selector,
        ltfs_path=name_space.ltfs_path,
        collection_path=name_space.collection_path,
        json_path=name_space.json_path,
        output_path=name_space.output_path,
        format_cartridge=name_space.format_cartridge)


def console_drive_mount(
        device_selector: str,
        ltfs_path: str,
        output_path: Optional[str] = None,
        format_cartridge: Optional[bool] = None,
        cartridge_code: Optional[str] = None) -> int:
    """Console function to mount a :emphasis:`Linear Tape File System` (LTFS) in a tape drive.

    :param device_selector: A tape drive :emphasis:`serial number` or :emphasis:`device path`.
    :type device_selector: str
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param output_path: An output (log) file path.
    :type output_path: str | None
    :param format_cartridge: Request formatting the cartridge.
    :type format_cartridge: bool | None
    :param cartridge_code: A :emphasis:`cartridge barcode` including
        a :emphasis:`volume serial number` and
        a :emphasis:`cartridge type`.
    :type cartridge_code: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if not device_selector:
        raise Exception(f"The 'device_selector' option is required.")

    if not ltfs_path:
        raise Exception(f"The 'ltfs_path' option is required.")

    if not os.path.isdir(ltfs_path):
        raise Exception(f"The 'ltfs_path' {ltfs_path!r} is not a directory.")

    if output_path:
        output_text_io = open(file=output_path, mode='wt')
    else:
        output_text_io = sys.stdout

    if format_cartridge and not cartridge_code:
        raise Exception("The 'cartridge_code' option is required if 'format_cartridge' is set.")

    # If the LTFS-specific virtual extended attribute 'ltfs.softwareProduct' does not exist,
    # the LTFS needs to be made and mounted before.

    if not is_ltfs(file_path=ltfs_path):
        if format_cartridge:
            # Create the LTFS if requested.

            module_logger.debug('Creating an LTFS on cartridge %r', cartridge_code)
            do_ltfs_drive_format(
                device_selector=device_selector,
                cartridge_code=cartridge_code,
                output_text_io=output_text_io)

        # Mount the LTFS.

        module_logger.debug('Mounting an LTFS on cartridge %r', cartridge_code)
        do_ltfs_drive_mount(
            device_selector=device_selector,
            ltfs_path=ltfs_path,
            output_text_io=output_text_io)

    # Since ltfs returns 0, even if LTFS mounting failed, the LTFS needs testing once more.

    if not is_ltfs(file_path=ltfs_path):
        raise Exception('Failed to mount a Linear Tape File System.')

    if output_path:
        output_text_io.close()

    return 0


def entry_point_drive_mount() -> int:
    """Console entry point to mount a :emphasis:`Linear Tape File System` (LTFS) in a tape drive.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Mount a Linear Tape File System (LTFS) in a tape drive.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--device-selector',
        default=LinearTapeFileSystem.get_drive_selector(),
        help=f'tape drive serial number or device path [{LinearTapeFileSystem.get_drive_selector()}]')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_drive_path(),
        help=f'LTFS mount point path [{LinearTapeFileSystem.get_drive_path()}]')

    argument_parser.add_argument(
        '--output-path',
        help='output (log) file path [sys.stout]')

    argument_parser.add_argument(
        '--format-cartridge',
        action='store_true',
        help='format the cartridge as LTFS [False]')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_drive_mount(
        device_selector=name_space.device_selector,
        ltfs_path=name_space.ltfs_path,
        output_path=name_space.output_path,
        format_cartridge=name_space.format_cartridge)


def console_drive_restore(
        device_selector: str,
        ltfs_path: str,
        target_path: str,
        json_path: Optional[str] = None,
        output_path: Optional[str] = None,
        file_path_pattern_list: Optional[list[str]] = None) -> int:
    """Console function to restore files from a :emphasis:`Linear Tape File System` (LTFS)
    in a tape drive in tape order.

    The script assumes that the LTFS volume name was set to the
    :emphasis:`cartridge barcode` including a :emphasis:`volume serial number` and a :emphasis:`cartridge type`,
    so that the :literal:`ltfs.volumeName` :emphasis:`virtual extended attribute` (VEA) is available.

    :literal:`$ mkltfs --tape-serial=AB1234 --volume-name=AB1234L8`

    - Check that the :literal:`device_selector` is a :emphasis:`character` device.
    - Check that the :literal:`ltfs_path` is a directory.
    - Check that the :literal:`target_path` is a directory.
    - Create a temporary log file unless the :literal:`output_path` option was provided.
    - Check the :literal:`ltfs_path` mount point if an LTFS is already mounted.
    - If it is not mounted already, mount the LTFS.
    - Read all LTFS virtual extended attributes (VEAs) for all objects on the tape.
    - Write all LTFS VEAs to a JSON file.
    - Filter for the :literal:`file_path_pattern` option and index by the :literal:`ltfs.partition` and
      :literal:`ltfs.startblock` VEAs.
    - Copy each matched file in tape order and convert :literal:`ltfs.hash.md5sum` VEAs into
      :emphasis:`GNU md5sum` files.
    - Unmount the LTFS.
    - If the :literal:`output_path` option was not provided, construct it from the
      :literal:`ltfs.volumeName` VEA and copy the temporary log file to the final one.

    :param device_selector: A tape drive :emphasis:`serial number` or :emphasis:`device path`.
    :type device_selector: str
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param target_path: A target directory path.
    :type target_path: str
    :param json_path: A JSON (data) file path.
    :type json_path: str | None
    :param output_path: An output (log) file path.
    :type output_path: str | None
    :param file_path_pattern_list: An optional Python :py:class:`list` object of
        Python :py:class:`str` (file path pattern) objects.
    :type file_path_pattern_list: list[str] | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if not device_selector:
        raise Exception(f"The 'device_selector' option is required.")

    if not ltfs_path:
        raise Exception(f"The 'ltfs_path' option is required.")

    if not os.path.isdir(ltfs_path):
        raise Exception(f"The 'ltfs_path' {ltfs_path!r} is not a directory.")

    if not target_path:
        raise Exception(f"The 'target_path' option is required.")

    if not os.path.isdir(target_path):
        raise Exception(f"The 'target_path' {target_path!r} is not a directory.")

    # Open the output path or create a NamedTemporaryFile object.

    if output_path:
        output_text_io = open(file=output_path, mode='at')
    else:
        output_text_io = NamedTemporaryFile(mode='wt', suffix='.txt', prefix='ltfs_bsf_', delete=False)

    # Mount the LTFS file system only if it is not already mounted.

    module_logger.debug('Mounting an LTFS.')
    do_ltfs_drive_mount(
        device_selector=device_selector,
        ltfs_path=ltfs_path,
        output_text_io=output_text_io)

    # Since ltfs returns 0, even if LTFS mounting failed, the LTFS needs testing once more.

    if not is_ltfs(file_path=ltfs_path):
        raise Exception('Failed to mount an Linear Tape File System.')

    module_logger.debug('Reading an LTFSCollection.')
    ltfs_collection: LTFSCollection = LTFSCollection.from_directory_path(directory_path=ltfs_path)

    # Copy the source file paths onto the target file system only if they do not exist already.

    module_logger.debug('Copying files')
    ltfs_collection.copy(
        target_directory_path=target_path,
        output_text_io=output_text_io,
        file_path_pattern_list=file_path_pattern_list)

    # Read all LTFS virtual extended attributes (VEAs).

    module_logger.debug('Extracting Virtual Extended Attributes.')
    vea_dict = get_ltfs_attributes(ltfs_path=ltfs_path)

    # Write the LTFS VEA to a JSON file.

    if not json_path:
        json_path = get_ltfs_attributes_file_name(cartridge_code=vea_dict['volume']['ltfs.volumeName'])

    with open(file=json_path, mode='wt') as json_text_io:
        json.dump(fp=json_text_io, obj=vea_dict, indent=2)

    # Unmount the LTFS.

    module_logger.debug('Unmounting the LTFS.')
    do_ltfs_unmount(ltfs_path=ltfs_path, output_text_io=output_text_io)

    print('All done.', file=output_text_io)
    output_text_io.close()

    # If no output path name was provided, create it on the basis of the LTFS volume name,
    # copy the temporary log file to the final one and remove the temporary file path.

    if not output_path:
        output_path = get_ltfs_output_file_name(cartridge_code=vea_dict['volume']['ltfs.volumeName'])
        shutil.copy2(src=output_text_io.name, dst=output_path)
        os.remove(output_text_io.name)

    module_logger.debug('All done.')

    return 0


def entry_point_drive_restore() -> int:
    """Console entry point to restore files from a :emphasis:`Linear Tape File System` (LTFS)
    in a tape drive in tape order.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Restore files from a Linear Tape File System (LTFS) in a tape drive.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--device-selector',
        default=LinearTapeFileSystem.get_drive_selector(),
        help=f'tape drive serial number or device path [{LinearTapeFileSystem.get_drive_selector()}]')

    argument_parser.add_argument(
        '--json-path',
        help='JSON (data) file path [<cartridge barcode>_<date>_<time>.json]')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_drive_path(),
        help=f'LTFS mount point path [{LinearTapeFileSystem.get_drive_path()}]')

    argument_parser.add_argument(
        '--output-path',
        help='output (log) file path [<cartridge barcode>_<date>_<time>_log.txt]')

    argument_parser.add_argument(
        '--target-path',
        default='.',
        help='target directory path [.]')

    argument_parser.add_argument(
        '--file-path-patterns',
        nargs='*',
        help='optional space-separated file path pattern(s)')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_drive_restore(
        device_selector=name_space.device_selector,
        ltfs_path=name_space.ltfs_path,
        target_path=name_space.target_path,
        json_path=name_space.json_path,
        output_path=name_space.output_path,
        file_path_pattern_list=name_space.file_path_patterns)


def console_drive_size(
        collection_path: str,
        mounted_cartridge: Optional[bool] = None,
        ltfs_path: Optional[str] = None) -> int:
    """Console function to summarise the sizes of files in a collection file to be written onto a
    :emphasis:`Linear Tape File System` (LTFS) in a tape drive.

    :param collection_path: A collection file path obeying a :literal:`<cartridge barcode>.txt` schema.
    :type collection_path: str
    :param mounted_cartridge: Request calculating the size on the basis of a mounted cartridge.
    :type mounted_cartridge: bool | None
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    ltfs_collection: LTFSCollection = LTFSCollection.from_collection_path(collection_path=collection_path)

    if mounted_cartridge:
        ltfs_collection.report_summary(ltfs_path=ltfs_path)
    else:
        ltfs_collection.report_summary()

    return 0


def entry_point_drive_size() -> int:
    """Console function to summarise the sizes of files in a collection file to be written onto a
    :emphasis:`Linear Tape File System` (LTFS) in a tape drive.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Linear Tape File System (LTFS) size calculation script for a tape drive.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--mounted-cartridge',
        action='store_true',
        help='calculate remaining space for a mounted cartridge')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_drive_path(),
        help=f'LTFS mount path [{LinearTapeFileSystem.get_drive_path()}] only required '
             f'for checking against a mounted cartridge')

    argument_parser.add_argument(
        'collection_path',
        help='collection file path obeying a <cartridge barcode>.txt schema (e.g., AB1234L8.txt)')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_drive_size(
        collection_path=name_space.collection_path,
        mounted_cartridge=name_space.mounted_cartridge,
        ltfs_path=name_space.ltfs_path)


def console_library_archive(
        ltfs_path: str,
        collection_path: str,
        json_path: Optional[str] = None,
        output_path: Optional[str] = None) -> int:
    """Console function to archive files on a :emphasis:`Linear Tape File System` (LTFS) in a tape library.

    - Check that the :literal:`ltfs_path` is a directory and on an LTFS.
    - Check that the :literal:`collection_path` is a regular file.
    - Read and check all file paths of the :literal:`collection_path` option.
    - Identify the cartridge path on the basis of the :emphasis:`cartridge barcode` in the :literal:`collection_path`.
    - Check that the cartridge path is a directory.
    - If the :literal:`output_path` was not provided, construct it from the :emphasis:`cartridge barcode`.
    - If the :literal:`json_path` was not provided, construct it from the :emphasis:`cartridge barcode`.
    - Copy all file paths of the :literal:`collection_path` to the cartridge path,
      but set :emphasis:`GNU md5sum` files as :literal:`ltfs.hash.md5sum` VEAs.
    - Write all virtual extended attributes (VEAs) of all LTFS objects to a JSON file.

    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param collection_path: A collection file path obeying a :literal:`<cartridge barcode>.txt` schema.
    :type collection_path: str
    :param json_path: A JSON (data) file path.
    :type json_path: str | None
    :param output_path: An output (log) file path.
    :type output_path: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if not ltfs_path:
        raise Exception(f"The 'ltfs_path' option is required.")

    if not os.path.isdir(ltfs_path):
        raise Exception(f"The 'ltfs_path' {ltfs_path!r} is not a directory.")

    if not is_ltfs(file_path=ltfs_path):
        raise Exception(f"The 'ltfs_path' {ltfs_path!r} is not on a LTFS.")

    if not collection_path:
        raise Exception(f"The 'collection_path' option is required.")

    if not os.path.isfile(collection_path):
        raise Exception(f"The 'collection_path' {collection_path!r} is not a regular file.")

    module_logger.debug('Reading an LTFSCollection.')
    ltfs_collection: LTFSCollection = LTFSCollection.from_collection_path(collection_path=collection_path)

    # The library uses <cartridge barcode>-<volume name> for tape-specific directories.
    # To support random volume names, the directory needs matching via its prefix.

    cartridge_path: Optional[str] = None
    for file_name in os.listdir(ltfs_path):
        if file_name.startswith(ltfs_collection.cartridge_code):
            cartridge_path = os.path.join(ltfs_path, file_name)
            break

    if not cartridge_path:
        raise Exception(
            f'An LTFS directory for cartridge barcode {ltfs_collection.cartridge_code!r} does not exist.')

    if not os.path.isdir(cartridge_path):
        raise Exception(f'The cartridge path {cartridge_path!r} is not a directory.')

    # Create the output file name, if it was not provided already.

    if not output_path:
        output_path = get_ltfs_output_file_name(cartridge_code=ltfs_collection.cartridge_code)

    output_text_io = open(file=output_path, mode='at')

    # Copy the source file paths onto the target file system only if they do not exist already.

    module_logger.debug('Copying files')
    ltfs_collection.copy(target_directory_path=cartridge_path, output_text_io=output_text_io)

    # Extract LTFS VEAs and store them as a JSON file.
    # Create the JSON file path, if it was not provided already.

    if not json_path:
        json_path = get_ltfs_attributes_file_name(cartridge_code=ltfs_collection.cartridge_code)

    module_logger.debug('Extracting Virtual Extended Attributes.')
    with open(file=json_path, mode='wt') as json_text_io:
        json.dump(fp=json_text_io, obj=get_ltfs_attributes(ltfs_path=cartridge_path), indent=2)

    print('All done.', file=output_text_io)
    output_text_io.close()

    module_logger.info('All done.')

    return 0


def entry_point_library_archive() -> int:
    """Console entry point to archive files on a :emphasis:`Linear Tape File System` (LTFS) in a tape library.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Archive files on a Linear Tape File System (LTFS) on a tape library.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'],
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--json-path',
        help='JSON (data) file path [<cartridge barcode>_<date>_<time>.json]')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_library_path(),
        help=f'LTFS mount point path [{LinearTapeFileSystem.get_library_path()}]')

    argument_parser.add_argument(
        '--output-path',
        help='output (log) file path [<cartridge barcode>_<date>_<time>_log.txt]')

    argument_parser.add_argument(
        'collection_path',
        help='collection file path obeying a <cartridge barcode>.txt schema (e.g., AB1234L8.txt)')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_library_archive(
        ltfs_path=name_space.ltfs_path,
        collection_path=name_space.collection_path,
        json_path=name_space.json_path,
        output_path=name_space.output_path)


def entry_point_library_mount() -> int:
    """Console entry point to mount a :emphasis:`Linear Tape File System` (LTFS) via a tape library changer.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Mount a Linear Tape File System (LTFS) in a tape library.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--device-selector',
        default=LinearTapeFileSystem.get_library_selector(),
        help=f'tape changer serial number or device path [{LinearTapeFileSystem.get_library_selector()}]')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_library_path(),
        help=f'LTFS mount point path [{LinearTapeFileSystem.get_library_path()}]')

    argument_parser.add_argument(
        '--output-path',
        help='output (log) file path [sys.stout]')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    # Since ltfs returns 0, even if LTFS mounting failed, the LTFS needs testing once more.

    do_ltfs_library_mount(
        device_selector=name_space.device_selector,
        ltfs_path=name_space.ltfs_path,
        output_text_io=sys.stdout)

    if not is_ltfs(file_path=name_space.ltfs_path):
        raise Exception('Failed to mount a Linear Tape File System.')

    return 0


def console_library_restore(
        cartridge_code: str,
        ltfs_path: str,
        target_path: str,
        json_path: Optional[str] = None,
        output_path: Optional[str] = None,
        file_path_pattern_list: Optional[list[str]] = None) -> int:
    """Console function to restore files from a :emphasis:`Linear Tape File System` (LTFS)
    in a tape library in tape order.

    - Check that the :literal:`ltfs_path` is a directory and on an LTFS.
    - Check that the :literal:`target_path` is a directory.
    - Create a temporary log file unless the :literal:`output_file` option was provided.
    - All :emphasis:`virtual extended attributes` (VEAs) of the LTFS are written to a JSON file.
    - Read all virtual extended attributes for all objects, filter for the :literal:`file_pattern` option
      and index by LTFS partition and LTFS start block.
    - Copy all files at once in tape order.
    - If the :literal:`output_path` option was not provided, construct it from the
      :literal:`ltfs.volumeName` VEA and copy the temporary log file to the final one.

    :param cartridge_code: A :emphasis:`cartridge barcode` including
        a :emphasis:`volume serial number` and
        a :emphasis:`cartridge type`.
    :type cartridge_code: str
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param target_path: A target directory path.
    :type target_path: str
    :param json_path: A JSON (data) file path.
    :type json_path: str | None
    :param output_path: An output (log) file path.
    :type output_path: str | None
    :param file_path_pattern_list: An optional Python :py:class:`list` object of
        Python :py:class:`str` (file path pattern) objects.
    :type file_path_pattern_list: list[str] | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if not ltfs_path:
        raise Exception(f"The 'ltfs_path' option is required.")

    if not os.path.isdir(ltfs_path):
        raise Exception(f"The 'ltfs_path' {ltfs_path!r} is not a directory.")

    if not is_ltfs(file_path=ltfs_path):
        raise Exception(f"The 'ltfs_path' {ltfs_path!r} is not on a LTFS.")

    if not target_path:
        raise Exception(f"The 'target_path' option is required.")

    if not os.path.isdir(target_path):
        raise Exception(f"The 'target_path' {target_path!r} is not a directory.")

    if not cartridge_code:
        raise Exception(f"The 'cartridge_code' option is required.")

    if cartridge_code.endswith('.txt'):
        cartridge_code = cartridge_code[:-4]

    # The library uses <cartridge barcode>-<volume name> for tape-specific directories.
    # To support random volume names, the directory needs matching via its prefix.

    cartridge_path: Optional[str] = None
    for file_name in os.listdir(ltfs_path):
        if file_name.startswith(cartridge_code):
            cartridge_path = os.path.join(ltfs_path, file_name)
            break

    if not cartridge_path:
        raise Exception(f'An LTFS directory for cartridge barcode {cartridge_code!r} does not exist.')

    if not os.path.isdir(cartridge_path):
        raise Exception(f'The cartridge path {cartridge_path} is not a directory.')

    # Open the output path or create a NamedTemporaryFile object.

    if output_path:
        output_text_io = open(file=output_path, mode='at')
    else:
        output_text_io = NamedTemporaryFile(mode='wt', suffix='.txt', prefix='ltfs_bsf_', delete=False)

    module_logger.debug('Reading an LTFSCollection.')
    ltfs_collection: LTFSCollection = LTFSCollection.from_directory_path(directory_path=cartridge_path)

    module_logger.debug('Copying files')
    ltfs_collection.copy(
        target_directory_path=target_path,
        output_text_io=output_text_io,
        file_path_pattern_list=file_path_pattern_list)

    # Read all LTFS virtual extended attributes (VEAs).

    module_logger.debug('Extracting Virtual Extended Attributes.')
    vea_dict = get_ltfs_attributes(ltfs_path=cartridge_path)

    # Write the LTFS VEA to a JSON file.

    if not json_path:
        json_path = get_ltfs_attributes_file_name(cartridge_code=cartridge_code)

    with open(file=json_path, mode='wt') as json_text_io:
        json.dump(fp=json_text_io, obj=vea_dict, indent=2)

    print('All done.', file=output_text_io)
    output_text_io.close()

    # If no output path name was provided, create it on the basis of the LTFS volume name,
    # copy the temporary log file to the final one and remove the temporary file path.

    if not output_path:
        output_path = get_ltfs_output_file_name(cartridge_code=cartridge_code)
        shutil.copy2(src=output_text_io.name, dst=output_path)
        os.remove(output_text_io.name)

    module_logger.info('All done.')

    return 0


def entry_point_library_restore() -> int:
    """Console entry point to restore files from a :emphasis:`Linear Tape File System` (LTFS)
    in a tape library in tape order.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Restore files from a Linear Tape File System (LTFS) in a tape library.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--json-path',
        help='JSON (data) file path [<cartridge barcode>_<date>_<time>.json]')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_library_path(),
        help=f'LTFS mount point path [{LinearTapeFileSystem.get_library_path()}]',
        dest='ltfs_path')

    argument_parser.add_argument(
        '--output-path',
        help='output (log) file path [<cartridge barcode>_<date>_<time>_log.txt]')

    argument_parser.add_argument(
        '--target-path',
        default='.',
        help='target directory path [.]')

    argument_parser.add_argument(
        '--file-path-patterns',
        nargs='*',
        help='optional space-separated file path pattern(s)')

    argument_parser.add_argument(
        'cartridge_code',
        help='cartridge barcode (e.g., AB1234L8)')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_library_restore(
        cartridge_code=name_space.cartridge_code,
        ltfs_path=name_space.ltfs_path,
        target_path=name_space.target_path,
        json_path=name_space.json_path,
        output_path=name_space.output_path,
        file_path_pattern_list=name_space.file_path_patterns)


def console_library_size(
        collection_path: str,
        mounted_cartridge: Optional[bool] = None,
        ltfs_path: Optional[str] = None) -> int:
    """Console function to summarise the sizes of files in a collection file to be written onto a
    :emphasis:`Linear Tape File System` (LTFS) in a tape library.

    :param collection_path: A collection file path obeying a :literal:`<cartridge barcode>.txt` schema.
    :type collection_path: str
    :param mounted_cartridge: Request calculating the size on the basis of a mounted cartridge.
    :type mounted_cartridge: bool | None
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    ltfs_collection: LTFSCollection = LTFSCollection.from_collection_path(collection_path=collection_path)

    if mounted_cartridge:
        # The library uses <cartridge barcode>-<volume name> for tape-specific directories.
        # To support random volume names, the directory needs matching via its prefix.

        cartridge_path: Optional[str] = None
        for file_name in os.listdir(ltfs_path):
            if file_name.startswith(ltfs_collection.cartridge_code):
                cartridge_path = os.path.join(ltfs_path, file_name)
                break

        ltfs_collection.report_summary(ltfs_path=cartridge_path)
    else:
        ltfs_collection.report_summary()

    return 0


def entry_point_library_size() -> int:
    """Console entry point to summarise the sizes of files on a cartridge list file to be written onto a
    :emphasis:`Linear Tape File System` (LTFS) in a tape library.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Linear Tape File System (LTFS) size calculation script for a tape library.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--mounted-cartridge',
        action='store_true',
        help='calculate remaining space for a mounted cartridge')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_library_path(),
        help=f'LTFS mount path [{LinearTapeFileSystem.get_library_path()}] only required '
             f'for checking against a mounted cartridge')

    argument_parser.add_argument(
        'collection_path',
        help='collection file path obeying a <cartridge barcode>.txt schema (e.g., AB1234L8.txt)')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_library_size(
        collection_path=name_space.collection_path,
        mounted_cartridge=name_space.mounted_cartridge,
        ltfs_path=name_space.ltfs_path)


def console_ltfscp(
        collection_path: str,
        target_path: str,
        total_buffer_size: str,
        buffer_size: str,
        log_level: str,
        mounted_cartridge: Optional[bool] = None,
        recursive: Optional[bool] = None,
        sparse: Optional[bool] = None) -> int:
    """Console function to prepare an XML batch file to copy files to or from a
    :emphasis:`Linear Tape File System` (LTFS) via the :emphasis:`Linear Tape File System Copy`
    (:literal:`ltfscp`) tool.

    This function requires a :emphasis:`cartridge barcode` including a :emphasis:`volume serial number` and
    a :emphasis:`cartridge type` and reads a corresponding text file with file paths from the current directory.
    The file sizes are determined and summed up to report the remaining or exceeded space depending on the
    :emphasis:`cartridge barcode` or a mounted cartridge. An XML file, which can be used as a batch file for the
    :emphasis:`Linear Tape File System Copy` (:literal:`ltfscp`) tool gets written into the current directory.

    :param collection_path: A collection file path obeying a :literal:`<cartridge barcode>.txt` schema.
    :type collection_path: str
    :param target_path: A target path.
    :type target_path: str
    :param total_buffer_size: A total buffer size.
        The default total buffer size is 400M, the maximum total buffer size is 4G.
    :type total_buffer_size: str
    :param buffer_size: A buffer size. The default buffer size is 128K, the maximum buffer size is 1G.
    :type buffer_size: str
    :param log_level: A log level.
    :type log_level: str
    :param mounted_cartridge: Request calculating remaining space for a mounted cartridge.
    :type mounted_cartridge: bool | None
    :param recursive: Request processing file system object recursively.
    :type recursive: bool | None
    :param sparse: Request not supporting sparse files.
    :type sparse: bool | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    ltfs_collection: LTFSCollection = LTFSCollection.from_collection_path(collection_path=collection_path)

    linear_tape_file_system_copy = LTFSCopy(
        total_buffer_size=total_buffer_size,
        buffer_size=buffer_size,
        log_level=log_level,
        recursive=recursive,
        sparse=sparse,
        default_target_path=target_path)

    if mounted_cartridge:
        ltfs_collection.report_summary(ltfs_path=linear_tape_file_system_copy.default_target_path)
    else:
        ltfs_collection.report_summary()

    # Transfer file paths from the LTFSCollection into tje LTFSCopy object.

    for ltfs_object in ltfs_collection.ltfs_object_dict.values():
        linear_tape_file_system_copy.add_source_file_path(source_path=ltfs_object.file_path)

    linear_tape_file_system_copy.write_batch_file(file_path=ltfs_collection.cartridge_code + '.xml')

    return 0


def entry_point_ltfscp() -> int:
    """Console entry point to prepare an XML batch file to copy files to or from a
    :emphasis:`Linear Tape File System` (LTFS) via the :emphasis:`Linear Tape File System Copy`
    (:literal:`ltfscp`) tool.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Linear Tape File System Copy tool batch file generator script.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    argument_parser.add_argument(
        '--total-buffer-size',
        default='4G',
        help='The default total buffer size is 400M. The maximum total buffer size is 4G. [4G]')

    argument_parser.add_argument(
        '--buffer-size',
        default='1G',
        help='The default buffer size is 128K. The maximum buffer size is 1G. [1G]')

    argument_parser.add_argument(
        '--no-sparse',
        action='store_false',
        help='do not support sparse files',
        dest='sparse')

    argument_parser.add_argument(
        '--recursive',
        action='store_true',
        help='process recursively')

    argument_parser.add_argument(
        '--log-level',
        default='INFO',
        choices=('URGENT', 'WARNING', 'INFO', 'DEBUG'),
        help='IBM ltfscp log level [INFO]')

    argument_parser.add_argument(
        '--target-path',
        default=LinearTapeFileSystem.get_drive_path(),
        help=f'target path [{LinearTapeFileSystem.get_drive_path()}]')

    argument_parser.add_argument(
        '--source-specification',
        default='',
        help='source specification pattern []')

    argument_parser.add_argument(
        '--mounted-cartridge',
        action='store_true',
        help='calculate remaining space for a mounted cartridge')

    argument_parser.add_argument(
        'collection_path',
        help='collection file path obeying a <cartridge barcode>.txt schema (e.g., AB1234L8.txt)')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_ltfscp(
        collection_path=name_space.collection_path,
        target_path=name_space.target_path,
        total_buffer_size=name_space.total_buffer_size,
        buffer_size=name_space.buffer_size,
        log_level=name_space.log_level,
        mounted_cartridge=name_space.mounted_cartridge,
        recursive=name_space.recursive,
        sparse=name_space.sparse)


class LTFSArchive(object):
    """The :py:class:`bsf.ltfs.LTFSArchive` class represents an archive of objects on LTFS tapes.

    :ivar archive_file_path: A :emphasis:`Linear Tape File System` archive file path.
    :type archive_file_path: str | None
    :ivar archive_dict: A Python :py:class:`dict` object of
        Python :py:class:`str` (file name) key and
        Python :py:class:`list` object of
        Python :py:class:`str` (cartridge barcode) objects.
    :type archive_dict: dict[str, list[str]]
    """

    def __init__(
            self,
            archive_file_path: Optional[str] = None,
            archive_dict: Optional[dict[str, list[str]]] = None) -> None:
        """Initialise a :py:class:`bsf.ltfs.LTFSArchive` object.

        :param archive_file_path: A :emphasis:`Linear Tape File System` archive file path.
        :type archive_file_path: str | None
        :param archive_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (file name) key and
            Python :py:class:`list` object of
            Python :py:class:`str` (cartridge barcode) objects.
        :type archive_dict: dict[str, list[str]]
        """
        self.archive_file_path = archive_file_path

        if archive_dict is None:
            self.archive_dict = {}
        else:
            self.archive_dict = archive_dict

        return

    def __repr__(self) -> str:
        """Return a printable representation of a :py:class:`bsf.ltfs.LTFSArchive` object.

        :return: A printable representation.
        :rtype: str
        """
        return \
            f'{self.__class__.__name__}(' \
            f'archive_file_path={self.archive_file_path!r}, ' \
            f'archive_dict={self.archive_dict!r})'

    @classmethod
    def from_file_path(cls, file_path: str) -> LTFSArchiveType:
        """Create a :py:class:`bsf.ltfs.LTFSArchive` object from an archive file path.

        :param file_path: A :emphasis:`Linear Tape File System` archive file path.
        :type file_path: str
        :return: A :py:class:`bsf.ltfs.LTFSArchive` object.
        :rtype: LTFSArchive
        """
        if not os.path.exists(file_path):
            raise Exception(f'File path {file_path!r} does not exist.')

        ltfs_archive = cls()

        with open(file=file_path, mode='rt') as text_io:
            for line_str in text_io:
                line_list = line_str.split()

                if line_list[0] not in ltfs_archive.archive_dict:
                    ltfs_archive.archive_dict[line_list[0]] = list()

                tape_list = ltfs_archive.archive_dict[line_list[0]]

                for tape_id in line_list[1].split(','):
                    if tape_id not in tape_list:
                        tape_list.append(tape_id)

        return ltfs_archive

    def to_file_path(
            self,
            file_path: Optional[str] = None) -> None:
        """Write a :py:class:`bsf.ltfs.LTFSArchive` object to an archive file path
        in :emphasis:`tab-separated value` (TSV) format.

        :param file_path: A :emphasis:`Linear Tape File System` archive file path
            if different from the one already stored in the :py:class:`bsf.ltfs.LTFSArchive` object.
        :type file_path: str | None
        """
        if not file_path:
            file_path = self.archive_file_path

        with open(file=file_path, mode='wt') as text_io:
            for irf_file_name in sorted(self.archive_dict):
                print(irf_file_name, ','.join(sorted(self.archive_dict[irf_file_name])), sep='\t', file=text_io)

        return

    def process_json_files(
            self,
            directory_path: str,
            file_pattern: str) -> None:
        """Process a directory of JSON files of LTFS :emphasis:`virtual extended attributes` (VEA)
        matching a file name pattern.

        :param directory_path: A directory path.
        :type directory_path: str
        :param file_pattern: A regular expression pattern matching JSON file names.
        :type file_pattern: str | None
        """
        re_pattern = re.compile(pattern=file_pattern)

        for base_path, directory_name_list, file_name_list in os.walk(top=directory_path):
            logging.debug('directory_path: %r', base_path)
            logging.debug('directory_name_list: %r', directory_name_list)
            logging.debug('file_name_list: %r', file_name_list)

            for file_name in file_name_list:
                re_match = re_pattern.search(string=file_name)
                if re_match is None:
                    logging.debug('Excluding file_name: %r', file_name)
                    continue

                # The tape name could be parsed from the file name via a regular expression, or better, retrieved from
                # the corresponding LTFS virtual extended attribute (VEA).
                # volume_name = re_match.group(1)
                with open(file=os.path.join(base_path, file_name), mode='rt') as text_io:
                    vea_dict = json.load(fp=text_io)

                    volume_name = vea_dict['volume']['ltfs.volumeName']

                    for vea_object_dict in vea_dict['objects']:
                        file_name = vea_object_dict['file.path']

                        if file_name not in self.archive_dict:
                            self.archive_dict[file_name] = list()

                        tape_list = self.archive_dict[file_name]

                        if volume_name not in tape_list:
                            tape_list.append(volume_name)

        return

    def process_ltfscp_files(
            self,
            directory_path: str,
            file_pattern: str) -> None:
        """Process a directory of :emphasis:`Linear Tape File System Copy` (:literal:`ltfscp`) tool output files,
        particularly the :literal:`ILT30505I` informational entry, that match a file name pattern.

        :param directory_path: A directory path.
        :type directory_path: str
        :param file_pattern: A regular expression pattern matching :emphasis:`Linear Tape File System Copy`
            (:literal:`ltfscp`) tool output file names.
        :type file_pattern: str
        """
        re_pattern = re.compile(pattern=file_pattern)

        for base_path, directory_name_list, file_name_list in os.walk(top=directory_path):
            logging.debug('directory_path: %r', base_path)
            logging.debug('directory_name_list: %r', directory_name_list)
            logging.debug('file_name_list: %r', file_name_list)

            for file_name in file_name_list:
                re_match = re_pattern.search(string=file_name)
                if re_match is None:
                    logging.debug('Excluding file_name: %r', file_name)
                    continue

                # Parse the volume name from the file name (e.g., BS0048L6_log.txt)
                volume_name = re_match.group(1)

                with open(file=os.path.join(base_path, file_name), mode='rt') as input_file:
                    for line_str in input_file:
                        # ILT30505I Copy /scratch/lab_bsf/archive_irf/BS0048L6/BSF_0648_H7WC2BBXY_3.bam to
                        #   /mnt/ltfs/BSF_0648_H7WC2BBXY_3.bam
                        if not line_str.startswith('ILT30505I'):
                            continue

                        line_list = line_str.split()
                        base_name = line_list[2].split('/')[-1]

                        logging.debug('File name: %r', base_name)

                        if base_name not in self.archive_dict:
                            self.archive_dict[base_name] = list()

                        tape_list = self.archive_dict[base_name]

                        if volume_name not in tape_list:
                            tape_list.append(volume_name)

        return

    def match_file_pattern(self, file_pattern: str) -> Optional[tuple[str, list[str]]]:
        """Match a file from an :py:class:`bsf.ltfs.LTFSArchive` and return the Python :py:class:`list` of
        Python :py:class:`str` (cartridge barcode) objects.

        :param file_pattern: A file name pattern.
        :type file_pattern: list[str]
        :return: A Python :py:class:`tuple` of a
            Python :py:class:`str` (file name) and a
            Python :py:class:`list` object of Python :py:class:`str` (cartridge barcode) objects.
        :rtype: tuple[str, list[str]] | None
        """
        re_pattern = re.compile(pattern=file_pattern)

        for file_name in self.archive_dict:
            re_match = re_pattern.search(string=file_name)

            if re_match is not None:
                # Find the correct tape via another match.

                return file_name, self.archive_dict[file_name]

        return None


def console_archive_update(
        directory_path: Optional[str] = None,
        archive_file_path: Optional[str] = None,
        json_pattern: Optional[str] = None,
        ltfscp_pattern: Optional[str] = None) -> int:
    """Console entry point to update a :emphasis:`Linear Tape File System` (LTFS) archive file,
    which lists archived file names and a comma-separated list of :emphasis:`cartridge barcode labels`
    in :emphasis:`tab-separated value` (TSV) format.

    A directory tree (:literal:`directory_path`) gets scanned for LTFS JSON files or
    legacy :emphasis:`Linear Tape File System Copy` (:literal:`ltfscp`) tool output files via
    regular expression patterns.
    A list of archived files is read from each file and merged into an (existing) LTFS archive file
    (:literal:`--file-path`).

    :param directory_path: A :emphasis:`Linear Tape File System` output directory path.
    :type directory_path: str | None
    :param archive_file_path: A :emphasis:`Linear Tape File System` archive file path.
    :type archive_file_path: str | None
    :param json_pattern: A regular expression pattern matching JSON file names.
    :type json_pattern: str | None
    :param ltfscp_pattern: A regular expression pattern matching :emphasis:`Linear Tape File System Copy`
        (:literal:`ltfscp`) tool output file names.
    :type ltfscp_pattern: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    ltfs_archive: LTFSArchive = LTFSArchive.from_file_path(file_path=archive_file_path)

    if json_pattern:
        ltfs_archive.process_json_files(directory_path=directory_path, file_pattern=json_pattern)

    if ltfscp_pattern:
        ltfs_archive.process_ltfscp_files(directory_path=directory_path, file_pattern=ltfscp_pattern)

    ltfs_archive.to_file_path(file_path=archive_file_path)

    return 0


def entry_point_archive_update() -> int:
    """Console entry point to update a :emphasis:`Linear Tape File System` (LTFS) archive file,
    which lists archived file names and a comma-separated list of :emphasis:`cartridge barcode labels`
    in :emphasis:`tab-separated value` (TSV) format.

    A directory tree (:literal:`directory_path`) gets scanned for LTFS JSON files or
    legacy :emphasis:`Linear Tape File System Copy` (:literal:`ltfscp`) tool output files via
    regular expression patterns.
    A list of archived files is read from each file and merged into an (existing) LTFS archive file
    (:literal:`--file-path`).

    The following defaults are used for the regular expression patterns:

    - :literal:`--json-pattern` :literal:`^([0-9A-Z]{6,6}L[5-8])_[0-9]{8,8}_[0-9]{6,6}\\.json$`
    - :literal:`--ltfscp-pattern` :literal:`^([0-9A-Z]{6,6}L[5-8])_[0-9]{8,8}_[0-9]{6,6}\\.txt$`

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Update an LTFS archive file')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    argument_parser.add_argument(
        '--json-pattern',
        default=r'^([0-9A-Z]{6,6}L[5-8])_[0-9]{8,8}_[0-9]{6,6}\.json$',
        help=r'regular expression matching JSON file names [^([0-9A-Z]{6,6}L[5-8])_[0-9]{8,8}_[0-9]{6,6}\.json$]')

    argument_parser.add_argument(
        '--ltfscp-pattern',
        default=r'^([0-9A-Z]{6,6}L[5-8])_[0-9]{8,8}_[0-9]{6,6}\.txt$',
        help=r'regular expression matching LTFSCP file names [^([0-9A-Z]{6,6}L[5-8])_[0-9]{8,8}_[0-9]{6,6}\.txt$]')

    argument_parser.add_argument(
        '--file-path',
        help='LTFS archive file path')

    argument_parser.add_argument(
        'directory_path',
        help='LTFS output directory path')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_archive_update(
        directory_path=name_space.directory_path,
        archive_file_path=name_space.file_path,
        json_pattern=name_space.json_pattern,
        ltfscp_pattern=name_space.ltfscp_pattern)


def console_library_multi_restore(
        archive_file_path: str,
        pattern_file_path: str,
        cartridge_pattern: str,
        ltfs_path: str,
        target_path: str,
        output_path: Optional[str] = None,
        dry_run: Optional[bool] = None) -> int:
    """Console function to restore multiple files from multiple :emphasis:`Linear Tape File System` (LTFS)
    cartridges in a tape library in tape order.

    An LTFS archive file of file paths and lists of cartridge barcodes  (:literal:`archive_file_path`)
    is searched with each pattern of a file of file path patterns (:literal:`archive_file_path`).
    For each matched file path the first cartridge barcode that matches the
    cartridge barcode pattern (:literal:`cartridge_pattern`) is selected.
    File paths are grouped by cartridge barcode and then retrieved from each cartridge in tape order via the
    :py:func:`bsf.ltfs.console_library_restore` function.

    :param archive_file_path: A :emphasis:`Linear Tape File System` archive file path.
    :type archive_file_path: str
    :param pattern_file_path: A file of file patterns to match from the LTFS archive.
    :type pattern_file_path: str
    :param cartridge_pattern: A cartridge barcode pattern.
    :type cartridge_pattern: str
    :param ltfs_path: A :emphasis:`Linear Tape File System` (LTFS) mount point path.
    :type ltfs_path: str
    :param target_path: A target directory path.
    :type target_path: str
    :param output_path: An output (log) file path.
    :type output_path: str | None
    :param dry_run: Request a dry-run.
    :type dry_run: bool | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    ltfs_archive: LTFSArchive = LTFSArchive.from_file_path(file_path=archive_file_path)

    re_pattern = re.compile(pattern=cartridge_pattern)

    # Create a dict with tape barcode key and file object value data.

    cartridge_dict: dict[str, list[str]] = {}

    with open(file=pattern_file_path, mode='rt') as text_io:
        for file_pattern_str in text_io:
            file_pattern_str = file_pattern_str.strip()

            logging.debug('Processing file pattern %r', file_pattern_str)

            archive_tuple = ltfs_archive.match_file_pattern(file_pattern=file_pattern_str)

            if archive_tuple is None:
                logging.info('Could not match a file name for pattern %r', file_pattern_str)
            else:
                logging.info('Matched file name %r for pattern %r', archive_tuple[0], file_pattern_str)

                for barcode in sorted(archive_tuple[1]):
                    logging.debug('Processing cartridge barcode %r', barcode)

                    re_match = re_pattern.search(string=barcode)

                    if re_match is not None:
                        logging.debug(
                            'Matched cartridge barcode %r for cartridge pattern %r',
                            barcode,
                            cartridge_pattern)

                        if barcode not in cartridge_dict:
                            cartridge_dict[barcode] = []

                        file_name_list = cartridge_dict[barcode]

                        if archive_tuple[0] not in file_name_list:
                            file_name_list.append(archive_tuple[0])

                        # Keep only the first tape barcode that matches.
                        break

    # Now, restore all files for each cartridge barcode.
    # Create one common output file for each cartridge.

    if not output_path:
        output_path = f'bsf_ltfs_multi_restore_{datetime.now(tz=timezone.utc):%Y%m%d_%H%M%S}.txt'

    for barcode in cartridge_dict:
        if dry_run:
            logging.info('Dry-run: Processing cartridge %r and files %r', barcode, cartridge_dict[barcode])
        else:
            logging.info('Processing cartridge %r', barcode)

            console_library_restore(
                cartridge_code=barcode,
                ltfs_path=ltfs_path,
                target_path=target_path,
                json_path=None,
                output_path=output_path,
                file_path_pattern_list=cartridge_dict[barcode])

    return 0


def entrypoint_library_multi_restore() -> int:
    """Console entry point to restore multiple files from multiple :emphasis:`Linear Tape File System` (LTFS)
    cartridges in a tape library in tape order.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='BSF utility script to restore multiple files matched from an LTFS archive file.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    argument_parser.add_argument(
        '--dry-run',
        action='store_true',
        help='dry run')

    argument_parser.add_argument(
        '--archive-file-path',
        help='LTFS archive file path')

    argument_parser.add_argument(
        '--pattern-file-path',
        help='file name of file path patterns')

    argument_parser.add_argument(
        '--cartridge-pattern',
        help='cartridge barcode pattern')

    argument_parser.add_argument(
        '--ltfs-path',
        default=LinearTapeFileSystem.get_library_path(),
        help=f'LTFS mount point path [{LinearTapeFileSystem.get_library_path()}]',
        dest='ltfs_path')

    argument_parser.add_argument(
        '--output-path',
        help='output (log) file path [bsf_ltfs_multi_restore_<date>_<time>_log.txt]')

    argument_parser.add_argument(
        '--target-path',
        default='.',
        help='target directory path [.]')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_library_multi_restore(
        archive_file_path=name_space.archive_file_path,
        pattern_file_path=name_space.pattern_file_path,
        cartridge_pattern=name_space.cartridge_pattern,
        ltfs_path=name_space.ltfs_path,
        target_path=name_space.target_path,
        output_path=name_space.output_path,
        dry_run=name_space.dry_run)
