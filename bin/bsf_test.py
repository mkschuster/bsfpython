#! /usr/bin/env python
#
# BSF Python script to test library functions.
#
#
# Copyright 2013 Michael K. Schuster
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

import csv

from bsf import Executable, Runnable
from bsf.database import DatabaseConnection, ProcessSLURM, ProcessSLURMAdaptor

sacct = Executable(name='sacct', program='sacct', stdout_path='sacct_mschuster.csv', stderr_path='sacct_mschuster.err')
sacct.add_option_long(key='user', value='mschuster')
sacct.add_option_long(key='starttime', value='2014-04-19')
sacct.add_switch_long(key='long')
sacct.add_switch_long(key='parsable')

return_code = Runnable.run(executable=sacct)

Runnable.evaluate_return_code(executable=sacct, return_code=return_code)

process_slurm_list = list()

dict_file = open('sacct_mschuster.csv', 'r')
dict_reader = csv.DictReader(f=dict_file, delimiter='|')

for row in dict_reader:
    process_slurm = ProcessSLURM(
        job_id=row['JobID'],
        job_name=row['JobName'],
        partition=row['Partition'],
        max_vm_size=row['MaxVMSize'],
        max_vm_size_node=row['MaxVMSizeNode'],
        max_vm_size_task=row['MaxVMSizeTask'],
        average_vm_size=row['AveVMSize'],
        max_rss=row['MaxRSS'],
        max_rss_node=row['MaxRSSNode'],
        max_rss_task=row['MaxRSSTask'],
        average_rss=row['AveRSS'],
        max_pages=row['MaxPages'],
        max_pages_node=row['MaxPagesNode'],
        max_pages_task=row['MaxPagesTask'],
        average_pages=row['AvePages'],
        min_cpu=row['MinCPU'],
        min_cpu_node=row['MinCPUNode'],
        min_cpu_task=row['MinCPUTask'],
        average_cpu=row['AveCPU'],
        number_tasks=row['NTasks'],
        allocated_cpus=row['AllocCPUS'],
        elapsed=row['Elapsed'],
        state=row['State'],
        exit_code=row['ExitCode'],
        average_cpu_frequency=row['AveCPUFreq'],
        requested_cpu_frequency=row['ReqCPUFreq'],
        requested_memory=row['ReqMem'],
        consumed_energy=row['ConsumedEnergy'],
        max_disk_read=row['MaxDiskRead'],
        max_disk_read_node=row['MaxDiskReadNode'],
        max_disk_read_task=row['MaxDiskReadTask'],
        average_disk_read=row['AveDiskRead'],
        max_disk_write=row['MaxDiskWrite'],
        max_disk_write_node=row['MaxDiskWriteNode'],
        max_disk_write_task=row['MaxDiskWriteTask'],
        average_disk_write=row['AveDiskWrite'])

    process_slurm_list.append(process_slurm)

dict_file.close()

print 'Number of objects on process_slurm_list: {}'.format(len(process_slurm_list))

database_connection = DatabaseConnection(file_path='bsfpython_test.db')
database_connection.create_schema()

process_slurm_adaptor = ProcessSLURMAdaptor(database_connection=database_connection)

for process_slurm in process_slurm_list:
    temporary_list = process_slurm_adaptor.select_by_job_id(process_slurm.job_id)
    if not len(temporary_list):
        process_slurm_adaptor.insert(data_object=process_slurm)
    else:
        print "Stored ProcessSLURM object job_id: {} job_name: {}".format(process_slurm.job_id, process_slurm.job_name)

database_connection.connection.commit()

process_slurm_list_new = process_slurm_adaptor.select_all()

# process_slurm_list_new = process_slurm_adaptor.select_all_by_job_name(
#    name='variant_calling_align_lane_BSF_0038_C274CACXX_4_AD1')

for process_slurm in process_slurm_list_new:
    print "Job identifier: {!r} name: {!r}".format(process_slurm.job_id, process_slurm.job_name)

database_connection.connection.close()
