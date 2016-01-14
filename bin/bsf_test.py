#! /usr/bin/env python
#
# BSF Python script to test library functions.
#
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

import csv

from bsf.database import DatabaseConnection, ProcessSLURM, ProcessSLURMAdaptor
from bsf.process import Executable

sacct = Executable(name='sacct', program='sacct', stdout_path='sacct_mschuster.csv', stderr_path='sacct_mschuster.err')
sacct.add_option_long(key='user', value='mschuster')
sacct.add_option_long(key='starttime', value='2014-04-19')
sacct.add_switch_long(key='long')
sacct.add_switch_long(key='parsable')

return_code = sacct.run()

sacct.evaluate_return_code(return_code=return_code)

process_slurm_list = list()

dict_file = open('sacct_mschuster.csv', 'r')
dict_reader = csv.DictReader(f=dict_file, delimiter='|')

for row_dict in dict_reader:
    process_slurm = ProcessSLURM(
            job_id=row_dict['JobID'],
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
            requested_cpu_frequency=row_dict['ReqCPUFreq'],
            requested_memory=row_dict['ReqMem'],
            consumed_energy=row_dict['ConsumedEnergy'],
            max_disk_read=row_dict['MaxDiskRead'],
            max_disk_read_node=row_dict['MaxDiskReadNode'],
            max_disk_read_task=row_dict['MaxDiskReadTask'],
            average_disk_read=row_dict['AveDiskRead'],
            max_disk_write=row_dict['MaxDiskWrite'],
            max_disk_write_node=row_dict['MaxDiskWriteNode'],
            max_disk_write_task=row_dict['MaxDiskWriteTask'],
            average_disk_write=row_dict['AveDiskWrite'])

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
