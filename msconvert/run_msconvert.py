__author__ = 'Jeroen'

import os
import subprocess
import time


def get_raw_file_locations(spectras_folder):
	raw_file_locations = []
	for file in os.listdir(spectras_folder):
		if file.endswith('.RAW'):
			raw_file_locations.append(spectras_folder + file)

	return raw_file_locations


def delete_RAW_files(spectras_folder):
	for RAW_file in get_raw_file_locations(spectras_folder):
		os.remove(RAW_file)


def check_processes(proclist):
	procstatus = [process.poll() for process in proclist]
	if procstatus.count(None) < len(procstatus):
		for status in procstatus:
			if status is not None:
				# print(proclist)
				proclist.pop(procstatus.index(status))
				# print(procstatus)
				# print(proclist)

	return proclist


def convert_all_raw_files(msconvert_exec, spectras_folder, nr_threads):
	RAW_files = get_raw_file_locations(spectras_folder)
	proclist = []

	if len(RAW_files) > 0:
		# print("Converting RAW spectra files to mzXML files using msconvert.exe")
		while len(RAW_files) > 0:
			while len(proclist) < nr_threads and len(RAW_files) > 0:
				proclist = check_processes(proclist)
				if len(RAW_files) > 0:
					new_raw_file = RAW_files.pop()

					new_process = subprocess.Popen(msconvert_exec + " %s -o %s --mzML -v" % (new_raw_file, spectras_folder))
					proclist.append(new_process)

			proclist = check_processes(proclist)


def main(msconvert_exec, spectras_folder, nr_threads):
	convert_all_raw_files(msconvert_exec, spectras_folder, nr_threads)
	# delete_RAW_files(spectras_folder)
