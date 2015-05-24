__author__ = 'Jeroen'

import os


def get_raw_file_names(spectras_folder):
    raw_file_names = []
    for fileName in os.listdir(spectras_folder):
        if fileName.endswith('.RAW'):
            raw_file_names.append(spectras_folder + fileName)

    return raw_file_names


def convert_all_raw_files(msconvert_exec, raw_files_folder):

    raw_file_names = get_raw_file_names(raw_files_folder)

    if len(raw_file_names) > 0:
        print("Converting RAW spectra files to mzML files using msconvert.exe")
        for raw_file in get_raw_file_names(raw_files_folder):
                # Voorlopig schrijven we de mzXML weg in dezelfde folder als waar de RAW files staan
                os.system(msconvert_exec + " %s -o %s --mzML -v"
                          % (raw_file, raw_files_folder))
    else:
        print("The given folder does not contain .RAW files to convert")
