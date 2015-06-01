import sys
import os
import timeit

from gui import wx_gui
from analysis_pipeline import preprocess_db
from Comet_20151 import run_comet as comet
from analysis_pipeline import find_metap_activity
from msconvert import run_msconvert as msconvert
from xinteract import run_xinteract as xinteact
from EMBOSS_sixpack_65 import run_sixpack as sixpack

__author__ = 'Jeroen'


def delete_previous_results():
    out_folders = [out_folder for out_folder in os.listdir('..\\' + 'analysis_pipeline')
                   if out_folder.endswith('_out')]
    for out_folder in out_folders:
        files = os.listdir(out_folder)
        for file_name in files:
            os.remove(out_folder + '/' + file_name)


def get_mzxmls(mzxml_folder):
    mzxml_file_names = []
    for fileName in os.listdir(mzxml_folder):
        if fileName.endswith('.mzML'):
            mzxml_file_names.append(mzxml_folder + fileName)

    return mzxml_file_names


def create_out_folder(out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)


def main(clear_prev_results, ms_run_code, on_sixframe, codon_table, on_digest, spectras_folder, prot_db_path,
         genome_db_path,
         min_pep_length, min_prob, cleavage_loc, motif_range_start, motif_range_end, nr_threads):
    # Windows executables
    on_os = sys.platform
    if on_os == 'win32':
        sixpack_executable = 'cd ..\\EMBOSS_sixpack_65& windows_sixpack.exe'
        comet_executable = 'cd ..\\Comet_20151& comet.2015011.win64.exe'
        xinteract_executable = 'xinteract\\xinteract.exe'
        msconvert_executable = 'cd ..\\msconvert& msconvert.exe'
    # Linux and Mac executables
    else:
        sixpack_executable = './EMBOSS_sixpack_65/linux_sixpack'
        comet_executable = '../Comet_20151/comet.2015011.linux.exe'
        xinteract_executable = 'xxx'
        msconvert_executable = 'xxx'

    ## Create all output folders ##
    create_out_folder("comet_out\\")
    create_out_folder("sixpack_out\\")
    create_out_folder("xinteract_out\\")
    create_out_folder("protxml_parsing_out\\")
    create_out_folder('weblogo_out\\')
    create_out_folder('preprocessed_db_out\\')

    ## Clear previous results ##
    if clear_prev_results:
        print("Clearing previous results")
        delete_previous_results()

    ## Convert any raw files to mzxml files ##
    msconvert.convert_all_raw_files(msconvert_executable, spectras_folder)

    ## SixPack ##
    if on_sixframe:
        prot_db_path = sixpack.run_sixpack(genome_db_path, sixpack_executable, on_os, min_pep_length + 1, codon_table)

    ## Process protein db ##
    if on_digest:
        processed_prot_db = preprocess_db.cleave_m_and_digest(prot_db_path, min_pep_length)
    else:
        processed_prot_db = preprocess_db.cleave_m_only(prot_db_path, min_pep_length)

    ## Comet ##
    comet_pep_xmls = []
    mzxmls = get_mzxmls(spectras_folder)
    for mzxml in mzxmls:
        comet_pep_xmls.append(comet.run_comet(comet_executable, processed_prot_db, mzxml))

    ## xinteract ##
    xinteact.run_xinteract(ms_run_code, xinteract_executable, comet_pep_xmls, min_pep_length, nr_threads)

    ## Find MetAp activity ##
    find_metap_activity.run_metap_pipeline(ms_run_code,
                                           "xinteract_out\\" + ms_run_code + '.interact.prot.xml',
                                           prot_db_path, min_prob, cleavage_loc, motif_range_start,
                                           motif_range_end)


if __name__ == '__main__':
    wx_gui.run()

    start_time = timeit.default_timer()
    main(True, "M_marium_spectra", True, 11, True, "C:\\Users\\Jeroen\\Desktop\\a\\",
         '..\\' + 'GitHub_test_files\\M_marium_proteome.fasta',
         '..\\' + 'GitHub_test_files\\M_tuberculosis_genome.fasta', 7, float(0.95), 1, 1, 6, 8)
    stop_time = timeit.default_timer()
    print("The Pipeline took %i seconds to run\n" % stop_time)
