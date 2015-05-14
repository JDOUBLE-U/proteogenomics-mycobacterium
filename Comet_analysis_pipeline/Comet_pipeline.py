import sys
import os

from Comet_analysis_pipeline import preprocess_db
from Comet_20151 import run_comet as comet
from Comet_analysis_pipeline import metap_motif_analysis_pipeline
from msconvert import run_msconvert as msconvert
from xinteract import run_xinteract as xinteact
from EMBOSS_sixpack_65 import run_sixpack as sixpack


__author__ = 'Jeroen'


def clear_previous_results():
    out_folders = [out_folder for out_folder in os.listdir('../' + 'Comet_analysis_pipeline')
                   if out_folder.endswith('_out')]
    for out_folder in out_folders:
        files = os.listdir(out_folder)
        for file_name in files:
            os.remove(out_folder + '/' + file_name)


def analyse_on_sixframe():
    print('Analyse on the six-frame translation of the genome? [y/n]')
    invoer = input(">>> ")
    if invoer.lower() == 'q':
        sys.exit("Tot ziens")
    elif not invoer.lower().startswith("y") and not invoer.lower().startswith("n"):
        print("The question can only be answered with y or n")
    elif invoer.lower().startswith("y"):
        return True
    else:
        return False


def get_mzxmls(mzxml_folder):
    mzxml_file_names = []
    for fileName in os.listdir(mzxml_folder):
        if fileName.endswith('.mzXML') or fileName.endswith('.mzML'):
            mzxml_file_names.append(mzxml_folder + fileName)

    return mzxml_file_names


def create_out_folder(out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    return out_folder


def main(genome_db, prot_db, spectras_folder, on_sixframe, min_pep_length):
    # Windows executables
    on_os = sys.platform
    if on_os == 'win32':
        sixpack_executable = 'cd ../EMBOSS_sixpack_65& windows_sixpack.exe'
        comet_executable = 'cd ../Comet_20151& comet.2015011.win64.exe'
        xinteract_executable = 'xinteract\\xinteract.exe'
        msconvert_executable = 'cd ../msconvert& msconvert.exe'
    # Linux and Mac executables
    else:
        sixpack_executable = './EMBOSS_sixpack_65/linux_sixpack'
        comet_executable = '../Comet_20151/comet.2015011.linux.exe'
        xinteract_executable = 'xxx'
        msconvert_executable = 'xxx'

    ## Create all output folders ##
    create_out_folder("comet_out/")
    create_out_folder("sixpack_out/")
    create_out_folder("xinteract_out/")
    create_out_folder("protxml_parsing_out/")
    create_out_folder('weblogo_out/')
    create_out_folder('preprocessed_db_out/')

    ## Clear previous results ##
    clear_previous_results()

    ## Convert any raw files to mzxml files ##
    msconvert.convert_all_raw_files(msconvert_executable, spectras_folder)

    ## SixPack ##
    if on_sixframe:
        prot_db = sixpack.run_sixpack(genome_db, sixpack_executable, on_os)

    ## Comet ##
    processed_prot_db = preprocess_db.cleave_m_only(prot_db)
    comet_pep_xmls = []
    mzxmls = get_mzxmls(spectras_folder)
    for mzxml in mzxmls:
        comet_pep_xmls.append(comet.run_comet(comet_executable, processed_prot_db, mzxml))

    ## Get ms run code ##
    first_xml = comet_pep_xmls[0]
    ms_run_code = first_xml[:-len(".pep.xml")]

    ## xinteract ##
    xinteact.run_xinteract(ms_run_code, xinteract_executable, comet_pep_xmls, min_pep_length)

    ## Find MetAp activity ##
    metap_motif_analysis_pipeline.run_metap_pipeline(ms_run_code, "..\\xinteract_out\\" + ms_run_code + '.prot.xml', prot_db)


if __name__ == '__main__':
    if analyse_on_sixframe():
        main('../' + 'GitHub_test_files/Mt_genome.fasta',
             None,
             '../' + 'GitHub_test_files/',
             True,
             5)
    else:
        main(None,
             '../' + 'GitHub_test_files/Mt_proteome.fasta',
             '../' + 'GitHub_test_files/',
             False,
             5)
