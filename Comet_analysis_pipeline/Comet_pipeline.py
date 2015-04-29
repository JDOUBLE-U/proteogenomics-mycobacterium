import sys
import os

from Comet_analysis_pipeline import preprocess_db
from Comet_20151 import run_comet as comet
from Comet_analysis_pipeline import metap_motif_analysis_pipeline
from xinteract import run_xinteract as xinteact
from EMBOSS_sixpack_65 import run_sixpack as sixpack


__author__ = 'Jeroen'


def clear_previous_results():
    out_folder = [out_folder for out_folder in os.listdir('../' + 'Comet_analysis_pipeline')
                  if out_folder.endswith('_out')]
    for out_folder in out_folder:
        files = os.listdir(out_folder)
        for file_name in files:
            os.remove(out_folder + '/' + file_name)


def analyse_on_sixframe():
    # TODO make elephant-proof
    print('Search on six-frame?')
    invoer = input(">>> ")
    if invoer.lower() == 'q':
        sys.exit("Tot ziens")
    elif invoer.lower().startswith("y") == False and invoer.lower().startswith("n") == False:
        print("Antwoord alstublieft allen met de tekens y of n")
    elif invoer.lower().startswith("y"):
        return True
    else:
        return False


def get_mzxmls(mzxml_folder):
    mzxml_names = []
    for fileName in os.listdir(mzxml_folder):
        if fileName.endswith('.mzXML'):
            mzxml_names.append(mzxml_folder + fileName)

    return mzxml_names


def main(genome_db, prot_db, mzxmls, on_sixframe, min_pep_length):
    on_os = sys.platform
    if on_os == 'win32':
        sixpack_executable = 'cd ../EMBOSS_sixpack_65& windows_sixpack.exe'
        comet_executable = 'cd ../Comet_20151& comet.2015011.win64.exe'
        xinteract_executable = 'cd ../xinteract& xinteract.exe'
    else:
        sixpack_executable = './EMBOSS_sixpack_65/linux_sixpack'
        comet_executable = '../Comet_20151/comet.2015011.linux.exe'
        xinteract_executable = 'xxx'

    ## Clear previous results
    clear_previous_results()

    ## SixPack ##
    if on_sixframe:
        prot_db = sixpack.run_sixpack(genome_db, sixpack_executable, on_os)

    ## Comet ##
    processed_prot_db = preprocess_db.cleave_m_only(prot_db)
    comet_pep_xmls = []
    for mzxml in mzxmls:
        comet_pep_xmls.append(comet.run_comet(comet_executable, processed_prot_db, mzxml))

    # ## Get ms run code
    # first_xml = comet_pep_xmls[0]
    # ms_run_code = first_xml[first_xml.rfind('/') + 1:first_xml.rfind(',')].upper()
    #
    # ## xinteract ##
    # xinteract_out = xinteact.run_xinteract(ms_run_code, xinteract_executable, comet_pep_xmls, min_pep_length)
    #
    # ## Find MetAp activity ##
    # metap_motif_analysis_pipeline.run_metap_pipeline(ms_run_code, xinteract_out[:-len('pep.xml')] + 'prot.xml', prot_db)
    #

if __name__ == '__main__':
    if analyse_on_sixframe():
        main('../' + 'GitHub_test_files/Mt_genome.fasta',
             None,
             get_mzxmls('../' + 'GitHub_test_files/'),
             True,
             5)
    else:
        main(None,
             '../' + 'GitHub_test_files/Mt_proteome.fasta',
             get_mzxmls('../' + 'GitHub_test_files/'),
             False,
             5)
