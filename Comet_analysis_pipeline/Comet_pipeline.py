import sys
import os

__author__ = 'Jeroen'

from Comet_analysis_pipeline import preprocess_db
from Comet_analysis_pipeline.run_comet import run_comet
from Comet_analysis_pipeline import metap_motif_analysis_pipeline
from xinteract import run_xinteract
from EMBOSS_sixpack_65 import sixpack


def get_mzxmls(folder):
    mzxml_names = []
    for fileName in os.listdir(folder):
        if fileName.endswith(".mzXML"):
            mzxml_names.append(folder + fileName)

    return mzxml_names


def main(prot_db, mzxmls, on_sixframe, min_pep_length):
    on_platform = sys.platform
    if on_platform == 'win32':
        sixpack_executable = "cd ../EMBOSS_sixpack_65& windows_sixpack.exe"
        comet_executable = "cd ../Comet_executables& comet.2015011.win64.exe"
        xinteract_executable = "cd ../xinteract& xinteract.exe"
    else:
        sixpack_executable = "./sixpack"
        comet_executable = "./Comet_executables/comet.2015011.linux.exe"

    ## SixPack ##
    if on_sixframe:
        prot_db = sixpack.run_sixpack(prot_db, sixpack_executable, on_platform)

    ## Comet ##
    processed_prot_db = preprocess_db.cleave_M_only(prot_db)
    comet_pep_xmls = []
    for mzxml in mzxmls:
        comet_pep_xmls.append(run_comet(comet_executable, processed_prot_db, mzxml))

    ## xinteract ##
    prot_xml = run_xinteract.xinteract(xinteract_executable, comet_pep_xmls, min_pep_length)

    ## Find MetAp activity ##
    metap_motif_analysis_pipeline.main(prot_xml, prot_db)


if __name__ == '__main__':
    main("../" + "GitHub_test_files/Uniprot_Mt_proteome.fasta",
         get_mzxmls("../Local_test_files/"),
         False,
         5)