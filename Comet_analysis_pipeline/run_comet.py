__author__ = 'Jan-Willem'

import sys
import os

# TODO make possible to chooose the type of Comet parameter file


def run_comet(protein_database, mzxml_file):
    params_file = "../" + "Comet_analysis_pipeline/Comet_params/comet.params.high-low.velos-orbitrap"
    out_file = "../" + "Comet_analysis_pipeline/comet_out" + mzxml_file[mzxml_file.rfind("/"):-len(".mzXML")]

    on_platform = sys.platform
    if on_platform == 'win32':
        comet_command = "cd ../Comet_executables& comet.2015011.win64.exe"
    else:
        comet_command = "./Comet_executables/comet.2015011.linux.exe"

    os.system(comet_command + " %s"
                              " -P%s"
                              " -N%s"
                              " -D%s"
              % (mzxml_file, params_file, out_file, protein_database))


if __name__ == '__main__':
    run_comet(protein_database="../" + "GitHub_test_files/Uniprot_Mt_proteome.fasta",
              mzxml_file="../" + "GitHub_test_files/Tuberculosis_test_spectra.mzXML")