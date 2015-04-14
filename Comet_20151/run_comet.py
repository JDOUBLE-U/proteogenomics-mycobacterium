__author__ = 'Jeroen'

import os


def run_comet(comet_exec, prot_db, mzxml):
    print("Starting Comet")

    params_file = "../" + "Comet_analysis_pipeline/Comet_params/comet.params.high-low.velos-orbitrap"
    # params_file = "../" + "Comet_analysis_pipeline/Comet_params/comet.params.low-low.ion-trap"

    out_file = "../" + "Comet_analysis_pipeline/comet_out" + mzxml[mzxml.rfind("/"):-len(".mzXML")]

    os.system(comet_exec + " %s"
                           " -P%s"
                           " -N%s"
                           " -D%s"
              % (mzxml, params_file, out_file, prot_db))

    return out_file + ".pep.xml"