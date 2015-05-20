__author__ = 'Jeroen'

import os


def run_comet(comet_exec, prot_db, mzxml):
    # params_file = "../" + "Comet_analysis_pipeline/Comet_params/comet.params.high-low.velos-orbitrap"
    params_file = "../" + "Comet_analysis_pipeline/Comet_params/comet.params.low-low.ion-trap"

    out_name = mzxml[mzxml.rfind("/") + 1:-len(".mzXML")]

    out_file = "../" + "Comet_analysis_pipeline/comet_out/" + out_name

    os.system(comet_exec + " %s"
                           " -P%s"
                           " -N%s"
                           " -D%s"
              % (mzxml, params_file, out_file, prot_db))

    return out_name + ".pep.xml"
