__author__ = 'Jeroen'

import os


def xinteract(pep_proph_exec, pep_xmls, min_pep_length):
    print("Starting xinteract")

    first_xml = pep_xmls[0]
    spectra_run_code = first_xml[first_xml.rfind("/"):first_xml.rfind(",")]

    interact_pep_xml_file = "../" + "Comet_analysis_pipeline/xinteract_out" + spectra_run_code + ".interact.prot.xml"
    # print("pep_proph_out +", interact_pep_xml_file)

    os.system(pep_proph_exec + " -p0.05"
                               " -l%i"
    # " -PPM" Accurate Mass Mode ?
                               " -Op"
                               " -N%s"
                               " %s"
              % (min_pep_length, interact_pep_xml_file, " ".join(pep_xmls)))

    return interact_pep_xml_file