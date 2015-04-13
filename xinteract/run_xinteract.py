__author__ = 'Jeroen'

import os


def xinteract(pep_proph_exec, pep_xmls, min_pep_length):
    print("Starting xinteract")

    first_xml = pep_xmls[0]
    ms_run_code = first_xml[first_xml.rfind("/") + 1:first_xml.rfind(",")]

    interact_absolute_out_path = os.getcwd() + "\\xinteract_out\\" + ms_run_code + ".pep.xml"

    os.system(pep_proph_exec + " -THREADS=7"
                               " -OpA"
                               " -l%i"
                               " -N%s"
                               " %s"
              % (min_pep_length, interact_absolute_out_path, " ".join(pep_xmls)))

    return interact_absolute_out_path