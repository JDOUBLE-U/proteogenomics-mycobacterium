__author__ = 'Jeroen'

import os


def run_xinteract(ms_run_code, pep_proph_exec, pep_xmls, min_pep_length):
    print("Starting xinteract")

    interact_absolute_out_path = os.getcwd() + "\\xinteract_out\\" + ms_run_code + ".pep.xml"

    os.system(pep_proph_exec + " -THREADS=7"
                               " -OpA"
                               " -l%i"
                               " -N%s"
                               " %s"
              % (min_pep_length, interact_absolute_out_path, " ".join(pep_xmls)))

    # glitch of xinteract, doesnt save in upper-case
    return interact_absolute_out_path.lower()