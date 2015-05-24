__author__ = 'Jeroen'

import os


def run_xinteract(ms_run_code, xinteract_exec, comet_pep_xmls, min_pep_length, nr_threads):
    interact_absolute_out_path = "..\\analysis_pipeline\\xinteract_out"

    for nr, pep in enumerate(comet_pep_xmls):
        comet_pep_xmls[nr] = "..\\comet_out\\" + pep

    # xinteract exe relative to interact_absolute_out_path
    os.system("cd %s& %s -N%s.interact.pep.xml -p0.05 -l%i -OpA %s -THREADS=%i" %
              (interact_absolute_out_path, "..\\..\\" + xinteract_exec, ms_run_code, min_pep_length,
               " ".join(comet_pep_xmls), nr_threads))
