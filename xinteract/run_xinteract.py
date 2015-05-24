__author__ = 'Jeroen'

import os


def run_xinteract(ms_run_code, xinteract_exec, pep_xmls, min_pep_length, nr_threads):
    interact_absolute_out_path = "..\\Comet_analysis_pipeline\\xinteract_out"

    for nr, pep in enumerate(pep_xmls):
        pep_xmls[nr] = "..\\comet_out\\" + pep

    # xinteract exe relative to interact_absolute_out_path
    os.system("cd %s& %s -N%s.interact.pep.xml -p0.05 -l%i -OpA %s -THREADS=%i" %
              (interact_absolute_out_path, "..\\..\\" + xinteract_exec, ms_run_code, min_pep_length, " ".join(pep_xmls),
               nr_threads))
