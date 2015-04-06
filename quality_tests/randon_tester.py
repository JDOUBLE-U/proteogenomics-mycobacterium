# coding=utf-8

__author__ = 'Jeroen'

import radon.complexity as compl
import radon.metrics as mi

sort_attrib = "LINES"

script_loc = "../" + "Comet/motif_analysis_pipeline.py"
script_name = script_loc[script_loc.rfind("/") + 1:]

raw_code = open(script_loc).read()

MI = mi.mi_visit(raw_code, False)

sort_order = getattr(compl, sort_attrib)

print('Quality analasys of "%s":' % script_name)
print("")
print('Maintainability Index equals %f -> ranked %s' % (MI, mi.mi_rank(MI)))
print("")
# print("""
#     ============= =====================================================
#         1 - 5        A (low risk - simple block)
#         6 - 10       B (low risk - well structured and stable block)
#         11 - 20      C (moderate risk - slightly complex block)
#         21 - 30      D (more than moderate risk - more complex block)
#         31 - 40      E (high risk - complex block, alarming)
#         41+          F (very high risk - error-prone, unstable block)
#     ============= =====================================================
# """)
print("Cyclomatic Complexity score per block, sorted on %s:" % sort_attrib)
for block in compl.sorted_results(compl.cc_visit(raw_code), order=sort_order):
    print("%s -> ranked %s" % (block, compl.cc_rank(block.complexity)))


