#!/usr/bin/python
# coding=utf-8
import os

from lxml import objectify

"""
Jeroen Merks
Instructions:
- Install lxml 3.4.2: https://pypi.python.org/pypi/lxml/
...
"""


class ProtProphXMLParser():
    """
    :param parse_results_folder:
    :param prot_prophet_xml:
    """

    def __init__(self, parse_results_folder, prot_prophet_xml):
        if not os.path.exists(parse_results_folder):
            os.mkdir(parse_results_folder)

        self.parse_results_folder = parse_results_folder
        self.prot_prophet_xml = prot_prophet_xml
        self.prot_groups = self.get_prot_groups()
        self.statistics_dict = self.get_error_prob_pairs()

    def get_prot_groups(self):
        """
        Extracts groups of proteins by looping through the "protein_group"
        XML-elements of the prot.xml file.

        :return: list
        """
        prot_groups = []
        prot_prophet_xml_object = open(self.prot_prophet_xml, encoding="utf-8")
        prot_prophet_xml_object.readline()  # skip XML encoding declaration

        for prot_group in objectify.fromstring(prot_prophet_xml_object.read()).protein_group:
            prot_groups.append(ProteinGroup(prot_group))

        prot_prophet_xml_object.close()

        return prot_groups

    def get_error_prob_pairs(self):
        """
 
        :return:
        """
        error_prob_pairs = []

        prot_prophet_xml_object = open(self.prot_prophet_xml, encoding="utf-8")
        prot_prophet_xml_object.readline()  # skip XML encoding

        element_name = 'protein_summary_data_filter'
        for data_filter in objectify.fromstring(prot_prophet_xml_object.read()).protein_summary_header.program_details. \
                proteinprophet_details.getchildren():

            # We are only interested in protein_summary_data in the protein_summary_data_filter element
            if data_filter.tag[-len(element_name):] == element_name:
                min_probability = float(data_filter.attrib["min_probability"])
                error_rate = float(data_filter.attrib["false_positive_error_rate"])
                error_prob_pairs.append((error_rate, min_probability))

        prot_prophet_xml_object.close()

        return error_prob_pairs


class ProteinGroup():
    """
 
    :param prot_group:
    """

    def __init__(self, prot_group):
        self.prob = float(prot_group.attrib["probability"])
        self.group_nr = int(prot_group.attrib["group_number"])

        # Not all protein_groups have pseudo names
        try:
            self.pseudo_name = int(prot_group.attrib["pseudo_name"])
        except KeyError:
            self.pseudo_name = None

        self.prots = self.set_prots(prot_group)

    @staticmethod
    def set_prots(prot_group):
        prots = []
        element_name = 'protein'
        for prot_group_child in prot_group.getchildren():
            if prot_group_child.tag[-len(element_name):] == element_name:
                prots.append(Protein(prot_group_child))

        return prots

    def get_prots(self):
        return self.prots


class Protein():
    def __init__(self, prot):
        self.prot_name = prot.attrib["protein_name"]
        self.desc = prot.annotation.attrib["protein_description"]
        self.prob = float(prot.attrib["probability"])
        self.conf = float(prot.attrib["confidence"])
        self.nr_peptides = int(prot.attrib["total_number_distinct_peptides"])
        self.group_sibling_id = prot.attrib["group_sibling_id"]

        # No peptides -> no coverage, duh!
        try:
            self.cov = prot.attrib["percent_coverage"]
        except KeyError:
            self.cov = 0

        # There's not always a spectrum id?
        try:
            self.spec_ids = float(prot.attrib["pct_spectrum_ids"])
        except KeyError:
            self.spec_ids = None

        self.peptides = self.set_peptides(prot)

        element_name = "indistinguishable_protein"
        if element_name in [prot_child.tag[-len(element_name):] for prot_child in prot.getchildren()]:
            self.indistinguishable_proteins = self.set_indistinguishable_proteins(prot)

    @staticmethod
    def child_is_peptide(peptides, prot_child):
        # We are only interested in peptides with an initial_probability
        # of at least > 0
        element_name = 'peptide'
        if prot_child.tag[-len(element_name):] == element_name and float(prot_child.attrib["initial_probability"]) > 0:

            # To avoid duplicate peptides detected on different charges
            if prot_child.attrib["peptide_sequence"] not in [peptide.get_seq() for peptide in peptides]:
                return True

    def set_peptides(self, prot):
        peptides = []
        for prot_child in prot.getchildren():
            if self.child_is_peptide(peptides, prot_child):
                peptides.append(Peptide(prot_child))

        return peptides

    @staticmethod
    def set_indistinguishable_proteins(protein):
        indistinguishable_proteins = []
        element_name = 'indistinguishable_protein'
        for prot_child in protein.getchildren():
            if prot_child.tag[-len(element_name):] == element_name:
                indistinguishable_proteins.append(prot_child.attrib)

        return indistinguishable_proteins

    def get_prot_name(self):
        return self.prot_name

    def get_descr(self):
        return self.desc

    def get_seq(self, protein_db):
        return str(protein_db[self.prot_name].seq)

    def get_prob(self):
        return self.prob

    def get_peptides(self):
        return self.peptides


class Peptide():
    """
 
    :param peptide:
    """

    def __init__(self, peptide):
        self.seq = peptide.attrib["peptide_sequence"]
        self.charge = peptide.attrib["charge"]
        self.initial_prob = peptide.attrib["initial_probability"]
        self.nsp_prob = peptide.attrib["nsp_adjusted_probability"]
        self.fpkm_prob = peptide.attrib["nsp_adjusted_probability"]

        if peptide.attrib["is_contributing_evidence"] == "Y":
            self.is_contributing_evidence = True
        else:
            self.is_contributing_evidence = False

        element_name = "peptide_parent_protein"
        if element_name in [pep_child.tag[-len(element_name):] for pep_child in peptide.getchildren()]:
            self.peptide_parent_proteins = self.set_peptide_parent_proteins(
                peptide)

    @staticmethod
    def set_peptide_parent_proteins(peptide):
        parent_proteins = []
        element_name = 'peptide_parent_protein'
        for peptide_child in peptide.getchildren():
            if peptide_child.tag[-len(element_name):] == element_name:
                parent_proteins.append(peptide_child.attrib)

        return parent_proteins

    def get_seq(self):
        return self.seq