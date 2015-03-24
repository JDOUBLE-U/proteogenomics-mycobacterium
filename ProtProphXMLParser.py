#!/usr/bin/python
# coding=utf-8
import os

from lxml import objectify

"""
Jeroen Merks

...
"""


class ProtXMLParser(object):
    """

    :param parse_results_folder:
    :param prot_prophet_xml:
    """

    def __init__(self, parse_results_folder, prot_prophet_xml):
        if not os.path.exists(parse_results_folder):
            os.mkdir(parse_results_folder)

        self.PARSE_RESULTS_FOLDER = parse_results_folder
        self.PROT_PROPHET_XML = prot_prophet_xml
        self.prots = self.get_prots()
        self.statistics_dict = self.get_statistics_dict()

    def get_prots(self):
        """
        Extracts groups of proteins by looping through the "protein_group" XML-elements of the prot.xml file.

        :return: list
        """
        prots = []
        prot_prophet_xml_object = open(self.PROT_PROPHET_XML, encoding="utf-8")
        prot_prophet_xml_object.readline()  # skip XML encoding declaration
        for prot in objectify.fromstring(prot_prophet_xml_object.read()).protein_group:
            prots.append(ProtProphHit(prot))

        prot_prophet_xml_object.close()

        return prots

    def get_statistics_dict(self):
        """

        :return:
        """
        statistics_dict = {}
        prot_prophet_xml_object = open(self.PROT_PROPHET_XML, encoding="utf-8")
        prot_prophet_xml_object.readline()  # skip XML encoding declaration
        for data_filter in objectify.fromstring(
                prot_prophet_xml_object.read()). \
                protein_summary_header.program_details.proteinprophet_details.getchildren():

            # We are only interested in protein_summary_data
            if data_filter.tag[-27:] == 'protein_summary_data_filter':
                min_probability = float(data_filter.attrib["min_probability"])
                error_rate = float(data_filter.attrib["false_positive_error_rate"])

                statistics_dict[error_rate] = min_probability

        prot_prophet_xml_object.close()

        return statistics_dict


class ProtProphHit(object):
    """

    :param prot:
    """

    def __init__(self, prot):
        self.prob = float(prot.attrib["probability"])
        self.conf = float(prot.protein.attrib["confidence"])
        self.group_nr = int(prot.attrib["group_number"])
        self.id = prot.protein.attrib["protein_name"]
        self.description = prot.protein.annotation.attrib["protein_description"]
        self.nr_peptides = int(prot.protein.attrib["total_number_distinct_peptides"])

        # No peptides -> no coverage, duh!
        try:
            self.cov = prot.protein.attrib["percent_coverage"]
        except KeyError:
            self.cov = 0

        # There's not always a spectrum id?
        try:
            self.spec_ids = float(prot.protein.attrib["pct_spectrum_ids"])
        except KeyError:
            self.spec_ids = None

        self.peptides = self.set_peptides(prot)

    def get_goup_nr(self):
        return self.group_nr

    @staticmethod
    def set_peptides(prot_group):
        peptides = []
        for pep in prot_group.protein.getchildren():
            # We are only interested in peptides with an initial_probability higher then 0
            if pep.tag[-7:] == 'peptide' and float(pep.attrib["initial_probability"]) > 0:

                # To avoid duplicate peptides detected on different charges
                if pep.attrib["peptide_sequence"] not in [peptide.get_seq() for peptide in peptides]:
                    peptides.append(Peptide(pep.attrib))

        return peptides

    def get_descr(self):
        return self.description

    def get_id(self):
        return self.id

    def get_seq(self, protein_db):
        return str(protein_db[self.id].seq)

    def get_prob(self):
        return self.prob

    def get_peptides(self):
        return self.peptides


class Peptide(object):
    """

    :param peptide_attributes:
    """

    def __init__(self, peptide_attributes):
        self.seq = peptide_attributes["peptide_sequence"]
        self.charge = peptide_attributes["charge"]
        self.initial_prob = peptide_attributes["initial_probability"]
        self.nsp_prob = peptide_attributes["nsp_adjusted_probability"]
        self.fpkm_prob = peptide_attributes["nsp_adjusted_probability"]

        self.is_contrib_evidence = peptide_attributes["is_contributing_evidence"]
        if self.is_contrib_evidence == "Y":
            self.is_contributing_evidence = True
        else:
            self.is_contributing_evidence = False

    def get_seq(self):
        return self.seq

    def get_charge(self):
        return self.charge
