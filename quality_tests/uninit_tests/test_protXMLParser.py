from unittest import TestCase

from Bio import SeqIO

from Comet_analysis_pipeline.ProtProphXMLParser import ProtProphXMLParser
from Comet_analysis_pipeline.ProtProphXMLParser import ProteinGroup

# TODO make file references relative

class TestProtXMLParser(TestCase):
    # Check if all protein groups get properly fetched and return correct group numbers and probabilities
    def test_get_prot_groups(self):
        protxmlparser = ProtProphXMLParser
        protxmlparser.prot_prophet_xml = "./GitHub_test_files/prophet_test_xml.xml"

        prot_groups = protxmlparser.get_prot_groups(protxmlparser)

        for prot_group in prot_groups:
            self.assertIsInstance(prot_group, ProteinGroup)

            self.assertIn(prot_group.group_nr, [1335, 1336])
            self.assertIn(prot_group.prob, [1.0000, 0.8])

        # Check if all protein groups are fetched
        self.assertEqual(len(prot_groups), 2)

    # Check if protein_summary_data_filter gets fetched and the proper false_positive_error rates are getting tied to the
    # proper min_probability rates
    def test_get_statistics_dict(self):
        protxmlparser = ProtProphXMLParser
        protxmlparser.prot_prophet_xml = "./GitHub_test_files/prophet_test_xml.xml"

        statistics_dict = protxmlparser.get_error_prob_pairs(protxmlparser)

        self.assertEqual(statistics_dict[0.002], 0.99)
        self.assertEqual(statistics_dict[0.075], 0.3)


class TestProteinGroup(TestCase):
    def test_set_and_get_prot_groups(self):

        protxmlparser = ProtProphXMLParser
        protxmlparser.prot_prophet_xml = "./GitHub_test_files/prophet_test_xml.xml"

        # Check if nr of proteins in group 0 is correct
        self.assertEqual(len(protxmlparser.get_prot_groups(protxmlparser)[0].get_prots()), 2)
        # Check if all proteins are properly fetched
        for prot in protxmlparser.get_prot_groups(protxmlparser)[0].get_prots():
            self.assertIn(prot.get_prot_name(), ["sp|P9WK95|LEUD_MYCTU", "tr|O53945|O53945_MYCTU"])

        # Check if nr of proteins in group 1 is correct
        self.assertEqual(len(protxmlparser.get_prot_groups(protxmlparser)[1].get_prots()), 3)
        # Check if all proteins are properly fetched
        for prot in protxmlparser.get_prot_groups(protxmlparser)[1].get_prots():
            self.assertIn(prot.get_prot_name(),
                ["sp|P9WNJ5|ESXL_MYCTU", "tr|O53945|O53945_MYCTU", "sp|P9WNJ3|ESXN_MYCTU",
                 "sp|P9WNI7|ES6L6_MYCTU"])


class TestProtein(TestCase):
    def test_test_get_peptide(self):
        protxmlparser = ProtProphXMLParser
        test_prot_db = SeqIO.index("./GitHub_test_files/prot_test_db", format='fasta')
        protxmlparser.prot_prophet_xml = "./GitHub_test_files/prophet_test_xml.xml"
        prot_groups = protxmlparser.get_prot_groups(protxmlparser)

        group0prot0 = prot_groups[0].get_prots()[0]
        self.assertEqual(group0prot0.get_seq(test_prot_db),
                         "MEAFHTHSGIGVPLRRSNVDTDQIIPAVFLKRVTRTGFEDGLFAGWRSDPAFVLNLSPFDRGSVLVAGPDFGTGSSREHAVWALMDYGFRVVISSRFGDIFRGNAGKAGLLAAEVAQDDVELLWKLIEQSPGLEITANLQDRIITAATVVLPFKIDDHSAWRLLEGLDDIALTLRKLDEIEAFEGACAYWKPRTLPAP")

        group0prot1 = prot_groups[0].get_prots()[1]
        self.assertEqual(group0prot1.get_seq(test_prot_db),
                         "MQRFGTGSSRSWCGRAGTATIAAVLLASGALTGLPPAYAISPPTIDPGALPPDGPPGPLAPMKQNAYCTEVGVLPGTDFQLQPKYMEMLNLNEAWQFGRGDGVKVAVIDTGVTPHPRLPRLIPGGDYVMAGGDGLSDCDAHGTLVASMIAAVPANGAVPLPSVPRRPVTIPTTETPPPPQTVTLSPVPPQTVTVIPAPPPEEGVPPGAPVPGPEPPPAPGPQPPAVDRGGGTVTVPSYSGGRKIAPIDNPRNPHPSAPSPALGPPPDAFSGIAPGVEIISIRQSSQAFGLKDPYTGDEDPQTAQKIDNVETMARAIVHAANMGASVINISDVMCMSARNVIDQRALGAAVHYAAVDKDAVIVAAAGDGSKKDCKQNPIFDPLQPDDPRAWNAVTTVVTPSWFHDYVLTVGAVDANGQPLSKMSIAGPWVSISAPGTDVVGLSPRDDGLINAIDGPDNSLLVPAGTSFSAAIVSGVAALVRAKFPELSAYQIINRLIHTARPPARGVDNQVGYGVVDPVAALTWDVPKGPAEPPKQLSAPLVVPQPPAPRDMVPIWVAAGGLAGALLIGGAVFGTATLMRRSRKQQ")

        group1prot0 = prot_groups[1].get_prots()[0]
        self.assertEqual(group1prot0.get_seq(test_prot_db),
                         "MTINYQFGDVDAHGAMIRAQAGLLEAEHQAIIRDVLTASDFWGGAGSAACQGFITQLGRNFQVIYEQANAHGQKVQAAGNNMAQTDSAVGSSWA")

        group1prot1 = prot_groups[1].get_prots()[1]
        self.assertEqual(group1prot1.get_seq(test_prot_db),
                         "MTINYQFGDVDAHGAMIRAQAASLEAEHQAIVRDVLAAGDFWGGAGSVACQEFITQLGRNFQVIYEQANAHGQKVQAAGNNMAQTDSAVGSSWA")


        def test_set_peptides(self):
            # Test if all peptides of a given protein are fetched
            self.assertEqual(len(prot_groups[0].get_prots()[0].get_peptides()), 2)
            self.assertEqual(len(prot_groups[0].get_prots()[1].get_peptides()), 2)
            self.assertEqual(len(prot_groups[2].get_prots()[0].get_peptides()), 6)
            self.assertEqual(len(prot_groups[2].get_prots()[0].get_peptides()), 5)
