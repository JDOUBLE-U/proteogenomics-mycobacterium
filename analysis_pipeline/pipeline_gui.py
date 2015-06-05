# import modules
import sys
sys.path.append('..\\')
import os

import wx
from analysis_pipeline import preprocess_db
from Comet_20151 import run_comet as comet
from msconvert import run_msconvert as msconvert
from EMBOSS_sixpack_65 import run_sixpack as sixpack
from xinteract import run_xinteract as xinteact
from analysis_pipeline import find_metap_activity


def delete_previous_results():
    out_folders = [out_folder for out_folder in os.listdir('..\\analysis_pipeline')
                   if out_folder.endswith('_out')]
    for out_folder in out_folders:
        files = os.listdir(out_folder)
        for file_name in files:
            os.remove(out_folder + '\\' + file_name)


def get_mzxmls(mzxml_folder):
    mzxml_file_names = []
    for fileName in os.listdir(mzxml_folder):
        if fileName.endswith('.mzML'):
            mzxml_file_names.append(mzxml_folder + fileName)

    return mzxml_file_names


def create_out_folder(out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)


def run_pipeline(clear_prev_results, spectra_folder_location, output_folder, proteome_path, on_sixframe, genome_db_path,
                 codon_table, motif_range_end, min_prob, min_pep_length, cleavage_loc, on_digest, nr_threads):
    # Windows executables
    on_os = sys.platform
    if on_os == 'win32':
        sixpack_executable = 'cd ..\\EMBOSS_sixpack_65& windows_sixpack.exe'
        comet_executable = 'cd ..\\Comet_20151& comet.2015011.win64.exe'
        xinteract_executable = 'xinteract\\xinteract.exe'
        msconvert_executable = 'cd ..\\msconvert& msconvert.exe'
    # Linux and Mac executables
    else:
        sixpack_executable = '.\\EMBOSS_sixpack_65\\linux_sixpack'
        comet_executable = '..\\Comet_20151\\comet.2015011.linux.exe'
        xinteract_executable = 'xxx'
        msconvert_executable = 'xxx'

    ## Create all output folders ##
    create_out_folder("comet_out\\")
    create_out_folder("sixpack_out\\")
    create_out_folder("xinteract_out\\")
    create_out_folder("protxml_parsing_out\\")
    create_out_folder('weblogo_out\\')
    create_out_folder('preprocessed_db_out\\')

    ## Clear previous results ##
    if clear_prev_results:
        delete_previous_results()

    ## Convert any raw files to mzxml files ##
    msconvert.convert_all_raw_files(msconvert_executable, spectra_folder_location)

    ## SixPack ##
    if on_sixframe:
        proteome_path = sixpack.run_sixpack(genome_db_path, sixpack_executable, on_os, min_pep_length + 1,
                                            codon_table)

    ## Process protein db ##
    if on_digest:
        processed_prot_db = preprocess_db.cleave_m_and_digest(proteome_path, min_pep_length)
    else:
        processed_prot_db = preprocess_db.cleave_m_only(proteome_path, min_pep_length)

    ## Comet ##
    comet_pep_xmls = []
    mzxmls = get_mzxmls(spectra_folder_location)
    for mzxml in mzxmls:
        comet_pep_xmls.append(comet.run_comet(comet_executable, processed_prot_db, mzxml, nr_threads))

    ## xinteract ##
    xinteact.run_xinteract(output_folder, xinteract_executable, comet_pep_xmls, min_pep_length, nr_threads)

    ## Find MetAp activity ##
    find_metap_activity.run_metap_pipeline(output_folder,
                                           "xinteract_out\\" + output_folder + '.interact.prot.xml',
                                           proteome_path, min_prob, cleavage_loc,
                                           motif_range_end)


class ProgPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id, style=wx.BORDER_SIMPLE)
        self.font = wx.Font(13, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        self.SetFont(self.font)
        self.mainbox = wx.BoxSizer()
        self.progtext = wx.TextCtrl(self, -1, "", style=wx.TE_READONLY | wx.TE_MULTILINE)
        self.mainbox.Add(self.progtext, 1, wx.EXPAND)
        sys.stdout = self.progtext
        self.SetSizer(self.mainbox)


class InputPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self.font = wx.Font(13, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        self.SetFont(self.font)
        #### Create UI elements
        ##create sizers
        self.mainbox = wx.BoxSizer()
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.mainbox.Add(self.vbox, 1, wx.EXPAND | wx.ALL, 7)
        self.SetSizer(self.mainbox)
        self.mainbox.AddSpacer(3)

        ## Select spectra folder
        self.vbox.Add(wx.StaticText(self, -1, "MS/MS Spectra folder location"))
        self.specfol = wx.DirPickerCtrl(self, -1, "", "Select MS/MS spectra folder")
        self.vbox.Add(self.specfol, 0, wx.EXPAND)

        ## Name of output files
        self.vbox.Add(wx.StaticText(self, -1, "Prefix name of output files"))
        self.namectrl = wx.TextCtrl(self, -1, "")
        self.vbox.Add(self.namectrl, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)

        ## Six frame or not
        self.vbox.Add(wx.StaticText(self, -1, "Analyse on the six-frame translation of the genome?"))
        self.radio_no = wx.RadioButton(self, -1, "No", style=wx.RB_GROUP)
        self.radio_yes = wx.RadioButton(self, -1, "Yes")
        self.sixchoice = False
        self.sixbox = wx.BoxSizer()
        self.sixbox.Add(self.radio_no, 0, wx.EXPAND)
        self.sixbox.Add(wx.StaticText(self, -1, "    "), 0, wx.EXPAND)
        self.sixbox.Add(self.radio_yes, 0, wx.EXPAND)
        self.vbox.Add(self.sixbox, 0, wx.EXPAND)
        self.radio_no.Bind(wx.EVT_RADIOBUTTON, self.onSixRadio)
        self.radio_yes.Bind(wx.EVT_RADIOBUTTON, self.onSixRadio)
        self.vbox.AddSpacer(5)

        ## Proteome location
        self.prottext = wx.StaticText(self, -1, "Proteome file location")
        self.vbox.Add(self.prottext, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.protfol = wx.FilePickerCtrl(self, -1, "", "Select proteome file")
        self.vbox.Add(self.protfol, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

        ## Genome location if using six frame
        self.genotext = wx.StaticText(self, -1, "Genome file location")
        self.vbox.Add(self.genotext, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.genotext)
        self.genofol = wx.FilePickerCtrl(self, -1, "", "Select Genome file")
        self.vbox.Add(self.genofol, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.genofol)
        self.Layout()

        ## startcodons
        self.codontext = wx.StaticText(self, -1, "Genetics code used for the translation")
        self.vbox.Add(self.codontext, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.codontext)
        self.choices = ["Standard", "Standard (with alternative initiation codons)",
                        "Vertebrate Mitochondrial", "Yeast Mitochondrial",
                        "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma",
                        "Invertebrate Mitochondrial", "Ciliate Macronuclear and Dasycladacean",
                        "Echinoderm Mitochondrial", "Euplotid Nuclear", "Bacterial",
                        "Alternative Yeast Nuclear", "Ascidian Mitochondrial", "Flatworm Mitochondrial",
                        "Blepharisma Macronuclear", "Chlorophycean Mitochondrial", "Trematode Mitochondrial",
                        "Scenedesmus obliquus", "Thraustochytrium Mitochondrial"]
        self.codonmenu = wx.ComboBox(self, -1, "Standard", choices=self.choices,
                                     style=wx.CB_READONLY)
        self.codonchoice = 0
        self.constructCodontable(0)
        self.vbox.Add(self.codonmenu, 0, wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.codonmenu)
        self.codonmenu.Bind(wx.EVT_COMBOBOX, self.onCombo)

        ## Weblogo motif range
        self.vbox.AddSpacer(5)

        ## Weblogo motif range
        self.vbox.Add(wx.StaticText(self, -1, "WebLogo motif range"))
        self.wlspin = wx.SpinCtrl(self, -1, "5", min=1, max=7)
        self.vbox.Add(self.wlspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Maximum error rate
        self.vbox.Add(wx.StaticText(self, -1, "Maximum peptide error rate"))
        self.radio_1 = wx.RadioButton(self, -1, "1%")
        self.radio_5 = wx.RadioButton(self, -1, "5%", style=wx.RB_GROUP)
        self.maxerrbox = wx.BoxSizer()
        self.maxerrbox.Add(self.radio_1, 0, wx.EXPAND)
        self.maxerrbox.Add(wx.StaticText(self, -1, "   "), 0)
        self.maxerrbox.Add(self.radio_5, 0, wx.EXPAND)
        self.radio_1.Bind(wx.EVT_RADIOBUTTON, self.onErrorRadio)
        self.radio_5.Bind(wx.EVT_RADIOBUTTON, self.onErrorRadio)
        self.vbox.Add(self.maxerrbox, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Minimum length of peptides
        self.vbox.Add(wx.StaticText(self, -1, "Minimum peptide length"))
        self.peplenspin = wx.SpinCtrl(self, -1, "5", min=5, max=45)
        self.vbox.Add(self.peplenspin, 0, wx.SHAPED)
        self.peplenspin.Bind(wx.EVT_SPINCTRL, self.onPepLenSpin)
        self.vbox.AddSpacer(5)

        ## cleavage location
        self.vbox.Add(wx.StaticText(self, -1, "Specific starting location of peptide on protein"))
        self.cleavespin = wx.SpinCtrl(self, -1, "1", min=1)
        self.vbox.Add(self.cleavespin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Digest before starting
        self.vbox.Add(wx.StaticText(self, -1, "Digest proteome database with Trypsin before starting?"))
        self.digradio_y = wx.RadioButton(self, -1, "Yes", style=wx.RB_GROUP)
        self.digradio_n = wx.RadioButton(self, -1, "No")
        self.digestchoice = False
        self.digestbox = wx.BoxSizer()
        self.digestbox.Add(self.digradio_n, 0, wx.EXPAND)
        self.digestbox.Add(wx.StaticText(self, -1, "    "), 0)
        self.digestbox.Add(self.digradio_y, 0, wx.EXPAND)
        self.vbox.Add(self.digestbox, 0, wx.EXPAND | wx.ALL, 1)
        self.digradio_n.Bind(wx.EVT_RADIOBUTTON, self.onDigestRadio)
        self.digradio_y.Bind(wx.EVT_RADIOBUTTON, self.onDigestRadio)
        self.vbox.AddSpacer(5)

        ## Number of threads
        self.vbox.Add(wx.StaticText(self, -1, "Number of CPU-threads to use"))
        self.threadspin = wx.SpinCtrl(self, -1, "8", min=1, max=100)
        self.vbox.Add(self.threadspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Start button
        self.start = wx.Button(self, -1, "Start")
        self.vbox.Add(self.start, 0, wx.SHAPED | wx.ALIGN_CENTER_HORIZONTAL)
        self.start.Bind(wx.EVT_BUTTON, self.onStart)

    def onPepLenSpin(self, event):
        self.wlspin.SetRange(1, 35)
        self.cleavespin.SetRange(1, 10000)

    def onSixRadio(self, event):
        if event.GetEventObject().GetLabel() == "Yes":
            self.sixchoice = True
            self.vbox.Show(self.genofol)
            self.vbox.Show(self.genotext)
            self.vbox.Show(self.codontext)
            self.vbox.Show(self.codonmenu)
            self.vbox.Hide(self.prottext)
            self.vbox.Hide(self.protfol)
        else:
            self.sixchoice = False
            self.vbox.Hide(self.genofol)
            self.vbox.Hide(self.genotext)
            self.vbox.Hide(self.codontext)
            self.vbox.Hide(self.codonmenu)
            self.vbox.Show(self.prottext)
            self.vbox.Show(self.protfol)
        self.Layout()

    def onCombo(self, event):
        self.codonchoice = self.indexToNCBINumber(self.codonmenu.GetSelection())

    def indexToNCBINumber(self, index_in):
        if index_in > 6:
            ncbinums = {7: 9, 8: 10, 9: 11, 10: 12, 11: 13, 12: 14, 13: 15, 14: 16, 15: 21, 16: 22, 17: 23}
            return ncbinums[index_in]
        else:
            return index_in

    def onErrorRadio(self, event):
        if event.GetEventObject().GetLabel() == "1%":
            self.max_error = 0.01
        else:
            self.max_error = 0.05

    def onDigestRadio(self, event):
        if event.GetEventObject().GetLabel() == "Yes":
            self.digestchoice = True
        else:
            self.digestchoice = False

    def GetSpectraFolder(self):
        return self.specfol.GetPath() + "\\"

    def GetOutputFolder(self):
        return self.outfol.GetPath() + "\\"

    def GetOutputNames(self):
        return self.namectrl.GetValue()

    def GetProteomeFile(self):
        return self.protfol.GetPath()

    def GetOnSixframe(self):
        return self.sixchoice

    def GetGenomeFile(self):
        if self.sixchoice == 1:
            return self.genofol.GetPath()

    def GetCodonTable(self):
        if self.sixchoice == 1:
            return self.codonmenu.GetSelection()

    def GetMotifLen(self):
        return self.wlspin.GetValue()

    def GetMaxError(self):
        return self.max_error

    def GetMinPepLen(self):
        return self.peplenspin.GetValue()

    def GetCleavageLoc(self):
        return self.cleavespin.GetValue()

    def GetOnDigest(self):
        return self.digestchoice

    def GetNrThreads(self):
        return self.threadspin.GetValue()

    def checkInputs(self):
        if self.namectrl.GetValue().isalnum():
            return True
        else:
            self.errorDialog("Output filename may only contain letters and numbers and can not be empty")
            return False

    def errorDialog(self, message):
        dial = wx.MessageBox(message, 'Error', wx.OK | wx.ICON_ERROR)

    def constructCodontable(self, index):
        self.codondict = {}
        linenum = 0
        with open("codontables.txt", "r") as f:
            for line in f:
                if linenum == index:
                    for pair in line.split(","):
                        pair = pair.split(":")
                        self.codondict[pair[0]] = pair[1]
                linenum += 1

    def onStart(self, event):
        if self.checkInputs():
            run_pipeline(True, self.GetSpectraFolder(), self.GetOutputNames(),
                         self.GetProteomeFile(), self.GetOnSixframe(),
                         self.GetGenomeFile(), self.GetCodonTable(), self.GetMotifLen(), float(0.95),
                         self.GetMinPepLen(), self.GetCleavageLoc(), self.GetOnDigest(), self.GetNrThreads())
            sys.exit()


class Mainframe(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(1000, 680))
        panel = wx.Panel(self, -1)
        self.inpanel = InputPanel(panel, -1)
        self.progpanel = ProgPanel(panel, -1)
        mainbox = wx.BoxSizer()
        mainbox.Add(self.inpanel, 3, wx.EXPAND)
        mainbox.Add(self.progpanel, 2, wx.EXPAND | wx.LEFT, 1)
        panel.SetSizer(mainbox)
        self.Centre()
        self.Show(True)


if __name__ == "__main__":
    # delete_previous_results()
    app = wx.App(False)
    program = Mainframe(None, -1, "title")
    app.MainLoop()

# start_time = timeit.default_timer()
# stop_time = timeit.default_timer()
# print("The Pipeline took %i seconds to run\n" % stop_time)
