# !/usr/bin/python
# coding=utf-8
import sys

sys.path.append('..\\')
import os
import subprocess

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


def run_pipeline(clear_prev_results, spectra_folder_location, output_folder, proteome_path, enzyme_choice,
                 allowed_missed_cleavage,
                 motif_len, min_prob, min_pep_length, mass_tolerance, pep_mass_unit, isotope_error, nr_threads):
    # on_os = sys.platform

    # sixpack_executable = 'cd ..\\EMBOSS_sixpack_65& windows_sixpack.exe'
    comet_executable = 'cd ..\\Comet_20151& comet.2015011.win64.exe'
    xinteract_executable = 'xinteract\\xinteract.exe'
    msconvert_executable = 'msconvert.exe'


    ## Create all output folders ##
    create_out_folder("comet_out\\")
    # create_out_folder("sixpack_out\\")
    create_out_folder("xinteract_out\\")
    create_out_folder("protxml_parsing_out\\")
    create_out_folder('weblogo_out\\')
    create_out_folder('preprocessed_db_out\\')

    ## Clear previous results ##
    if clear_prev_results:
        delete_previous_results()

    ## Convert any raw files to mzxml files ##
    msconvert.main(msconvert_executable, spectra_folder_location, nr_threads)

    # ## Process protein db ##
    # if on_digest:
    #     processed_prot_db = preprocess_db.cleave_m_and_digest(proteome_path, min_pep_length)
    # else:
    processed_prot_db = preprocess_db.cleave_m_only(proteome_path, min_pep_length)

    ## Comet ##
    comet_pep_xmls = []
    mzxmls = get_mzxmls(spectra_folder_location)
    for mzxml in mzxmls:
        comet_pep_xmls.append(
            comet.run_comet(comet_executable, processed_prot_db, mzxml, mass_tolerance, pep_mass_unit, isotope_error,
                            enzyme_choice, allowed_missed_cleavage,
                            nr_threads))

    ## xinteract ##
    xinteact.run_xinteract(output_folder, xinteract_executable, comet_pep_xmls, min_pep_length, nr_threads)

    ## Find proteins with peptides on the given location ##
    find_metap_activity.main(output_folder, "xinteract_out\\" + output_folder + '.interact.prot.xml',
                             processed_prot_db, min_prob,
                             motif_len)


class InputPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self.font = wx.Font(13, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        self.SetFont(self.font)
        #### Create UI elements
        ##create sizers
        self.mainbox = wx.BoxSizer()
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.mainbox.Add(self.vbox, 1, wx.EXPAND | wx.ALL, 5)
        self.SetSizer(self.mainbox)
        self.mainbox.AddSpacer(3)

        ## Select spectra folder
        self.vbox.Add(wx.StaticText(self, -1, "MS/MS spectra folder location"))
        self.specfol = wx.DirPickerCtrl(self, -1, "", "Select MS/MS spectra folder")
        self.vbox.Add(self.specfol, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)

        ## Name of output files
        self.vbox.Add(wx.StaticText(self, -1, "Prefix name of output files"))
        self.namectrl = wx.TextCtrl(self, -1, "")
        self.vbox.Add(self.namectrl, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)

        ## Proteome location
        self.prottext = wx.StaticText(self, -1, "Proteome file location")
        self.vbox.Add(self.prottext, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.protfol = wx.FilePickerCtrl(self, -1, "", "Select proteome file")
        self.vbox.Add(self.protfol, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

        ## digestion enzyme
        self.enzymetext = wx.StaticText(self, -1, "Digestion enzyme")
        self.vbox.Add(self.enzymetext, 0, wx.EXPAND)
        self.choices = ["No enzyme", "Trypsin", "Trypsin/P", "Lys C", "Lys N", "Arg C", "Asp N",
                        "CNBr", "Glu C", "PepsinA", "Chymotrypsin"]
        self.enzymemenu = wx.ComboBox(self, -1, "Trypsin", choices=self.choices,
                                      style=wx.CB_READONLY)
        self.enzymechoice = 0
        self.vbox.Add(self.enzymemenu, 0)
        self.enzymemenu.Bind(wx.EVT_COMBOBOX, self.onCombo)
        self.vbox.AddSpacer(5)

        ## max number of overcleavages
        self.vbox.Add(wx.StaticText(self, -1, "Maximum number of allowed overcleavages"))
        self.maxoverspin = wx.SpinCtrl(self, -1, "2", min=0, max=5)
        self.vbox.Add(self.maxoverspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Weblogo motif range
        self.vbox.Add(wx.StaticText(self, -1, "WebLogo motif range"))
        self.wlspin = wx.SpinCtrl(self, -1, "5", min=2, max=10)
        self.vbox.Add(self.wlspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Minnimum probability protein
        self.vbox.Add(wx.StaticText(self, -1, "Minnimum initial probability on protein level"))
        self.errpercspin = wx.SpinCtrlDouble(self, -1, "0.95", min=float(0.05), max=float(1.00))
        self.errpercspin.SetIncrement(0.01)
        self.vbox.Add(self.errpercspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Minimum length of peptides
        self.vbox.Add(wx.StaticText(self, -1, "Minimum peptide length"))
        self.peplenspin = wx.SpinCtrl(self, -1, "5", min=5, max=45)
        self.vbox.Add(self.peplenspin, 0, wx.SHAPED)
        self.peplenspin.Bind(wx.EVT_SPINCTRL, self.onPepLenSpin)
        self.vbox.AddSpacer(5)

        ## peptide mass tolerance
        self.vbox.Add(wx.StaticText(self, -1, "Peptide mass tolerance"))
        self.peptolspin = wx.SpinCtrlDouble(self, -1, "3.00", min=1, max=10)
        self.peptolspin.SetIncrement(0.01)
        self.vbox.Add(self.peptolspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## peptide mass units
        self.vbox.Add(wx.StaticText(self, -1, "Peptide mass units (0=amu, 1=mmu, 2=ppm)"))
        self.pepunitspin = wx.SpinCtrl(self, -1, "0", min=0, max=2)
        self.vbox.Add(self.pepunitspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Isotope error
        self.vbox.Add(wx.StaticText(self, -1, "Isotope error "))
        self.choices = ["Off", "On, -1/0/1/2/3 (standard C13 error)", "On, -8/-4/0/4/8 (for +4/+8 labeling)"]
        self.isomenu = wx.ComboBox(self, -1, "Off", choices=self.choices, style=wx.CB_READONLY)
        self.isochoice = 0
        self.vbox.Add(self.isomenu, 0)
        self.isomenu.Bind(wx.EVT_COMBOBOX, self.onIsoCombo)
        self.vbox.AddSpacer(5)

        ## Number of threads
        self.vbox.Add(wx.StaticText(self, -1, "Number of CPU-threads to use"))
        self.threadspin = wx.SpinCtrl(self, -1, str(os.cpu_count()), min=1, max=100)
        self.vbox.Add(self.threadspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## delete previous results
        self.vbox.Add(wx.StaticText(self, -1, "Delete previous results?"))
        self.prevdelradio_n = wx.RadioButton(self, -1, "No", style=wx.RB_GROUP)
        self.prevdelradio_y = wx.RadioButton(self, -1, "Yes")
        self.prevdelchoice = False
        self.prevdelbox = wx.BoxSizer()
        self.prevdelbox.Add(self.prevdelradio_n, 0, wx.EXPAND)
        self.prevdelbox.Add(wx.StaticText(self, -1, "    "), 0)
        self.prevdelbox.Add(self.prevdelradio_y, 0, wx.EXPAND)
        self.vbox.Add(self.prevdelbox, 0, wx.EXPAND | wx.ALL, 1)
        self.prevdelradio_n.Bind(wx.EVT_RADIOBUTTON, self.onPrevdelRadio)
        self.prevdelradio_y.Bind(wx.EVT_RADIOBUTTON, self.onPrevdelRadio)
        self.vbox.AddSpacer(5)

        ## Start button
        self.start = wx.Button(self, -1, "Start")
        self.vbox.Add(self.start, 0, wx.SHAPED | wx.ALIGN_CENTER_HORIZONTAL)
        self.start.Bind(wx.EVT_BUTTON, self.onStart)
        self.vbox.AddSpacer(5)

        ## open results folder
        self.openresfol = wx.Button(self, -1, "View results folder")
        self.vbox.Add(self.openresfol, 0, wx.SHAPED | wx.ALIGN_CENTER_HORIZONTAL)
        self.openresfol.Bind(wx.EVT_BUTTON, self.onOpenResFol)

    def onPepLenSpin(self, event):
        self.wlspin.SetRange(1, 35)
        self.cleavespin.SetRange(1, 5)

    def onCombo(self, event):
        self.enzymechoice = self.enzymemenu.GetSelection()

    def onIsoCombo(self, event):
        self.isochoice = self.isomenu.GetSelection()

    def onPrevdelRadio(self, event):
        if event.GetEventObject().GetLabel() == "Yes":
            self.prevdelchoice = True
        else:
            self.prevdelchoice = False

    def GetSpectraFolder(self):
        return self.specfol.GetPath() + "\\"

    def GetOutputNames(self):
        return self.namectrl.GetValue()

    def GetProteomeFile(self):
        return self.protfol.GetPath()

    def GetDigestEnzyme(self):
        return self.enzymechoice

    def GetAllowedMissedCleavages(self):
        return self.maxoverspin.GetValue()

    def GetMotifLen(self):
        return self.wlspin.GetValue()

    def GetMinPercentage(self):
        return self.errpercspin.GetValue()

    def GetMinPepLen(self):
        return self.peplenspin.GetValue()

    def GetNrThreads(self):
        return self.threadspin.GetValue()

    def GetPepMassTol(self):
        return self.peptolspin.GetValue()

    def GetPepMassUnit(self):
        return self.pepunitspin.GetValue()

    def GetIsotopeError(self):
        return self.isochoice

    def GetDeletePrevRes(self):
        return self.prevdelchoice

    def checkInputs(self):
        outname = self.namectrl.GetValue()
        if len(set(['?', '<', '>', ':', '*', '|', '"', '^']) - set(outname)) >= 8:
            return True
        else:
            self.errorDialog("Filename contains illegal characters or is empty")
            return False

    def errorDialog(self, message):
        dial = wx.MessageBox(message, 'Error', wx.OK | wx.ICON_ERROR)

    def onOpenResFol(self, event):
        subprocess.Popen('explorer /open,"' + str(os.getcwd()) + '\\protxml_parsing_out"')

    def onStart(self, event):
        if self.checkInputs():
            run_pipeline(self.GetDeletePrevRes(), self.GetSpectraFolder(), self.GetOutputNames(),
                         self.GetProteomeFile(), self.GetDigestEnzyme(), self.GetAllowedMissedCleavages(),
                         self.GetMotifLen(), self.GetMinPercentage(),
                         self.GetMinPepLen(), self.GetPepMassTol(), self.GetPepMassUnit(), self.GetIsotopeError(),
                         self.GetNrThreads())


class Mainframe(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(750, 790))
        panel = wx.Panel(self, -1)
        self.inpanel = InputPanel(panel, -1)
        mainbox = wx.BoxSizer()
        mainbox.Add(self.inpanel, 3, wx.EXPAND)
        panel.SetSizer(mainbox)
        self.Centre()
        self.Show(True)


if __name__ == "__main__":
    app = wx.App(False)
    program = Mainframe(None, -1, "Proteogenomics MetAp activity finder")
    app.MainLoop()
