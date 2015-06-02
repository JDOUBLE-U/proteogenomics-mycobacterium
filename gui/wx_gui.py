# import modules
import sys

import wx


class CodonDialog(wx.Dialog):
    def __init__(self, parent, id, title, label, codondict):
        wx.Dialog.__init__(self, parent, id, title, size=(500, 500))
        box = wx.BoxSizer()
        grid = wx.GridSizer(16, 8, 0, 0)
        for key in sorted(codondict.keys()):
            grid.Add(self.makeSText(key), 1, wx.EXPAND | wx.ALL, 0)
            grid.Add(self.makeSText(codondict[key]), 1, wx.EXPAND | wx.ALL, 0)
        box.Add(grid, 1, wx.EXPAND)
        self.SetSizer(box)

    def makeSText(self, label):
        if len(label) == 3:
            text = wx.StaticText(self, -1, label, style=wx.BORDER_SIMPLE | wx.ALIGN_CENTRE_HORIZONTAL)
        else:
            text = wx.StaticText(self, -1, label, style=wx.BORDER_SIMPLE | wx.ALIGN_CENTRE_HORIZONTAL)
            text.SetBackgroundColour("WHITE")
        return text


class ProgPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id, style=wx.BORDER_SIMPLE)
        self.font = wx.Font(13, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        self.SetFont(self.font)
        self.mainbox = wx.BoxSizer()
        self.progtext = wx.TextCtrl(self, -1, "", style=wx.TE_READONLY | wx.TE_MULTILINE)
        #        self.progresstext.SetBackgroundColour("BLACK")
        #        self.progresstext.SetForegroundColour("WHITE")
        self.mainbox.Add(self.progtext, 1, wx.EXPAND)
        sys.stdout = self.progtext
        self.SetSizer(self.mainbox)


class InputPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self.font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        self.SetFont(self.font)
        #### Create UI elements
        ##create sizers
        self.mainbox = wx.BoxSizer()
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.mainbox.Add(self.vbox, 1, wx.EXPAND | wx.ALL, 7)
        self.SetSizer(self.mainbox)
        self.mainbox.AddSpacer(3)

        ## Select spectra folder
        self.vbox.Add(wx.StaticText(self, -1, "Spectra folder location"))
        self.specfol = wx.DirPickerCtrl(self, -1, "", "Select spectra folder")
        self.vbox.Add(self.specfol, 0, wx.EXPAND)

        ## Select output folder
        self.vbox.Add(wx.StaticText(self, -1, "Output folder location"))
        self.outfol = wx.DirPickerCtrl(self, -1, "", "Select output folder")
        self.vbox.Add(self.outfol, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)

        ## Name of output files
        self.vbox.Add(wx.StaticText(self, -1, "Name of output files"))
        self.namectrl = wx.TextCtrl(self, -1, "")
        self.vbox.Add(self.namectrl, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)

        ## Proteome location
        self.vbox.Add(wx.StaticText(self, -1, "Proteome location"))
        self.protfol = wx.DirPickerCtrl(self, -1, "", "Select proteome folder")
        self.vbox.Add(self.protfol, 0, wx.EXPAND)
        self.vbox.AddSpacer(5)

        ## Six frame or not
        self.vbox.Add(wx.StaticText(self, -1, "Use six-frame translation?"))
        self.radio_no = wx.RadioButton(self, -1, "No", style=wx.RB_GROUP)
        self.radio_yes = wx.RadioButton(self, -1, "Yes")
        self.sixchoice = 0
        self.sixbox = wx.BoxSizer()
        self.sixbox.Add(self.radio_no, 0, wx.EXPAND)
        self.sixbox.Add(wx.StaticText(self, -1, "    "), 0, wx.EXPAND)
        self.sixbox.Add(self.radio_yes, 0, wx.EXPAND)
        self.vbox.Add(self.sixbox, 0, wx.EXPAND)
        self.radio_no.Bind(wx.EVT_RADIOBUTTON, self.onSixRadio)
        self.radio_yes.Bind(wx.EVT_RADIOBUTTON, self.onSixRadio)

        ## Genome location if using six frame
        self.genotext = wx.StaticText(self, -1, "Genome location")
        self.vbox.Add(self.genotext, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.genotext)
        self.genofol = wx.DirPickerCtrl(self, -1, "", "Select Genome location")
        self.vbox.Add(self.genofol, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.genofol)
        self.Layout()

        ## startcodons
        self.codontext = wx.StaticText(self, -1, "Codons")
        self.vbox.Add(self.codontext, 0, wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.codontext)
        self.choices = ["The Standard Code", "The Vertebrate Mitochondrial Code",
                        "The Yeast Mitochondrial Code",
                        "The Mold, Protozoan, Coelenterate Mitochondrial Code, Myco/Spiroplasma Code",
                        "The Invertebrate Mitochondrial Code", "The Ciliate, Dasycladacean and Hexamita Nuclear Code",
                        "The Echinoderm and Flatworm Mitochondrial Code", "The Euplotid Nuclear Code",
                        "The Bacterial, Archaeal and Plant Plastid Code", "The Alternative Yeast Nuclear Code",
                        "The Ascidian Mitochondrial Code", "The Alternative Flatworm Mitochondrial Code",
                        "Chlorophycean Mitochondrial Code", "Trematode Mitochondrial Code",
                        "Scenedesmus obliquus Mitochondrial Code", "Thraustochytrium Mitochondrial Code",
                        "Pterobranchia Mitochondrial Code", "Candidate Division SR1 and Gracilibacteria Code"]
        self.codonmenu = wx.ComboBox(self, -1, "The Standard Code", choices=self.choices,
                                     style=wx.CB_READONLY)
        self.codonkeuze = 0
        self.constructCodontable(0)
        self.vbox.Add(self.codonmenu, 0, wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.codonmenu)
        self.codonmenu.Bind(wx.EVT_COMBOBOX, self.onCombo)
        self.viewbutton = wx.Button(self, -1, "View codon table")
        self.vbox.Add(self.viewbutton, 0, wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
        self.vbox.Hide(self.viewbutton)
        self.viewbutton.Bind(wx.EVT_BUTTON, self.onView)
        self.vbox.AddSpacer(5)

        ## Weblogo motif range
        self.vbox.Add(wx.StaticText(self, -1, "WebLogo motif range"))
        self.wlspin = wx.SpinCtrl(self, -1, "1", min=1, max=5)
        self.vbox.Add(self.wlspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Maximum error rate
        self.vbox.Add(wx.StaticText(self, -1, "Maximum error rate"))
        self.radio_1 = wx.RadioButton(self, -1, "1%", style=wx.RB_GROUP)
        self.radio_5 = wx.RadioButton(self, -1, "5%")
        self.maxerrbox = wx.BoxSizer()
        self.maxerrbox.Add(self.radio_1, 0, wx.EXPAND)
        self.maxerrbox.Add(wx.StaticText(self, -1, "   "), 0)
        self.maxerrbox.Add(self.radio_5, 0, wx.EXPAND)
        self.radio_1.Bind(wx.EVT_RADIOBUTTON, self.onErrorRadio)
        self.radio_5.Bind(wx.EVT_RADIOBUTTON, self.onErrorRadio)
        self.vbox.Add(self.maxerrbox, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Minimum length of peptides
        self.vbox.Add(wx.StaticText(self, -1, "Minimum length of peptides"))
        self.peplenspin = wx.SpinCtrl(self, -1, "5", min=5, max=100)
        self.vbox.Add(self.peplenspin, 0, wx.SHAPED)
        self.peplenspin.Bind(wx.EVT_SPINCTRL, self.onPepLenSpin)
        self.vbox.AddSpacer(5)

        ## cleavage location
        self.vbox.Add(wx.StaticText(self, -1, "Cleavage location"))
        self.cleavespin = wx.SpinCtrl(self, -1, "1", min=1, max=self.peplenspin.GetValue() - 1)
        self.vbox.Add(self.cleavespin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Digest before starting
        self.vbox.Add(wx.StaticText(self, -1, "Digest before starting?"))
        self.digradio_n = wx.RadioButton(self, -1, "No", style=wx.RB_GROUP)
        self.digradio_y = wx.RadioButton(self, -1, "Yes")
        self.digestchoice = 0
        self.digestbox = wx.BoxSizer()
        self.digestbox.Add(self.digradio_n, 0, wx.EXPAND)
        self.digestbox.Add(wx.StaticText(self, -1, "    "), 0)
        self.digestbox.Add(self.digradio_y, 0, wx.EXPAND)
        self.vbox.Add(self.digestbox, 0, wx.EXPAND | wx.ALL, 1)
        self.digradio_n.Bind(wx.EVT_RADIOBUTTON, self.onDigestRadio)
        self.digradio_y.Bind(wx.EVT_RADIOBUTTON, self.onDigestRadio)
        self.vbox.AddSpacer(5)

        ## Number of threads
        self.vbox.Add(wx.StaticText(self, -1, "Number of threads"))
        self.threadspin = wx.SpinCtrl(self, -1, "1", min=1, max=25)
        self.vbox.Add(self.threadspin, 0, wx.SHAPED)
        self.vbox.AddSpacer(5)

        ## Start button
        self.start = wx.Button(self, -1, "Start pipeline")
        self.vbox.Add(self.start, 0, wx.SHAPED | wx.ALIGN_CENTER_HORIZONTAL)
        self.start.Bind(wx.EVT_BUTTON, self.onStart)

    def onView(self, event):
        codondia = CodonDialog(self, -1, "Codon table", self.codonmenu.GetStringSelection(),
                               self.codondict)
        codondia.ShowModal()
        codondia.Destroy()

    def onPepLenSpin(self, event):
        self.wlspin.SetRange(1, self.peplenspin.GetValue())
        self.cleavespin.SetRange(1, self.peplenspin.GetValue() - 1)

    def onSixRadio(self, event):
        if event.GetEventObject().GetLabel() == "Yes":
            self.sixchoice = 1
            self.vbox.Show(self.genofol)
            self.vbox.Show(self.genotext)
            self.vbox.Show(self.codontext)
            self.vbox.Show(self.codonmenu)
            self.vbox.Show(self.viewbutton)
        else:
            self.sixchoice = 0
            self.vbox.Hide(self.genofol)
            self.vbox.Hide(self.genotext)
            self.vbox.Hide(self.codontext)
            self.vbox.Hide(self.codonmenu)
            self.vbox.Hide(self.viewbutton)
        self.Layout()

    def onCombo(self, event):
        self.codonkeuze = self.codonmenu.GetSelection()

    def onErrorRadio(self, event):
        if event.GetEventObject().GetLabel() == "1%":
            self.max_error = 0.01
        else:
            self.max_error = 0.05

    def onDigestRadio(self, event):
        if event.GetEventObject().GetLabel() == "Yes":
            self.digestchoice = 1
        else:
            self.digestchoice = 0
        print(event.GetEventObject().GetLabel())

    def onStart(self, event):
        self.checkInputs()

    def checkInputs(self):

        if self.namectrl.GetValue().isalnum():
            print("ja")
        else:
            self.errorDialog("Output filename may only contain letters and numbers and can not be empty")

    def errorDialog(self, message):
        dial = wx.MessageBox(message, 'Error', wx.OK | wx.ICON_ERROR)

    def constructCodontable(self, index):
        self.codondict = {}
        linenum = 0
        with open("..\\gui\\codontables.txt", "r") as f:
            for line in f:
                if linenum == index:
                    for pair in line.split(","):
                        pair = pair.split(":")
                        self.codondict[pair[0]] = pair[1]
                linenum += 1

    def get_spectra_folder(self):
        return self.specfol.GetPath()

    def get_output_folder(self):
        return self.outfol.GetPath()

    def get_output_name(self):
        return self.namectrl.GetValue()

    def get_proteome_folder(self):
        return self.protfol.GetPath()

    def get_sixframe_choice(self):
        return self.sixchoice

    def get_genome_folder(self):
        if self.sixchoice == 1:
            return self.genofol.GetPath()

    def get_codon_table(self):
        if self.sixchoice == 1:
            return self.codonmenu.GetSelection()

    def get_motif_range(self):
        return [1, self.wlspin.GetValue()]

    def get_max_error(self):
        return self.max_error

    def get_min_pep_pength(self):
        return self.peplenspin.GetValue()

    def get_cleavage_location(self):
        return self.cleavespin.GetValue()

    def get_digest_choice(self):
        return self.digestchoice

    def get_nr_threads(self):
        return self.threadspin.GetValue()


class Mainframe(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(1000, 750))
        panel = wx.Panel(self, -1)
        self.inpanel = InputPanel(panel, -1)
        self.progpanel = ProgPanel(panel, -1)
        mainbox = wx.BoxSizer()
        mainbox.Add(self.inpanel, 3, wx.EXPAND)
        mainbox.Add(self.progpanel, 2, wx.EXPAND | wx.LEFT, 1)
        panel.SetSizer(mainbox)
        self.Centre()
        self.Show(True)


def run_gui():
    app = wx.App(False)
    program = Mainframe(None, -1, "Proteogenomics Mycobacterium")
    app.MainLoop()

