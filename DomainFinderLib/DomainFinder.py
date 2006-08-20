# -*- coding: iso-8859-1 -*-

#
# DomainFinder 2.0.2
# (c) 1998-2001 by Konrad Hinsen
#
# For more information, see
#         http://dirac.cnrs-orleans.fr/DomainFinder/
#
# Last revision: 2005-5-2
#

from Tkinter import *
from GUIUtility import *
import tkFileDialog
from Scientific.TkWidgets.TkPlotCanvas import PolyLine, PlotGraphics, \
                                              PlotCanvas
from Scientific.TkWidgets.TkVisualizationCanvas import VisualizationCanvas, \
                                                       VisualizationGraphics, \
                                                       PolyLine3D
from MMTK.Tk.ProteinVisualization import ProteinBackboneGraphics
from Scientific.Visualization.Color import ColorByName, ColorScale
from DataGraph import DataGraph
import Numeric; N = Numeric
import operator, os, string

#
# Since the DomainFinder code is rather long, it has been split into
# two classes, Calculation and GUI. This division is for readability
# only; the two classes can't be used separately.
#

#
# The class Calculation implements all the low-level calculations.
#
class Calculation:

    def __init__(self):
        self.data = DataGraph(pdb = None, pdb2 = None,
                              sequence = None, seq2 = None, filter = 0.01,
                              nmodes_calc = 0, nmodes_keep = 0,
                              mode_selection = None,
                              rigid_limit_c = 500., cluster_factor_c = 10.,
                              rigid_limit_m = 500., cluster_factor_m = 10.,
                              universe = None, comp_conf = None,
                              test_conf = None, rms = None,
                              modes = None,
                              deformation_function = None,
                              mode_deformations = None,
                              deformation_c = None, deformation_m = None,
                              rigid_regions_c = None, clusters_c = None,
                              rigid_regions_m = None, clusters_m = None)

        self.data.defineFunction('universe', ['sequence'], self.calcUniverse)
        self.data.defineFunction('comp_conf', ['universe', 'seq2'],
                                 self.calcCompConf)
        self.data.defineFunction('test_conf',
                                 ['universe', 'comp_conf', 'filter'],
                                 self.calcTestConf)
        self.data.defineFunction('rms', ['universe', 'comp_conf', 'test_conf'],
                                 self.calcRMS)
        self.data.defineFunction('modes', ['universe',
                                           'nmodes_calc', 'nmodes_keep'],
                                 self.calcModes)
        self.data.defineFunction('deformation_function', ['universe'],
                                 self.calcDeformationFunction)
        self.data.defineFunction('mode_deformations', ['modes'],
                                 self.calcModeDeformations)
        self.data.defineFunction('deformation_c',
                                 ['deformation_function', 'test_conf'],
                                 self.calcDeformationC)
        self.data.defineFunction('deformation_m',
                                 ['universe', 'mode_deformations',
                                  'mode_selection'],
                                 self.calcDeformationM)
        self.data.defineFunction('rigid_regions_c',
                                 ['test_conf','deformation_c','rigid_limit_c'],
                                 self.calcRigidRegionsC)
        self.data.defineFunction('rigid_regions_m',
                                 ['mode_deformations', 'modes',
                                  'mode_selection', 'rigid_limit_m'],
                                 self.calcRigidRegionsM)
        self.data.defineFunction('clusters_c',
                                 ['rigid_regions_c','cluster_factor_c'],
                                 self.calcClustersC)
        self.data.defineFunction('clusters_m',
                                 ['rigid_regions_m','cluster_factor_m'],
                                 self.calcClustersM)

    def calcUniverse(self, sequence):
        self.status.set('Building model.')
        import MMTK, MMTK.Proteins
        from MMTK.ForceFields import CalphaForceField
        force_field = CalphaForceField()
        universe = MMTK.InfiniteUniverse(force_field)
        chains = []
        for chain in sequence:
            model = MMTK.Proteins.PeptideChain(chain, model='calpha')
            for i in range(len(chain)):
                model[i].name = string.capwords(chain[i].name) + \
                                `chain[i].number`
                model[i].number = chain[i].number
                chains.append(model)
        universe.protein = MMTK.Proteins.Protein(chains)
        self.status.clear()
        return universe

    def calcCompConf(self, universe, seq2):
        import MMTK
        ref_conf = MMTK.copy(universe.configuration())
        for i in range(len(universe.protein)):
            seq2[i].applyTo(universe.protein[i])
        tr, rms = universe.protein.findTransformation(ref_conf)
        universe.protein.applyTransformation(tr)
        comp_conf = MMTK.copy(universe.configuration())
        universe.setConfiguration(ref_conf)
        del ref_conf
        return comp_conf

    def calcTestConf(self, universe, comp_conf, filter):
        if filter > 0.:
            self.status.set('Applying noise filter.')
            from MMTK.Deformation import FiniteDeformationReducer
            red = FiniteDeformationReducer(universe, form='calpha',
                                          factor=1., cutoff=2.5)
            test_conf = red(comp_conf, filter)
            self.status.clear()
        else:
            test_conf = comp_conf
        return test_conf

    def calcRMS(self, universe, comp_conf, test_conf):
        rms1 = universe.rmsDifference(comp_conf)
        rms2 = universe.rmsDifference(test_conf)
        rms3 = universe.rmsDifference(comp_conf, test_conf)
        return rms1, rms2, rms3

    def calcModes(self, universe, nmodes_calc, nmodes_keep):
        self.status.set('Calculating normal modes.')
        from MMTK.NormalModes import NormalModes, \
             SparseMatrixSubspaceNormalModes
        from MMTK.FourierBasis import FourierBasis, countBasisVectors
        natoms = universe.numberOfCartesianCoordinates()
        if nmodes_calc > natoms:
            nmodes_calc = 3*natoms
        if nmodes_calc == 3*natoms:
            modes = NormalModes(universe)
            modes.calculated = 3*natoms
            modes.cutoff = None
        else:
            p1, p2 = universe.boundingBox()
            cutoff_max = (p2-p1).length()
            cutoff = 0.5*cutoff_max
            nmodes_opt = nmodes_calc
            nmodes_calc = countBasisVectors(universe, cutoff)
            while nmodes_calc > nmodes_opt:
                cutoff = cutoff + 0.1
                if cutoff > cutoff_max:
                    cutoff = cutoff_max
                    break
                nmodes_calc = countBasisVectors(universe, cutoff)
            while nmodes_calc < nmodes_opt:
                cutoff = cutoff - 0.1
                if cutoff < 0.1:
                    cutoff = 0.1
                    break
                nmodes_calc = countBasisVectors(universe, cutoff)
            basis = FourierBasis(universe, cutoff)
            basis.may_modify = 1
            modes = SparseMatrixSubspaceNormalModes(universe, basis)
            modes.calculated = len(modes)
            modes.cutoff = cutoff
        if nmodes_keep < nmodes_calc:
            modes.reduceToRange(0, nmodes_keep)
        modes.conformation = self.pdbfile.get()
        self.status.clear()
        return modes

    def calcDeformationFunction(self, universe):
        from MMTK.Deformation import FiniteDeformationFunction
        df = FiniteDeformationFunction(universe, form='calpha',
                                       factor=1., cutoff=2.5)
        return lambda v, df=df: 5.*df(v)

    def calcModeDeformations(self, modes):
        self.status.set('Calculating deformation.')
        from MMTK.Deformation import NormalizedDeformationFunction
        df = NormalizedDeformationFunction(modes.universe, factor=1.,
                                           form='calpha', cutoff=2.5)
        mode_deformations = 6*[None]
        for i in range(6, len(modes)):
            mode_deformations.append(df(modes[i]))
        self.status.clear()
        return mode_deformations

    def calcDeformationC(self, deformation_function, test_conf):
        self.status.set('Calculating deformation.')
        d = deformation_function(test_conf)
        self.status.clear()
        return d

    def calcDeformationM(self, universe, mode_deformations, mode_selection):
        if len(mode_selection) == 1:
            return mode_deformations[mode_selection[0]+6]
        else:
            import MMTK
            d = MMTK.ParticleScalar(universe)
            for i in mode_selection:
                N.add(d.array, mode_deformations[i+6].array, d.array)
            N.divide(d.array, len(mode_selection), d.array)
            return d

    def calcRigidRegionsC(self, test_conf, deformation, rigid_limit):
        self.status.set('Calculating deformation.')
        from DomainAnalysis import rigidRegionsFinite
        data, regions = rigidRegionsFinite(test_conf, deformation, rigid_limit)
        self.status.clear()
        return data, regions

    def calcRigidRegionsM(self, mode_deformations, modes, mode_selection,
                         rigid_limit):
        self.status.set('Finding rigid regions.')
        from DomainAnalysis import rigidRegions
        selected_modes = []
        selected_deformation = []
        for i in mode_selection:
            selected_modes.append(modes[i+6])
            selected_deformation.append(mode_deformations[i+6])
        data, regions = rigidRegions(selected_modes, selected_deformation,
                                     rigid_limit)
        self.status.clear()
        return data, regions

    def calcClustersC(self, rigid_regions, cluster_factor):
        self.status.set('Determining domains.')
        from DomainAnalysis import clusterAnalysis
        data, regions = rigid_regions
        if len(data) == 0:
            clusters = []
        else:
            t = data[:,:3]
            r = data[:,3:]
            dt = N.add.reduce((t[:, N.NewAxis, :]-t[N.NewAxis, :, :])**2,-1)
            dt = dt/N.add.reduce((t[:, N.NewAxis, :]+t[N.NewAxis, :, :])**2,-1)
            dr = N.add.reduce((r[:, N.NewAxis, :]-r[N.NewAxis, :, :])**2,-1)
            dr = dr/N.add.reduce((r[:, N.NewAxis, :]+r[N.NewAxis, :, :])**2,-1)
            diag = N.identity(len(data))
            dt = 1./(N.sqrt(dt)+diag)-diag
            dr = 1./(N.sqrt(dr)+diag)-diag
            clusters = clusterAnalysis(3.*dr+dt, cluster_factor)
        self.status.clear()
        return clusters

    def calcClustersM(self, rigid_regions, cluster_factor):
        self.status.set('Determining domains.')
        from DomainAnalysis import clusterAnalysis
        data, regions = rigid_regions
        if len(data) == 0:
            clusters = []
        else:
            d = 0.
            for i in range(data.shape[1]/6):
                t = data[:,6*i:6*i+3]
                r = data[:,6*i+3:6*i+6]
                dt = N.add.reduce((t[:, N.NewAxis, :]
                                   - t[N.NewAxis, :, :])**2, -1)
                dt = dt/N.add.reduce((t[:, N.NewAxis, :]
                                      + t[N.NewAxis, :, :])**2, -1)
                dr = N.add.reduce((r[:, N.NewAxis, :]
                                   - r[N.NewAxis, :, :])**2, -1)
                dr = dr/N.add.reduce((r[:, N.NewAxis, :]
                                      + r[N.NewAxis, :, :])**2, -1)
                diag = N.identity(len(data))
                dt = 1./(N.sqrt(dt)+diag)-diag
                dr = 1./(N.sqrt(dr)+diag)-diag
                d = d + 3.*dr+dt
            clusters = clusterAnalysis(d, cluster_factor)
        self.status.clear()
        return clusters

#
# The class GUI contains the Tk-based user interface
#
class GUI(Frame):

    def __init__(self, master):
        Frame.__init__(self, master)
        self.pdbfile = StringVar()
        self.pdbfile.set('')
        self.pdbfile2 = StringVar()
        self.pdbfile2.set('')
        self.rms = StringVar()
        self.rms.set('RMS distance:')
        self.mode = StringVar()
        self.mode.set('modes')
        self.ready_for_mode_analysis = 0
        self.ready_for_conformation_analysis = 0
        self.createMenu()
        self.setBindings()
        self.createMain()
        self.switchMode()
        self.activateMenus()

    def setBindings(self):
        self.setShortcut('o', lambda s=self: self.loadConformation(0))
        self.setShortcut('q', self.quit)

    def setShortcut(self, key, function):
        if mac_conventions:
            key = '<Command-%s>' % key
        else:
            key = '<Control-%s>' % key
        self.master.bind(key, lambda event, f=function: f())

    #
    # Set up the menus.
    #
    def createMenu(self):
        menu_bar = Menu(self)
        self.createFileMenu(menu_bar)
        self.createDeformationMenu(menu_bar)
        self.createDomainMenu(menu_bar)
        self.createHelpMenu(menu_bar)
        self.master.config(menu=menu_bar)

    def createFileMenu(self, menu_bar):
        if mac_conventions:
            command_key = 'Command-'
        else:
            command_key = 'Ctrl-'
        menu = Menu(menu_bar, tearoff=0)
        # entry 0
        menu.add_command(label='Load reference conformation...',
                         command=lambda object=self:
                         object.loadConformation(0),
                         accelerator=command_key+'O')
        # entry 1
        menu.add_command(label='Load comparison conformation...',
                         command=lambda object=self:
                         object.loadConformation(1))
        # entry 2
        menu.add_separator()
        # entry 3
        menu.add_command(label='Load modes...', command=self.loadModes)
        # entry 4
        menu.add_command(label='Save modes...', command=self.saveModes)
        # entry 5
        menu.add_separator()
        # entry 6
        menu.add_command(label='Show both conformations',
                         command = self.showCorrespondence)
        # entry 7
        menu.add_command(label='Show filter effect',
                         command = self.showFilter)
        # entry 8
        menu.add_command(label='RMS distances', command=self.displayRMS)
        # entry 9
        menu.add_command(label='Transition path', command=self.transitionPath)
        # entry 10
        menu.add_separator()
        # entry 11
        menu.add_command(label='Quit', command=self.quit,
                         accelerator=command_key+'Q')
        menu_bar.add_cascade(label='File', menu=menu)
        self.file_menu = menu

    def createDeformationMenu(self, menu_bar):
        menu = Menu(menu_bar, tearoff=0)
        menu.add_command(label='Show deformation',
                         command=self.showDeformation)
        menu.add_separator()
        menu.add_command(label='Write PDB file...',
                         command=self.writeDeformationPDB)
        menu.add_command(label='Write VRML file...',
                         command=self.writeDeformationVRML)
        menu_bar.add_cascade(label='Deformation', menu=menu)
        self.deformation_menu = menu

    def createDomainMenu(self, menu_bar):
        menu = Menu(menu_bar, tearoff=0)
        menu.add_command(label='Show domains', command=self.listDomains)
        menu.add_separator()
        menu.add_command(label='Write domain list...',
                         command=self.writeDomainList)
        menu.add_command(label='Write PDB file...',
                         command=self.writeDomainsPDB)
        menu.add_command(label='Write VRML file...',
                         command=self.writeDomainsVRML)
        menu.add_separator()
        menu.add_command(label='Save in MMTK format...',
                         command=self.saveDomains)
        menu_bar.add_cascade(label='Domains', menu=menu)
        self.domain_menu = menu

    def createHelpMenu(self, menu_bar):
        menu = Menu(menu_bar, tearoff=0)
        menu.add_command(label='About DomainFinder', command=self.about)
        menu_bar.add_cascade(label='Help', menu=menu)
        self.help_menu = menu

    #
    # Enable and disable menu entries according to the status.
    #
    def activateMenus(self):
        if self.mode.get() == 'modes':
            ready = self.ready_for_mode_analysis
            try: # this is for TclTkAqua, which doesn't understand "state"
                self.file_menu.entryconfig(1, state=DISABLED)
                self.file_menu.entryconfig(6, state=DISABLED)
                self.file_menu.entryconfig(7, state=DISABLED)
                self.file_menu.entryconfig(8, state=DISABLED)
                self.file_menu.entryconfig(9, state=DISABLED)
                self.file_menu.entryconfig(3, state=NORMAL)
                if ready:
                    self.file_menu.entryconfig(4, state=NORMAL)
                else:
                    self.file_menu.entryconfig(4, state=DISABLED)
            except SyntaxError:
                pass
        else:
            ready = self.ready_for_conformation_analysis
            try:
                if self.pdbfile.get():
                    self.file_menu.entryconfig(1, state=NORMAL)
                else:
                    self.file_menu.entryconfig(1, state=DISABLED)
                if ready:
                    self.file_menu.entryconfig(6, state=NORMAL)
                    self.file_menu.entryconfig(7, state=NORMAL)
                    self.file_menu.entryconfig(8, state=NORMAL)
                    self.file_menu.entryconfig(9, state=NORMAL)
                else:
                    self.file_menu.entryconfig(6, state=DISABLED)
                    self.file_menu.entryconfig(7, state=DISABLED)
                    self.file_menu.entryconfig(8, state=DISABLED)
                    self.file_menu.entryconfig(9, state=DISABLED)
                self.file_menu.entryconfig(3, state=DISABLED)
                self.file_menu.entryconfig(4, state=DISABLED)
            except TclError:
                pass
        try:
            if ready:
                self.deformation_menu.master.configure(state=NORMAL)
                self.domain_menu.master.configure(state=NORMAL)
            else:
                self.deformation_menu.master.configure(state=DISABLED)
                self.domain_menu.master.configure(state=DISABLED)
        except TclError:
            pass

    #
    # Construct the main window.
    #
    def createMain(self):
        frame = Frame(self)
        frame.pack(side=TOP, anchor=W, fill=X, padx=7, pady=7)
        Label(frame, text="Reference conformation: ").pack(side=LEFT)
        Label(frame, textvariable=self.pdbfile).pack(side=LEFT, fill=X)
        self.choice = Frame(self)
        self.choice.pack(side=TOP, fill=X)
        self.mode_switch = Radiobutton(self.choice,
                                       text="Normal mode analysis",
                                       variable = self.mode, value = "modes",
                                       borderwidth = 1, relief = RIDGE,
                                       command = self.switchMode)
        self.mode_switch.pack(side=LEFT, fill=X, ipadx=4, ipady=4)
        self.conf_switch = Radiobutton(self.choice,
                                       text="Conformation comparison",
                                       variable = self.mode,
                                       value = "conformations",
                                       borderwidth = 1, relief = RIDGE,
                                       command = self.switchMode)
        self.conf_switch.pack(side=RIGHT, fill=X, ipadx=4, ipady=4)
        self.createModesBox()
        self.createConfBox()
        self.status = StatusBar(self)
        self.status.pack(side=BOTTOM, anchor=W, fill=BOTH)

    def createModesBox(self):
        self.modes_box = Frame(self, borderwidth=2, relief=RIDGE)
        self.cluster_factor_m = FloatEntry(self.modes_box, 'domain coarseness',
                                           self.data.cluster_factor_m, 1.,None)
        self.cluster_factor_m.pack(side=BOTTOM, anchor=W, fill=X)
        self.rigid_limit_m = FloatEntry(self.modes_box,'deformation threshold',
                                        self.data.rigid_limit_m, 0., None)
        self.rigid_limit_m.pack(side=BOTTOM, anchor=W, fill=X)
        frame = Frame(self.modes_box)
        frame.pack(side=LEFT, fill=Y)
        Label(frame, text="Modes calculated:").pack(side=TOP, anchor=W,
                                                          pady=3)
        self.nmodes_calc = IntEntry(frame, '', self.data.nmodes_calc,
                                    30, None, 'modes calculated')
        self.nmodes_calc.pack(side=TOP, anchor=E)
        Label(frame, text="Modes kept for analysis:").pack(side=TOP,
                                                           anchor=W, pady=3)
        self.nmodes_keep = IntEntry(frame, '', self.data.nmodes_keep,
                                    7, None, 'modes kept for analysis')
        self.nmodes_keep.pack(side=TOP, anchor=E)
        self.mode_button = Button(frame, text="Calculate modes",
                                  command=self.calculateModes,
                                  state=DISABLED)
        self.mode_button.pack(side=TOP, padx=10, pady=15)
        frame = Frame(self.modes_box)
        frame.pack(side=RIGHT, fill=BOTH, expand=YES)
        scrollbar = Scrollbar(frame, orient=VERTICAL)
        self.modelist = Listbox(frame, relief=SUNKEN, selectmode=EXTENDED,
                                exportselection=0,
                                yscroll=scrollbar.set, width=20)
        scrollbar['command'] = self.modelist.yview
        self.modelist.pack(side=LEFT, expand=YES, fill=BOTH)
        scrollbar.pack(side=RIGHT, fill=Y)

    def createConfBox(self):
        self.conf_box = Frame(self, borderwidth=2, relief=RIDGE)
        self.cluster_factor_c = FloatEntry(self.conf_box, 'domain coarseness',
                                           self.data.cluster_factor_c, 1.,None)
        self.cluster_factor_c.pack(side=BOTTOM, anchor=W, fill=X)
        self.rigid_limit_c = FloatEntry(self.conf_box, 'deformation threshold',
                                        self.data.rigid_limit_c, 0., None)
        self.rigid_limit_c.pack(side=BOTTOM, anchor=W, fill=X)
        frame = Frame(self.conf_box)
        frame.pack(side=LEFT, fill=X)
        Label(frame, text="Comparison conformation:").pack(side=TOP, anchor=W,
                                                           pady=5)
        Label(frame, textvariable=self.pdbfile2).pack(side=TOP,anchor=W,fill=X)
        Label(frame, textvariable=self.rms).pack(side=TOP, anchor=W, pady=5)
        self.filter = Scale(self.conf_box, from_=0.1, to=0, resolution=0.001,
                            orient=VERTICAL, bd=0, label = 'Noise filter',
                            length=174)
        self.filter.pack(side=RIGHT, fill=BOTH)
        self.filter.set(self.data.filter)

    def switchMode(self):
        pass
        if self.mode.get() == 'modes':
            self.modes_box.pack(side=TOP, fill=BOTH, expand=YES)
            self.conf_box.forget()
        else:
            self.conf_box.pack(side=TOP, fill=BOTH, expand=YES)
            self.modes_box.forget()
        self.activateMenus()

    def calculateModes(self):
        self.data.nmodes_calc = self.nmodes_calc.get()
        self.data.nmodes_keep = self.nmodes_keep.get()
        self.displayModes()

    def displayModes(self):
        mode_deformations = self.data.mode_deformations
        self.modelist.delete(0, END)
        natoms = self.data.universe.numberOfCartesianCoordinates()
        for i in range(6, len(mode_deformations)):
            d_mean = N.add.reduce(mode_deformations[i].array)/natoms
            line = 'Mode %2d: %8.3f' % (i+1, d_mean)
            self.modelist.insert(END, line)
        self.nmodes_calc.set(self.data.modes.calculated)
        self.nmodes_keep.set(len(self.data.modes))
        self.ready_for_mode_analysis = 1
        self.activateMenus()

    #
    # Actions from the File menu.
    #
    def loadConformation(self, which, filename = None):
        if filename is None:
            open = tkFileDialog.Open(self,
                                     filetypes=[("PDB files", "*.pdb*"),
                                                 ("All files", "*")],
                                     title = "Choose PDB file")
            filename = open.show()
            if not filename:
                return
        self.status.set('Loading conformation.')
        from MMTK.PDB import PDBConfiguration
        sequence = PDBConfiguration(filename).peptide_chains
        self.status.clear()
        if which == 0:
            self.data.pdb = os.path.split(filename)[1]
            self.pdbfile.set(self.data.pdb)
            self.data.sequence = sequence
            natoms = reduce(operator.add, map(len, sequence))
            self.nmodes_calc.set(max(natoms/6, min(200, 3*natoms)))
            self.nmodes_keep.set(max(25, natoms/100))
            self.mode_button['state'] = NORMAL
            self.modelist.delete(0, END)
            self.data.seq2 = None
            self.data.pdb2 = None
            self.pdbfile2.set('')
            self.rms.set('RMS distance:')
            self.ready_for_mode_analysis = 0
            self.ready_for_conformation_analysis = 0
        else:
            if not self.checkSequenceCompatibility(self.data.sequence,
                                                   sequence):
                return
            self.data.seq2 = sequence
            self.data.pdb2 = os.path.split(filename)[1]
            self.pdbfile2.set(self.data.pdb2)
            rms = self.data.universe.rmsDifference(self.data.comp_conf)
            self.rms.set('RMS distance: %5.2f nm' % rms)
            self.ready_for_conformation_analysis = 1
        self.activateMenus()

    def checkSequenceCompatibility(self, seq1, seq2):
        len1 = map(len, seq1)
        len2 = map(len, seq2)
        if len1 != len2:
            if len(len1) != len(len2):
                Dialog.Dialog(self, title='Conformations incompatible',
                              text='The reference conformation has ' +
                                   `len(len1)` + ' chain(s), the comparison' +
                                   ' conformation has ' + `len(len2)` +
                                   ' chain(s).',
                              bitmap='warning', default=0,
                              strings = ('Cancel',))
            else:
                test = map(abs, map(cmp, len1, len2))
                index = test.index(1)
                Dialog.Dialog(self, title='Conformations incompatible',
                              text='Chain ' + `index+1` + ' has ' +
                                    `len1[index]` + ' residue(s) in the ' +
                                    'reference conformation, but ' +
                                    `len2[index]` + ' residue(s) in the ' +
                                    'comparison conformation.',
                              bitmap='warning', default=0,
                              strings = ('Cancel',))
            return 0
        return 1

    def loadModes(self, filename = None):
        if filename is None:
            open = tkFileDialog.Open(self,
                                     filetypes=[("Modes files", "*.modes"),
                                                 ("All files", "*")],
                                     title = "Choose modes file")
            filename = open.show()
            if not filename:
                return
        self.status.set('Loading modes.')
        import MMTK
        try:
            modes = MMTK.load(filename)
            universe = modes.universe
            calculated = modes.calculated
            conformation = modes.conformation
        except:
            Dialog.Dialog(self, title='File format error',
                          text=filename + ' is not a mode file.',
                          bitmap='warning', default=0,
                          strings = ('Cancel',))
            self.status.clear()
            return
        self.data.setValue('universe', universe)
        self.data.setValue('pdb', conformation)
        self.data.setValue('nmodes_calc', calculated)
        self.data.setValue('nmodes_keep', len(modes))
        self.data.setValue('modes', modes)
        self.status.clear()
        self.pdbfile.set(os.path.split(conformation)[-1])
        self.nmodes_calc.set(calculated)
        self.nmodes_keep.set(len(modes))
        self.mode_button['state'] = NORMAL
        self.displayModes()

    def saveModes(self, filename = None):
        if filename is None:
            open = tkFileDialog.SaveAs(filetypes=[("Mode files", "*.modes"),
                                                  ("All files", "*")],
                                       title = "Write modes file")
            filename = open.show()
            if not filename:
                return
        self.status.set('Saving modes.')
        import MMTK
        MMTK.save(self.data.modes, filename)
        self.status.clear()
        
    def showCorrespondence(self):
        self.processInput()
        universe = self.data.universe
        comp_conf = self.data.comp_conf
        self.status.set('Preparing display')
        window = Toplevel(self)
        window.title('Conformations')
        structure = VisualizationCanvas(window, width=400, height = 400,
                                        background='#BBB',
                                        relief=SUNKEN, border=2)
        structure.pack(side=TOP, fill=BOTH, expand=YES)
        gr1 = ProteinBackboneGraphics(universe.protein, None,'black')
        gr2 = ProteinBackboneGraphics(universe.protein, comp_conf,'red')
        graphics = [gr1, gr2]
        for a in universe.protein.atomList():
            p1 = a.position()
            p2 = a.position(comp_conf)
            graphics.append(PolyLine3D([p1.array, p2.array], color = 'blue'))
        structure.draw(VisualizationGraphics(graphics))
        Label(window,
              text = 'Black: reference, Red: comparison').pack(side=BOTTOM)
        self.status.clear()

    def showFilter(self):
        self.processInput()
        universe = self.data.universe
        comp_conf = self.data.comp_conf
        test_conf = self.data.test_conf
        self.status.set('Preparing display.')
        window = Toplevel(self)
        window.title('Filter')
        structure = VisualizationCanvas(window, width=400, height = 400,
                                        background='#BBB',
                                        relief=SUNKEN, border=2)
        structure.pack(side=TOP, fill=BOTH, expand=YES)
        gr1 = ProteinBackboneGraphics(universe.protein, comp_conf, 'blue')
        gr2 = ProteinBackboneGraphics(universe.protein, test_conf,'red')
        graphics = VisualizationGraphics([gr1, gr2])
        structure.draw(graphics)
        Label(window,
              text = 'Blue: not filtered, Red: filtered').pack(side=BOTTOM)
        self.status.clear()

    def displayRMS(self):
        try:
            self.processInput()
            message = ("Reference - comparison: %6.3f nm\n" +
                       "Reference - filtered comparison: %6.3f nm\n" +
                       "Comparison - filtered comparison: %6.3f nm") \
                       % (self.data.rms)
            Dialog.Dialog(self, title='RMS distances',
                          text=message, bitmap="", default=0,
                          strings = ('OK',))
        except GUIError:
            pass

    def transitionPath(self):
        import tempfile
        open = tkFileDialog.SaveAs(filetypes=[("netCDF files", "*.nc"),
                                              ("All files", "*")],
                                   title = "Write transition trajectory")
        filename = open.show()
        if not filename:
            return
        pdb1 = tempfile.mktemp('.pdb')
        pdb2 = tempfile.mktemp('.pdb')
        self.data.universe.writeToFile(pdb1)
        self.data.universe.writeToFile(pdb2, self.data.test_conf)
        os.system("TransitionPath %s %s %s 2000 1 1>/tmp/tp.log 2>&1 &"
                  % (pdb1, pdb2, filename))

    #
    # Actions from the Deformation menu
    #
    def showDeformation(self):
        try:
            self.prepareDeformation()
            from Scientific.Statistics.Histogram import Histogram
            deformation = self.prepareDeformation()
            nbins = max(25, len(deformation)/30)
            h = Histogram(deformation.array, nbins)
            self.status.set('Preparing display.')

            window = Toplevel(self)
            if self.mode.get() == 'modes':
                title = "Deformation from normal modes, deformation " + \
                        "threshold " + `self.data.rigid_limit_m`
            else:
                title = "Deformation from conformations, deformation " + \
                        "threshold "+ `self.data.rigid_limit_c`
            window.title(title)

            canvas = PlotCanvas(window, width = 300, height = 200,
                                background='#BBB', relief=SUNKEN, border=2)
            canvas.pack(side=BOTTOM, fill=BOTH, expand=YES)
            bin_width = 0.5*(h[1][0]-h[0][0])
            h2 = Numeric.repeat(h.array, len(h)*[2])
            h2[::2, 0] = h[:, 0]-bin_width
            h2[1::2, 0] = h[:, 0]+bin_width
            lines = PolyLine(h2, color='red')
            canvas.draw(lines, 'automatic', 'automatic')

            Label(window,text='Deformation histogram').pack(side=BOTTOM,fill=X)

            structure = VisualizationCanvas(window, width=500, height = 300,
                                            background='#BBB',
                                            relief=SUNKEN, border=2)
            structure.pack(side=LEFT, fill=BOTH, expand=YES)
            atom_color = {}
            for atom in self.data.universe.atomList():
                atom_color[atom.index] = rgbcolor(atom.color)
            graphics = ProteinBackboneGraphics(self.data.universe.protein,
                                               None, atom_color)
            structure.draw(graphics)

            self.status.clear()
        except GUIError:
            pass

    def writeDeformationPDB(self):
        try:
            self.prepareDeformation()
            open = tkFileDialog.SaveAs(filetypes=[("PDB files", "*.pdb*"),
                                                  ("All files", "*")],
                                       title = "Write deformation PDB file")
            filename = open.show()
            if filename:
                self.status.set('Writing PDB file.')
                self.data.universe.writeToFile(filename, None, 'pdb')
                self.status.clear()
        except GUIError:
            pass

    def writeDeformationVRML(self):
        try:
            self.prepareDeformation()
            open = tkFileDialog.SaveAs(filetypes=[("VRML files", "*.wrl*"),
                                                  ("All files", "*")],
                                       title = "Write deformation VRML file")
            filename = open.show()
            if filename:
                self.data.universe.writeToFile(filename, None,
                                               'vrml.wireframe')
        except GUIError:
            pass
        
    def prepareDeformation(self):
        self.processInput()
        if self.mode.get() == 'modes':
            deformation = self.data.deformation_m
            limit = self.data.rigid_limit_m
        else:
            deformation = self.data.deformation_c
            limit = self.data.rigid_limit_c
        factor = 99./(2.5*limit)
        scale = ColorScale(2.5*limit)
        for a in self.data.universe.atomList():
            v = deformation[a]
            if v > 2.*limit:
                v = 2.5*limit
            a.color = scale(v)
            a.temperature_factor = factor*v
        return deformation

    #
    # Actions from the Domains menu
    #
    def listDomains(self):
        try:
            domains, colors = self.getDomainsAndColors()
            if self.mode.get() == 'modes':
                title = "Domains from normal modes, coarseness " + \
                        `self.data.cluster_factor_m` + ", deformation " + \
                        "threshold " + `self.data.rigid_limit_m`
            else:
                title = "Domains from conformations, coarseness " + \
                        `self.data.cluster_factor_c` + ", deformation " + \
                        "threshold " + `self.data.rigid_limit_c`
            self.status.set('Preparing display.')
            DomainList(self, self.data.universe.protein,
                       domains, colors, title)
            self.status.clear()
        except GUIError:
            pass

    def writeDomainList(self):
        from Scientific.IO.TextFile import TextFile
        try:
            open = tkFileDialog.SaveAs(filetypes=[("Text files", "*.txt"),
                                                  ("All files", "*")],
                                       title = "Write domain list")
            filename = open.show()
            if filename:
                domains, colors = self.getDomainsAndColors()
                file = TextFile(filename, 'w')
                for i in range(len(domains)):
                    d = domains[i]
                    line = 'Domain ' + `i+1` + ' (' + colors[i] + '), ' + \
                           `len(d)` + ' residues, similarity ' + \
                           `d.similarity[-1]` + '\n'
                    file.write(line)
                    file.write((len(line)-1)*'=' + '\n\n')
                    for line in residueDescription(d,
                                                   self.data.universe.protein):
                        file.write(line + '\n')
                    file.write('\n\n')
                file.close()
        except GUIError:
            pass

    def writeDomainsPDB(self):
        try:
            self.prepareDeformation()
            self.prepareDomains()
            open = tkFileDialog.SaveAs(filetypes=[("PDB files", "*.pdb*"),
                                                  ("All files", "*")],
                                       title = "Write domain PDB file")
            filename = open.show()
            if filename:
                self.status.set('Writing PDB file.')
                self.data.universe.writeToFile(filename, None, 'pdb')
                self.status.clear()
        except GUIError:
            pass

    def writeDomainsVRML(self):
        try:
            self.prepareDomains()
            open = tkFileDialog.SaveAs(filetypes=[("VRML files", "*.wrl*"),
                                                  ("All files", "*")],
                                       title = "Write domain VRML file")
            filename = open.show()
            if filename:
                self.status.set('Writing VRML file.')
                self.data.universe.writeToFile(filename,None,'vrml.wireframe')
                self.status.clear()
        except GUIError:
            pass

    def saveDomains(self):
        try:
            domains, rigid_regions, rbdata = self.findDomains()
            open = tkFileDialog.SaveAs(filetypes=[("Domain files", "*.domains"),
                                                  ("All files", "*")],
                                       title = "Write domain file (MMTK format)")
            filename = open.show()
            if filename:
                self.status.set('Saving domains.')
                import MMTK
                domain_objects = []
                for d in domains:
                    do = MMTK.Collection(list(N.repeat(rigid_regions, d[1])))
                    do.similarity = d[0]
                    do.rbdata = N.repeat(rbdata, d[1])
                    domain_objects.append(do)
                data = {}
                data['conformation'] = self.data.pdb
                data['universe'] = self.data.universe
                data['domains'] = domain_objects
                if self.mode.get() == 'modes':
                    data['mode_selection'] = self.data.mode_selection
                    data['nmodes'] = self.data.nmodes_calc
                    data['deformation_threshold'] = self.data.rigid_limit_m
                    data['domain_coarseness'] = self.data.cluster_factor_m
                    data['deformation'] = self.data.deformation_m
                else:
                    data['comparison'] = self.data.pdb2
                    data['filter'] = self.data.filter
                    data['deformation_limit'] = self.data.rigid_limit_c
                    data['domain_coarseness'] = self.data.cluster_factor_c
                    data['deformation'] = self.data.deformation_c
                MMTK.save(data, filename)
                self.status.clear()
        except GUIError:
            pass

    def getDomainsAndColors(self):
        domains, regions, rbdata = self.findDomains()
        import MMTK
        domain_objects = []
        for i in range(len(domains)):
            d = domains[i]
            do = MMTK.Collection(list(N.repeat(regions, d[1])))
            do.similarity = d[0]
            try:
                do.similarity[0]
            except TypeError:
                do.similarity = [d[0]]
            do.rbdata = N.repeat(rbdata, d[1])
            domain_objects.append(do)
        col = self.domain_colors
        while len(col) < len(domains):
            col = col + col
        col = col[:len(domains)]
        return domain_objects, col

    def prepareDomains(self):
        domains, regions, rbdata = self.findDomains()
        col = self.domain_colors
        while len(col) < len(domains):
            col = col + col
        col = col[:len(domains)]
        for a in self.data.universe.atomList():
            a.color = 'black'
            a.occupancy = 0
        import MMTK
        for a in MMTK.Collection(list(regions)):
            a.color = 'white'
        for i in range(len(domains)):
            d = domains[i]
            c = col[i]
            do = MMTK.Collection(list(N.repeat(regions, d[1])))
            do.similarity = d[0]
            for a in do.atomList():
                a.color = c
                a.occupancy = i+1

    domain_colors = ['yellow', 'red', 'green', 'cyan', 'blue',
                     'magenta', 'orange', 'dark green', 'brown',
                     'violet']

    def findDomains(self):
        self.processInput()
        if self.mode.get() == 'modes':
            rbdata, regions = self.data.rigid_regions_m
            domains = self.data.clusters_m
        else:
            rbdata, regions = self.data.rigid_regions_c
            domains = self.data.clusters_c
        return domains, regions, rbdata

    #
    # Actions from the Help menu.
    #
    def about(self):
        window = Toplevel(self)
        window.title('About DomainFinder')
        text = 'DomainFinder 2.0.2\n' + \
               '\n' + \
               '(c) 1998-2005 by Konrad Hinsen\n' + \
               '\n' + \
               'Laboratoire Léon Brillouin (CEA-CNRS)\n' + \
               'CEA Saclay\n' + \
               '91191 Gif sur Yvette Cedex\n' + \
               'France\n' + \
               '\n' + \
               'E-Mail: khinsen@cea.fr\n' + \
               '\n' + \
               'http://dirac.cnrs-orleans.fr/DomainFinder/'
        Label(window, text=text).pack(side=TOP)

    #
    # Process input data before performing any action.
    #
    def processInput(self):
        self.master.update_idletasks()
        if not self.pdbfile.get():
            Dialog.Dialog(self, title='No input structure',
                          text = 'You must first load a ' +
                          'reference conformation', bitmap = 'warning',
                          default = 0, strings = ('Cancel',))
            raise GUIError
        if self.mode.get() == 'modes':
            mode_selection = self.modelist.curselection()
            try: mode_selection = map(string.atoi, mode_selection)
            except ValueError: pass
            if not mode_selection:
                Dialog.Dialog(self, title='No modes selected',
                              text='Please select at least one mode.',
                              bitmap='warning', default=0,
                              strings = ('Cancel',))
                raise GUIError
            self.data.mode_selection = mode_selection
            try:
                self.data.rigid_limit_m = self.rigid_limit_m.get()
                self.data.cluster_factor_m = self.cluster_factor_m.get()
            except ValueError:
                raise GUIError
        else:
            if not self.pdbfile2.get():
                Dialog.Dialog(self, title='No input structure',
                              text = 'You must first load a ' +
                              'comparison conformation', bitmap = 'warning',
                              default = 0, strings = ('Cancel',))
                raise GUIError
            self.data.filter = self.filter.get()
            try:
                self.data.rigid_limit_c = self.rigid_limit_c.get()
                self.data.cluster_factor_c = self.cluster_factor_c.get()
            except ValueError:
                raise GUIError


#
# The DomainFinder class is simply a combination of GUI and calculation
#
class DomainFinder(GUI, Calculation):

    def __init__(self, master):
        Calculation.__init__(self)
        GUI.__init__(self, master)

#
# The following code implements the separate windows used to display
# domain lists and residue lists.
#
def residueDescription(residues, protein):
    chain_number = {}
    for i in range(len(protein)):
        chain_number[protein[i]] = i+1
    residues = map(lambda r: r.parent, residues)
    real_numbers = 1
    try:
        index = map(lambda r, cn=chain_number: (cn[r.parent], r.number, r),
                    residues)
    except AttributeError:
        index = map(lambda r, cn=chain_number: (cn[r.parent],
                                                r.sequence_number, r),
                    residues)
        real_numbers = 0
    index.sort(lambda a, b: cmp(a[0], b[0]) or cmp(a[1], b[1]))
    lines = []
    if real_numbers:
        while index:
            n = 1
            while n < len(index) and index[n][0] == index[0][0] and \
                  index[n][1] == index[n-1][1]+1:
                n = n + 1
            line = 'Chain ' + `index[0][0]` + ', ' + index[0][2].name
            if n > 1:
                line = line + ' - ' + index[n-1][2].name
            lines.append(line)
            index = index[n:]
    else:
        for i, j, r in index:
            line = 'Chain ' + `i` + ', residue ' + `j` + \
                   ' (' + r.name[:3] + ')'
            lines.append(line)
    return lines


class DomainList(Toplevel):

    def __init__(self, master, protein, domains, colors, title = 'Domains'):
        Toplevel.__init__(self, master)
        self.title(title)

        self.protein = protein
        self.domains = domains
        self.colors = colors

        plot = PlotCanvas(self, width=300, height = 180,
                          background='#BBB', relief=SUNKEN, border=2)
        plot.pack(side=BOTTOM, fill=BOTH, expand=YES)

        self.structure = VisualizationCanvas(self, width=300, height = 300,
                                             background='#BBB',
                                             relief=SUNKEN, border=2)
        self.structure.pack(side=LEFT, fill=BOTH, expand=YES)

        scrollbar = Scrollbar(self, orient=VERTICAL)
        listbox = Listbox(self, relief=SUNKEN, yscroll=scrollbar.set,
                          width = 40)
        scrollbar['command'] = listbox.yview
        scrollbar.pack(side=RIGHT, fill=Y)
        listbox.pack(side=RIGHT, expand=1, fill=BOTH)

        atom_color = {}
        for atom in protein.atomList():
            atom_color[atom.index] = 'black'
        plotobjects = []
        lower = 0.
        upper = 0.
        for i in range(len(domains)):
            d = domains[i]
            line = 'Domain %d (%s), %d residues' % (i+1, colors[i], len(d))
            listbox.insert(END, line)
            color = rgbcolor(ColorByName(colors[i]))
            for atom in d.atomList():
                atom_color[atom.index] = color
            for point in d.rbdata:
                point = N.transpose(N.array([N.arange(len(point)), point]))
                plotobjects.append(PolyLine(point, color=color))
            lower = min(lower, N.minimum.reduce(N.ravel(d.rbdata)))
            upper = max(upper, N.maximum.reduce(N.ravel(d.rbdata)))
        listbox.bindtags(listbox.bindtags()+('DomainSelection',))
        listbox.bind_class('DomainSelection', '<1>', self.selectDomain)
        for i in range(domains[0].rbdata.shape[1]):
            plotobjects.append(PolyLine([(i,lower), (i,upper)], color='black'))
        plot.draw(PlotGraphics(plotobjects))
        self.plotrange = (lower, upper)
        self.structure.draw(ProteinBackboneGraphics(protein, None, atom_color))

    def selectDomain(self, event):
        self.configure(cursor='watch')
        self.update_idletasks()
        number = event.widget.nearest(event.y)
        ResidueList(None, self.protein, self.domains[number],
                    self.colors[number], self.plotrange,
                    self.structure, 'Domain '+`number+1`)
        self.configure(cursor='top_left_arrow')


class ResidueList(Toplevel):

    def __init__(self, master, protein, domain, color, plotrange,
                 structure_canvas, title):
        Toplevel.__init__(self, master)
        self.title(title)

        canvas2 = PlotCanvas(self, width=200, height = 50,
                             background='#BBB', relief=SUNKEN, border=2)
        canvas2.pack(side=BOTTOM, fill=X)
        coarseness = domain.similarity[-1]/domain.similarity[0]
        text = 'Similarity %.0f-%.0f, Coarseness 1-%.0f' % \
               (domain.similarity[-1], domain.similarity[0], coarseness)
        Label(self, text = text).pack(side=BOTTOM)

        canvas1 = PlotCanvas(self, width=200, height = 180,
                             background='#BBB', relief=SUNKEN, border=2)
        canvas1.pack(side=BOTTOM, fill=BOTH, expand=YES)

        scrollbar = Scrollbar(self, orient=VERTICAL)
        listbox = Listbox(self, relief=SUNKEN, yscroll=scrollbar.set, width=30)
        scrollbar['command'] = listbox.yview
        scrollbar.pack(side=RIGHT, fill=Y)
        listbox.pack(side=RIGHT, expand=1, fill=BOTH)
        for line in residueDescription(domain, protein):
            listbox.insert(AtEnd(), line)

        structure = VisualizationCanvas(self, width=300, height = 300,
                                        background='#BBB',
                                        relief=SUNKEN, border=2)
        structure.pack(side=LEFT, fill=BOTH, expand=YES)

        atom_color = {}
        for atom in protein.atomList():
            atom_color[atom.index] = 'black'
        for atom in domain.atomList():
            atom_color[atom.index] = color

        plotobjects = []
        color = rgbcolor(ColorByName(color))
        lower, upper = plotrange
        for point in domain.rbdata:
            point = N.transpose(N.array([N.arange(len(point)), point]))
            plotobjects.append(PolyLine(point, color=color))
        for i in range(domain.rbdata.shape[1]):
            plotobjects.append(PolyLine([(i,lower), (i,upper)], color='black'))
        canvas1.draw(PlotGraphics(plotobjects))
        plotobjects = []
        max_sim = domain.similarity[-1]
        for s in domain.similarity:
            x = max_sim/s
            plotobjects.append(PolyLine([(x, 0.), (x, 1.)], color=color))
        canvas2.draw(PlotGraphics(plotobjects), xaxis='automatic')
        structure.copyViewpointFrom(structure_canvas)
        structure.draw(ProteinBackboneGraphics(protein, None, atom_color))

#
# Convert color names from visualization to Tk TGB codes
#
def rgbcolor(color):
    rgb = map(lambda x: ('0'+hex(int(255*x))[2:])[-2:], color.rgb)
    return '#' + string.upper(string.join(rgb, ''))

#
# And finally: the main program...
#
def run(filename1 = None, filename2 = None):

    global mac_conventions
    
    root = Tk()
    try:
        root.tk.call('console','hide')
    except TclError:
        pass
    
    if sys.platform == 'darwin':
        mac_conventions = 'X11' not in root.winfo_server()
        if 0:
            # Unfortunately the Scale widget doesn't work with tile enabled.
            # Perhaps this will be fixed, until then we have to live with
            # standard Tk look.
            try:
                root.tk.call('package', 'require', 'tile')
                root.tk.call('namespace', 'import', '-force', 'ttk::*')
                root.tk.call('tile::setTheme', 'aqua')
            except TclError:
                pass
    else:
        mac_conventions = 0
    
    root.title('DomainFinder 2.0.2')
    root.resizable(width=NO,height=NO)
    app = DomainFinder(root)
    app.pack(side=TOP, anchor=W, fill=BOTH, expand=YES)

    if filename1:
        app.loadConformation(0, filename1)
    if filename2:
        app.loadConformation(1, filename2)
        app.conf_switch.select()
        app.switchMode()
    app.mainloop()
    
if __name__ == '__main__':

    run()
