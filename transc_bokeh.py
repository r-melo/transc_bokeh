#!/usr/bin/env python2
import platform
import sys
import threading
try:
        # these libs work with python2
        import Tkinter, tkFileDialog, tkSimpleDialog, tkMessageBox
except:
        # these libs work with python3
        import tkinter as Tkinter
        import tkinter.filedialog as tkFileDialog
        import tkinter.simpledialog as tkSimpleDialog
        import tkinter.messagebox as tkMessageBox
import ttk 
import subprocess

listOfFilesIdentifiers = ["ordering", "dictionary", "expression", "association", "orderingModularity", "associationModularity", "associationNetwork", "orderingMatrix", "associationMatrix", "orderingOccupation", "associationOccupation", "orderingEnrichment", "dictionaryEnrichment", "filesList", "expressionStatistics", "expressionGeneSetDiffExp", "filesListGeneSetDiffExp", "dictionaryGeneSetDiffExp", "expressionDendrogram", "orderingTranslate", "expressionTranslate", "dictionaryTranslate","pvalFigure","goFigure","keggFigure","dictFigure"]
filesNames = {}
labelFileName = {}
labelFileLines = {}

headerEntryStatistics = {}
headerLabelStatistics = {}
headerEntryGeneSetDiffExp = {}
headerLabelGeneSetDiffExp = {}


# Get the folder path where this .py file is and where the binaries must be
binPath = subprocess.os.path.dirname(subprocess.os.path.abspath(__file__))
if (subprocess.mswindows == True):
        binPath = binPath + '\\'
else:
        binPath = binPath + '/'


class Application(Tkinter.Frame):
        def __init__(self, root):
                self.is_plotting = False
                self.plot_th = None
                Tkinter.Frame.__init__(self, root)


                #################################
                ##            Tabs             ##
                #################################

                notebookFrame = ttk.Notebook(self)
                notebookFrame.grid()

                ordererFrame = ttk.Frame(notebookFrame)
                modularityFrame = ttk.Frame(notebookFrame)
                networkFrame = ttk.Frame(notebookFrame)
                matrixFrame = ttk.Frame(notebookFrame)
                occupationFrame = ttk.Frame(notebookFrame)
                transcriptogramerFrame = ttk.Frame(notebookFrame)
                termEnrichmentFrame = ttk.Frame(notebookFrame)
                statisticsFrame = ttk.Frame(notebookFrame)
                geneSetDiffExpFrame = ttk.Frame(notebookFrame)
                dendrogramFrame = ttk.Frame(notebookFrame)
                translateFrame = ttk.Frame(notebookFrame)
                figureFrame = ttk.Frame(notebookFrame)
                notebookFrame.add(ordererFrame, text="Orderer")
                notebookFrame.add(modularityFrame, text="Ordering properties")
#                notebookFrame.add(networkFrame, text="Network properties")
                notebookFrame.add(matrixFrame, text="Matrix")
                notebookFrame.add(occupationFrame, text="Occupation Level")
                notebookFrame.add(transcriptogramerFrame, text="Transcriptogramer")
                notebookFrame.add(termEnrichmentFrame, text="Term Enrichment")
                notebookFrame.add(statisticsFrame, text="Statistics")
                notebookFrame.add(geneSetDiffExpFrame, text="Gene Set Differential Expression")
                notebookFrame.add(dendrogramFrame, text="Dendrogram")
                notebookFrame.add(translateFrame, text="Translate")
                notebookFrame.add(figureFrame, text="Figure")

                ordererFrame.columnconfigure(0, weight=1)
                modularityFrame.columnconfigure(0, weight=1)
                networkFrame.columnconfigure(0, weight=1)
                matrixFrame.columnconfigure(0, weight=1)
                occupationFrame.columnconfigure(0, weight=1)
                transcriptogramerFrame.columnconfigure(0, weight=1)
                termEnrichmentFrame.columnconfigure(0, weight=1)
                statisticsFrame.columnconfigure(0, weight=1)
                geneSetDiffExpFrame.columnconfigure(0, weight=1)
                dendrogramFrame.columnconfigure(0, weight=1)
                translateFrame.columnconfigure(0, weigh=1)
                figureFrame.columnconfigure(0, weigh=1)

                for fileIdentifier in listOfFilesIdentifiers:
                        filesNames[fileIdentifier] = "NONE"
                        labelFileName[fileIdentifier] = Tkinter.Label()
                        labelFileLines[fileIdentifier] = Tkinter.Label()


                #################################
                ##           Orderer           ##
                #################################

                ## Association                  ##
                associationLabelFrame = Tkinter.LabelFrame(ordererFrame,  text=" 1 - Choose association file: ")
                associationLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("association", associationLabelFrame)

                ##                Isothermal steps                ##
                isothermalLabelFrame = Tkinter.LabelFrame(ordererFrame,  text=" 2 - Choose number of isothermal steps: ")
                isothermalLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                isothermalSteps=Tkinter.StringVar()
                isothermalSteps.set("100")
                isothermalEntry = Tkinter.Entry(isothermalLabelFrame, validate="key", textvariable=isothermalSteps)
                isothermalEntry['validatecommand'] = (self.register(self.validateIntegerKey), '%S', '%d')
                isothermalEntry.grid(padx=5, pady=5)

                ##          Monte Carlo Steps         ##
                monteCarloLabelFrame = Tkinter.LabelFrame(ordererFrame,  text=" 3 - Choose number of Monte Carlo steps: ")
                monteCarloLabelFrame.grid(row=2, padx=5, pady=5, sticky='WE')
                monteCarloSteps = Tkinter.StringVar()
                monteCarloSteps.set("2000")
                monteCarloEntry = Tkinter.Entry(monteCarloLabelFrame, validate="key", textvariable=monteCarloSteps)
                monteCarloEntry['validatecommand'] = (self.register(self.validateIntegerKey), '%S', '%d')
                monteCarloEntry.grid(padx=5, pady=5)

                ##          Cooling Factor         ##
                coolingLabelFrame = Tkinter.LabelFrame(ordererFrame,  text=" 4 - Choose cooling factor: ")
                coolingLabelFrame.grid(row=3, padx=5, pady=5, sticky='WE')
                coolingFactor = Tkinter.StringVar()
                coolingFactor.set("0.5")
                coolingEntry = Tkinter.Entry(coolingLabelFrame, validate="key", textvariable=coolingFactor)
                coolingEntry['validatecommand'] = (self.register(self.validateFloatKey), '%S', '%d')
                coolingEntry.grid(padx=5, pady=5)

                ##                                Alpha         ##
                alphaLabelFrame = Tkinter.LabelFrame(ordererFrame,  text=" 5 - Choose alpha value: ")
                alphaLabelFrame.grid(row=4, padx=5, pady=5, sticky='WE')
                alpha = Tkinter.StringVar()
                alpha.set("1.0")
                alphaEntry = Tkinter.Entry(alphaLabelFrame, validate="key", textvariable=alpha)
                alphaEntry['validatecommand'] = (self.register(self.validateFloatKey), '%S', '%d')
                alphaEntry.grid(padx=5, pady=5)

                ##          Initial Weigth                ##
                percentualWeigthFrame = Tkinter.LabelFrame(ordererFrame,  text=" 6 - Choose percentual energy for initial temperature: ")
                percentualWeigthFrame.grid(row=5, padx=5, pady=5, sticky='WE')
                percentualWeigth = Tkinter.StringVar()
                percentualWeigth.set("0.0001")
                percentualWeigthEntry = Tkinter.Entry(percentualWeigthFrame, validate="key", textvariable=percentualWeigth)
                percentualWeigthEntry['validatecommand'] = (self.register(self.validateIntegerKey), '%S', '%d')
                percentualWeigthEntry.grid(padx=5, pady=5)

                ##    Program    ##
                Tkinter.Button(ordererFrame, text='Order', command=lambda: self.runOrderer(isothermalEntry.get(), monteCarloEntry.get(), coolingEntry.get(), alphaEntry.get(), percentualWeigthEntry.get() ) ).grid()


                #################################
                ##      Ordering Properties    ##
                #################################

                ##    Ordering                         ##
                orderingModularityLabelFrame = Tkinter.LabelFrame(modularityFrame,  text=" 1 - Choose ordering file: ")
                orderingModularityLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("orderingModularity", orderingModularityLabelFrame)

                ## Association                  ##
                associationModularityLabelFrame = Tkinter.LabelFrame(modularityFrame,  text=" 2 - Choose association file: ")
                associationModularityLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                self.createFileInformation("associationModularity", associationModularityLabelFrame)

                ##                                Radius    ##
                radiusModularityLabelFrame = Tkinter.LabelFrame(modularityFrame,  text=" 3 - Choose radius: ")
                radiusModularityLabelFrame.grid(row=3, padx=5, pady=5, sticky='WE')
                radiusModularity = Tkinter.StringVar()
                radiusModularity.set("4.0")
                radiusModularityEntry = Tkinter.Entry(radiusModularityLabelFrame, validate="key", textvariable=radiusModularity)
                radiusModularityEntry['validatecommand'] = (self.register(self.validateFloatKey), '%S', '%d')
                radiusModularityEntry.grid(padx=5, pady=5)

                ##    Program    ##
                Tkinter.Button(modularityFrame, text='Calculate modularity', command=lambda: self.runProperties(radiusModularityEntry.get())).grid()


                #################################
                ##     Network Properties      ##
                #################################

                ## Association                  ##
                associationNetworkLabelFrame = Tkinter.LabelFrame(networkFrame,  text=" 1 - Choose association file: ")
                associationNetworkLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("associationNetwork", associationNetworkLabelFrame)

                ##    Program    ##
                Tkinter.Button(networkFrame, text='Calculate network properties', command=lambda: self.runNetwork()).grid()


                #################################
                ##             Matrix          ##
                #################################

                ##    Ordering                         ##
                orderingMatrixLabelFrame = Tkinter.LabelFrame(matrixFrame,  text=" 1 - Choose ordering file: ")
                orderingMatrixLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("orderingMatrix", orderingMatrixLabelFrame)

                ## Association                  ##
                associationMatrixLabelFrame = Tkinter.LabelFrame(matrixFrame,  text=" 2 - Choose association file: ")
                associationMatrixLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                self.createFileInformation("associationMatrix", associationMatrixLabelFrame)

                ##    Program    ##
                Tkinter.Button(matrixFrame, text='Calculate Matrix', command=lambda: self.runMatrix()).grid()



                #################################
                ##           Occupation        ##
                #################################

                ##    Ordering                         ##
                orderingOccupationLabelFrame = Tkinter.LabelFrame(occupationFrame,  text=" 1 - Choose ordering file: ")
                orderingOccupationLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("orderingOccupation", orderingOccupationLabelFrame)

                ## Association                  ##
                associationOccupationLabelFrame = Tkinter.LabelFrame(occupationFrame,  text=" 2 - Choose association file: ")
                associationOccupationLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                self.createFileInformation("associationOccupation", associationOccupationLabelFrame)

                ##    Program    ##
                Tkinter.Button(occupationFrame, text='Calculate occupation level', command=lambda: self.runOccupation()).grid()



                #################################
                ##     Transcriptogramer       ##
                #################################

                ##    Ordering ##
                orderingLabelFrame = Tkinter.LabelFrame(transcriptogramerFrame,  text=" 1 - Choose ordering file: ")
                orderingLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("ordering", orderingLabelFrame)

                ## Dictionary ##
                dictionaryLabelFrame = Tkinter.LabelFrame(transcriptogramerFrame,  text=" 2 - Choose dictionary file: ")
                dictionaryLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                self.createFileInformation("dictionary", dictionaryLabelFrame)

                ## Expression ##
                expressionLabelFrame = Tkinter.LabelFrame(transcriptogramerFrame,  text=" 3 - Choose expression file: ")
                expressionLabelFrame.grid(row=2, padx=5, pady=5, sticky='WE')
                self.createFileInformation("expression", expressionLabelFrame)

                ## Radius ##
                radiusLabelFrame = Tkinter.LabelFrame(transcriptogramerFrame,  text=" 4 - Choose radius: ")
                radiusLabelFrame.grid(row=3, padx=5, pady=5, sticky='WE')
                radius = Tkinter.StringVar()
                radius.set("4.0")
                radiusEntry = Tkinter.Entry(radiusLabelFrame, validate="key", textvariable=radius)
                radiusEntry['validatecommand'] = (self.register(self.validateFloatKey), '%S', '%d')
                radiusEntry.grid(padx=5, pady=5)

                ## Program ##
                Tkinter.Button(transcriptogramerFrame, text='Make transcriptogram', command=lambda: self.runTranscriptogramer(radiusEntry.get())).grid()


                #################################
                ##       Term Enrichment       ##
                #################################

                ## Ordering ##
                orderingEnrichmentLabelFrame = Tkinter.LabelFrame(termEnrichmentFrame,  text=" 1 - Choose ordering file: ")
                orderingEnrichmentLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("orderingEnrichment", orderingEnrichmentLabelFrame)

                ## Dictionary ##
                dictionaryEnrichmentLabelFrame = Tkinter.LabelFrame(termEnrichmentFrame,  text=" 2 - Choose dictionary file: ")
                dictionaryEnrichmentLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                self.createFileInformation("dictionaryEnrichment", dictionaryEnrichmentLabelFrame)

                ## Files List ##
                filesListLabelFrame = Tkinter.LabelFrame(termEnrichmentFrame,  text=" 3 - Choose files list: ")
                filesListLabelFrame.grid(row=2, padx=5, pady=5, sticky='WE')
                self.createFileInformation("filesList", filesListLabelFrame)

                ## Radius ##
                radiusEnrichmentLabelFrame = Tkinter.LabelFrame(termEnrichmentFrame,  text=" 4 - Choose radius: ")
                radiusEnrichmentLabelFrame.grid(row=3, padx=5, pady=5, sticky='WE')
                radiusEnrichment = Tkinter.StringVar()
                radiusEnrichment.set("4.0")
                radiusEnrichmentEntry = Tkinter.Entry(radiusEnrichmentLabelFrame, validate="key", textvariable=radiusEnrichment)
                radiusEnrichmentEntry['validatecommand'] = (self.register(self.validateFloatKey), '%S', '%d')
                radiusEnrichmentEntry.grid(padx=5, pady=5)

                ## Program ##
                Tkinter.Button(termEnrichmentFrame, text='Calculate term enrichment', command=lambda: self.runTermEnrichment(radiusEnrichmentEntry.get())).grid()

                #################################
                ##          Statistics         ##
                #################################

                ## Expression ##
                expressionStatisticsLabelFrame = Tkinter.LabelFrame(statisticsFrame,  text=" 1 - Choose transcriptogram file: ")
                expressionStatisticsLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("expressionStatistics", expressionStatisticsLabelFrame)

                ## Number ofpermutatins for frdr calculation ##
                permutationLabelFrame = Tkinter.LabelFrame(statisticsFrame,  text=" 2 - Choose the number of permutations for the false discovery rate calculation: ")
                permutationLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                permutation = Tkinter.StringVar()
                permutation.set("1000")
                permutationEntry = Tkinter.Entry(permutationLabelFrame, validate="key", textvariable=permutation)
                permutationEntry['validatecommand'] = (self.register(self.validateIntegerKey), '%S', '%d')
                permutationEntry.grid(padx=5, pady=5)

                ## Headers ##
                headerLabelStatisticsFrame = Tkinter.LabelFrame(statisticsFrame,  text=" 3 - Choose the identifiers of the groups: ")
                headerLabelStatisticsFrame.grid(row=2, padx=5, pady=5, sticky='WE')
                headerLabelStatistics[0] = Tkinter.Label(headerLabelStatisticsFrame,text="Group 1: ")
                headerLabelStatistics[0].grid(row=0,column=0,padx=5, pady=5)
                headerEntryStatistics[0] = Tkinter.Entry(headerLabelStatisticsFrame)
                headerEntryStatistics[0].grid(row=0,column=1,padx=5, pady=5)
                insertStatisticsButton = Tkinter.Button(headerLabelStatisticsFrame, text='Add', command=lambda: self.createField(headerLabelStatisticsFrame,headerEntryStatistics,headerLabelStatistics,insertStatisticsButton))
                insertStatisticsButton.grid()

                ## Program ##
                Tkinter.Button(statisticsFrame, text='Calculate statistics', command=lambda: self.runStatistics(permutationEntry.get())).grid()



                #################################
                ## Gene Set Differential Expression ##
                #################################

                ## Expression ##
                expressionGeneSetDiffExpLabelFrame = Tkinter.LabelFrame(geneSetDiffExpFrame, text=" 1 - Choose expression file: ")
                expressionGeneSetDiffExpLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("expressionGeneSetDiffExp", expressionGeneSetDiffExpLabelFrame)

                ## Files List ##
                filesListGeneSetDiffExpLabelFrame = Tkinter.LabelFrame(geneSetDiffExpFrame,  text=" 2 - Choose files list: ")
                filesListGeneSetDiffExpLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                self.createFileInformation("filesListGeneSetDiffExp", filesListGeneSetDiffExpLabelFrame)

                ## Dictionary ##
                dictionaryGeneSetDiffExpLabelFrame = Tkinter.LabelFrame(geneSetDiffExpFrame,  text=" 3 - Choose dictionary file: ")
                dictionaryGeneSetDiffExpLabelFrame.grid(row=2, padx=5, pady=5, sticky='WE')
                self.createFileInformation("dictionaryGeneSetDiffExp", dictionaryGeneSetDiffExpLabelFrame)

                ## Number of permutatins for frdr calculation ##
                permutationGeneSetDiffExpLabelFrame = Tkinter.LabelFrame(geneSetDiffExpFrame,  text=" 4 - Choose the number of permutations for the false discovery rate calculation: ")
                permutationGeneSetDiffExpLabelFrame.grid(row=3, padx=5, pady=5, sticky='WE')
                permutationGeneSetDiffExp = Tkinter.StringVar()
                permutationGeneSetDiffExp.set("1000")
                permutationGeneSetDiffExpEntry = Tkinter.Entry(permutationGeneSetDiffExpLabelFrame, validate="key", textvariable=permutationGeneSetDiffExp)
                permutationGeneSetDiffExpEntry['validatecommand'] = (self.register(self.validateIntegerKey), '%S', '%d')
                permutationGeneSetDiffExpEntry.grid(padx=5, pady=5)
                
                ## Headers ##
                headerLabelGeneSetDiffExpFrame = Tkinter.LabelFrame(geneSetDiffExpFrame,  text=" 5 - Choose the identifiers of the groups: ")
                headerLabelGeneSetDiffExpFrame.grid(row=4, padx=5, pady=5, sticky='WE')
                headerLabelGeneSetDiffExp[0] = Tkinter.Label(headerLabelGeneSetDiffExpFrame,text="Group 1: ")
                headerLabelGeneSetDiffExp[0].grid(row=0,column=0,padx=5, pady=5)
                headerEntryGeneSetDiffExp[0] = Tkinter.Entry(headerLabelGeneSetDiffExpFrame)
                headerEntryGeneSetDiffExp[0].grid(row=0,column=1,padx=5, pady=5)
                insertGeneSetDiffExpButton = Tkinter.Button(headerLabelGeneSetDiffExpFrame, text='Add', command=lambda: self.createField(headerLabelGeneSetDiffExpFrame,headerEntryGeneSetDiffExp,headerLabelGeneSetDiffExp,insertGeneSetDiffExpButton))
                insertGeneSetDiffExpButton.grid()

                ##    Program    ##
                Tkinter.Button(geneSetDiffExpFrame, text='Calculate', command=lambda: self.runGeneSetDiffExp(permutationEntry.get())).grid()


                #################################
                ##           Dendrogram        ##
                #################################

                ## Expression ##
                expressionDendrogramLabelFrame = Tkinter.LabelFrame(dendrogramFrame,  text=" 1 - Choose transcriptogram file: ")
                expressionDendrogramLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("expressionDendrogram", expressionDendrogramLabelFrame)

                ##    Program    ##
                Tkinter.Button(dendrogramFrame, text='Make dendrogram', command=lambda: self.runDendrogram()).grid()


                #################################
                ##          Translate          ##
                #################################

                ## Ordering ##
                orderingTranslateLabelFrame = Tkinter.LabelFrame(translateFrame,  text=" 1a - Choose ordering file... ")
                orderingTranslateLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                self.createFileInformation("orderingTranslate", orderingTranslateLabelFrame)

                ## Expression ##
                expressionTranslateLabelFrame = Tkinter.LabelFrame(translateFrame,  text=" 1b - or choose a transcriptogram file: ")
                expressionTranslateLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                self.createFileInformation("expressionTranslate", expressionTranslateLabelFrame)

                ## Dictionary ##
                expressionTranslateLabelFrame = Tkinter.LabelFrame(translateFrame,  text=" 2 - Choose a dictionary file: ")
                expressionTranslateLabelFrame.grid(row=2, padx=5, pady=5, sticky='WE')
                self.createFileInformation("dictionaryTranslate", expressionTranslateLabelFrame)

                ## Program ##
                Tkinter.Button(translateFrame, text='Translate', command=lambda: self.runTranslate()).grid()

                #################################
                ##           Figure            ##
                #################################
                
                signameFigureLabelFrame = Tkinter.LabelFrame(figureFrame,  text=" 1 - Choose y axis label: ")
                signameFigureLabelFrame.grid(row=0, padx=5, pady=5, sticky='WE')
                signameFigure=Tkinter.StringVar()
                signameFigure.set("Signal")
                signameFigureEntry = Tkinter.Entry(signameFigureLabelFrame, validate="key", textvariable=signameFigure)
                signameFigureEntry.grid(padx=5, pady=5)

                signalFigureLabelFrame = Tkinter.LabelFrame(figureFrame,  text=" 2 - Choose file to plot. ")
                signalFigureLabelFrame.grid(row=1, padx=5, pady=5, sticky='WE')
                self.createFileInformation("signalFigure", signalFigureLabelFrame)
                
                pvalFigureLabelFrame = Tkinter.LabelFrame(figureFrame,  text=" 3 - Choose file with pvalue (optional, only with relative averages). ")
                pvalFigureLabelFrame.grid(row=2, padx=5, pady=5, sticky='WE')
                self.createFileInformation("pvalFigure", pvalFigureLabelFrame)
                
                goFigureLabelFrame = Tkinter.LabelFrame(figureFrame,  text=" 4 - Choose file with GO Terms Enrichment (optional). ")
                goFigureLabelFrame.grid(row=3, padx=5, pady=5, sticky='WE')
                self.createFileInformation("goFigure", goFigureLabelFrame)
                
                keggFigureLabelFrame = Tkinter.LabelFrame(figureFrame,  text=" 5 - Choose file with KEGG Terms Enrichment (optional). ")
                keggFigureLabelFrame.grid(row=4, padx=5, pady=5, sticky='WE')
                self.createFileInformation("keggFigure", keggFigureLabelFrame)
                
                dictFigureLabelFrame = Tkinter.LabelFrame(figureFrame,  text=" 6 - Choose dictionary file (optional). ")
                dictFigureLabelFrame.grid(row=4, padx=5, pady=5, sticky='WE')
                self.createFileInformation("dictFigure", dictFigureLabelFrame)
                
                Tkinter.Button(figureFrame, text='Plot', command=lambda: self.runSimplePlot(signameFigure.get())).grid(row=5, padx=300, sticky='W')

                Tkinter.Button(figureFrame, text='Stop Plot Server', command=lambda: self.runStopPlot(signameFigure.get())).grid(row=5, padx=300, sticky='E')

                #Tkinter.Button(figureFrame, text='Plot', command=lambda: self.runSimplePlot(signameFigure.get())).grid()

                


        #################################
        ##  Calling external programs  ##
        #################################

        #################################
        def runOrderer(self, isothermalSteps, monteCarloSteps, coolingFactorText, alphaText, percentualWeigthText):
                if not self.validateFiles(["association"]):
                        return

                try: 
                        coolingFactor = float(coolingFactorText)
                except:
                        tkMessageBox.showwarning("Error","This cooling factor is not valid!")
                        return

                try: 
                        alpha = float(alphaText)
                except:
                        tkMessageBox.showwarning("Error","This alpha is not valid!")
                        return

                try: 
                        percentualWeigth = float(percentualWeigthText)
                except:
                        tkMessageBox.showwarning("Error","This percentual energy is not valid!")
                        return

                subprocess.call([binPath+'ordering1D', 'f='+filesNames["association"].encode(sys.getfilesystemencoding()), 'i='+isothermalSteps, 'm='+monteCarloSteps, 'c='+str(coolingFactor), 'a='+str(alpha), 'p='+str(percentualWeigth)])
                

        #################################
        def runProperties(self, radiusText):
                if not self.validateFiles(["orderingModularity","associationModularity"]):
                        return

                try: 
                        radius = float(radiusText)
                except:
                        tkMessageBox.showwarning("Error","This radius is not valid!")
                        return

                subprocess.call([binPath+'orderingProperties', 'o='+filesNames["orderingModularity"].encode(sys.getfilesystemencoding()), 'a='+filesNames["associationModularity"].encode(sys.getfilesystemencoding()), 'r='+str(radius)])


        #################################
        def runNetwork(self):
                if not self.validateFiles(["associationNetwork"]):
                        return

                subprocess.call([binPath+'networkProperties', 'a='+filesNames["associationNetwork"].encode(sys.getfilesystemencoding())])


        #################################
        def runMatrix(self):

                if not self.validateFiles(["orderingMatrix","associationMatrix"]):
                        return

                subprocess.call([binPath+'mountInteractionMatrix','o='+filesNames["orderingMatrix"].encode(sys.getfilesystemencoding()),'a='+filesNames["associationMatrix"].encode(sys.getfilesystemencoding())])

                # save gnuplot file
                print("Saving gnuplot script file ... ")
                try:
                        f = open('matrix_gnuplot.plt', 'w')
                        f.write("set tics scale 0\n"
                                                "# edit the next line to select the file format\n"
                                                "set terminal postscript\n"
                                                "set nokey\n"
                                                "set xlabel \"Gene ID\"\n"
                                                "set ylabel \"Gene ID\"\n"
                                                "set output \'matrix.ps\'\n"
                                                "# edit the next line to correspond to the ordering file name\n"
                                                "plot \'associationMatrix.dat\' w dots lt 0\n")
                        f.close()
                except (OSError, IOError) as e:
                        print("ERROR: Cannot save gnuplot script file !!!\n")
                print("Done !")


        #################################
        def runOccupation(self):

                if not self.validateFiles(["orderingOccupation","associationOccupation"]):
                        return

                subprocess.call([binPath+'occupationLevel','o='+filesNames["orderingOccupation"].encode(sys.getfilesystemencoding()),'a='+filesNames["associationOccupation"].encode(sys.getfilesystemencoding())])

                # save gnuplot file
                print("Saving gnuplot script file ... ")
                try:
                        f = open('occupation_gnuplot.plt', 'w')
                        f.write( "set size square\n"
                                                "set tics scale 0\n"
                                                "set terminal postscript\n"
                                                "set nokey\n"
                                                "set xlabel \"norm distance\"\n"
                                                "set logscale y\n"
                                                "set ylabel \"occupation level (log)\"\n"
                                                "set output \'occupation.ps\'\n"
                                                "plot \'occupation.dat\' using 1:2 with dots\n")
                        f.close()
                except (OSError, IOError) as e:
                        print("ERROR: Cannot save gnuplot script file !!!\n")
                print("Done !")

        #################################
        def runTranscriptogramer(self, radiusText):
                if not self.validateFiles(["ordering","dictionary","expression"]):
                        return

                try: 
                        radius = float(radiusText)
                except:
                        tkMessageBox.showwarning("Error","This radius is not valid!")
                        return

                subprocess.call([binPath+'transcriptogram', 'o='+filesNames["ordering"].encode(sys.getfilesystemencoding()), 'd='+filesNames["dictionary"].encode(sys.getfilesystemencoding()), 't='+filesNames["expression"].encode(sys.getfilesystemencoding()), 'r='+str(radius)])

                # save gnuplot file
                print("Saving gnuplot script file ... ")
                try:
                        f = open('transcriptogram_gnuplot.plt', 'w')
                        
                        f.write( "set term post eps enhanced color \"Arial\" 11\n"
						"set xlabel \"Gene ID\"\n"
						"# greek letter only work with eps format ! no png !\n"
						"set ylabel '{/Symbol t}'\n"
						"set output 'transcript.eps'\n"
						"# draw lines from columns 3 to column 6. To draw more lines,\n"
						"# for instance 20 lines, just change '[i=3:6]' to '[i=3:20]'.\n"
						"plot for [i=3:6] \"transcriptogram.dat\" using i with lines\n"
						"# the user might also want the following chart format. In this case,\n"
						"# it is suggested to sort the 2nd column of the .dat file. In Linux, ese command 'sort'.\n"
						"#plot for [i=3:6] \"transcriptogram.dat\" using i with filledcurves closed\n")
                        f.close()
                except (OSError, IOError) as e:
                        print("ERROR: Cannot save gnuplot script file !!!\n")
                print("Done !")

        #################################
        def runTermEnrichment(self, radiusText):
                if not self.validateFiles(["orderingEnrichment","dictionaryEnrichment","filesList"]):
                        return        
                try: 
                        radius = float(radiusText)
                except:
                        tkMessageBox.showwarning("Error","This radius is not valid!")
                        return

                subprocess.call([binPath+'termEnrichment', 'o='+filesNames["orderingEnrichment"].encode(sys.getfilesystemencoding()), 'd='+filesNames["dictionaryEnrichment"].encode(sys.getfilesystemencoding()), 'l='+filesNames["filesList"].encode(sys.getfilesystemencoding()), 'r='+str(radius)])


        #################################
        def runStatistics(self, permutation):
                if not self.validateFiles(["expressionStatistics"]):
                        return        

                headers = ','.join(headerEntryStatistics[entry].get() for entry in headerEntryStatistics)
                if not headers:
                        tkMessageBox.showwarning("Error","You must type at least one header identifier.")
                        return

                subprocess.call([binPath+'statistics', 't='+filesNames["expressionStatistics"].encode(sys.getfilesystemencoding()), 'p='+permutation, 'h='+headers])


        #################################
        def runGeneSetDiffExp(self, permutation):
                if not self.validateFiles(["expressionGeneSetDiffExp","dictionaryGeneSetDiffExp","filesListGeneSetDiffExp"]):
                        return        

                headers = ','.join(headerEntryGeneSetDiffExp[entry].get() for entry in headerEntryGeneSetDiffExp)
                if not headers:
                        tkMessageBox.showwarning("Error","You must type at least one header identifier.")
                        return

                subprocess.call([binPath+'geneSetDiffExp', 't='+filesNames["expressionGeneSetDiffExp"].encode(sys.getfilesystemencoding()), 'd='+filesNames["dictionaryGeneSetDiffExp"].encode(sys.getfilesystemencoding()), 'l='+filesNames["filesListGeneSetDiffExp"].encode(sys.getfilesystemencoding()), 'p='+permutation, 'h='+headers])


        #################################
        def runDendrogram(self):
                if not self.validateFiles(["expressionDendrogram"]):
                        return        

                subprocess.call([binPath+'dendrogram', 't='+filesNames["expressionDendrogram"].encode(sys.getfilesystemencoding())])

        #################################
        def runTranslate(self):
                if not self.validateFiles(["dictionaryTranslate"]):
                        return        

                if not self.validateFilesOptional(["orderingTranslate"]):
                        if not self.validateFilesOptional(["expressionTranslate"]):
                                tkMessageBox.showwarning("Error","You need to select the transcriptogram or the ordering file to translate!")
                                return
                        subprocess.call([binPath+'translate', 'd='+filesNames["dictionaryTranslate"].encode(sys.getfilesystemencoding()), 't='+filesNames["expressionTranslate"].encode(sys.getfilesystemencoding())])
                        return

                if not self.validateFilesOptional(["expressionTranslate"]):
                        if not self.validateFilesOptional(["orderingTranslate"]):
                                tkMessageBox.showwarning("Error","You need to select the transcriptogram or the ordering file to translate!")
                                return
                        subprocess.call(['binPath+/translate', 'd='+filesNames["dictionaryTranslate"].encode(sys.getfilesystemencoding()), 'o='+filesNames["orderingTranslate"].encode(sys.getfilesystemencoding())])
                        return

                subprocess.call(['binPath+/translate', 'd='+filesNames["dictionaryTranslate"].encode(sys.getfilesystemencoding()), 't='+filesNames["expressionTranslate"].encode(sys.getfilesystemencoding()), 'o='+filesNames["orderingTranslate"].encode(sys.getfilesystemencoding())])
        
        #################################
        def runPlot(self,signameFigure):
            command = 'bokeh serve --show ' + binPath+'plot.py --args ' + signameFigure +' ' + filesNames["signalFigure"].encode(sys.getfilesystemencoding()) 
            if self.validateFilesOptional(["pvalFigure"]):
                command += ' ' + filesNames["pvalFigure"].encode(sys.getfilesystemencoding())
                if self.validateFilesOptional(["goFigure"]):
                    command += ' ' + filesNames["goFigure"].encode(sys.getfilesystemencoding())
                    if self.validateFilesOptional(["keggFigure"]):
                        command += ' ' + filesNames["keggFigure"].encode(sys.getfilesystemencoding())
                            
#                        return        

#                subprocess.call(['bokeh serve --show ' + binPath + 'plot.py'], shell=True)
            subprocess.call([command], shell=True)

#################################
        def runSimplePlot(self,signameFigure):
        
            self.is_plotting = True
            command = 'python ' + binPath+'simple_plot.py \'' + signameFigure +'\' \'' + filesNames["signalFigure"].encode(sys.getfilesystemencoding()) + '\''
            if self.validateFilesOptional(["pvalFigure"]):
                command += ' \'' + filesNames["pvalFigure"].encode(sys.getfilesystemencoding()) + '\''
            else:
                command += ' NONE'
            if self.validateFilesOptional(["goFigure"]):
                command += ' \'' + filesNames["goFigure"].encode(sys.getfilesystemencoding()) + '\''
            else:
                command += ' NONE'
            if self.validateFilesOptional(["keggFigure"]):
                command += ' \'' + filesNames["keggFigure"].encode(sys.getfilesystemencoding()) + '\''
            else:
                command += ' NONE'
            if self.validateFilesOptional(["dictFigure"]):
                command += ' \'' + filesNames["dictFigure"].encode(sys.getfilesystemencoding()) + '\''
            else:
                command += ' NONE'
                            
#                        return        

#                subprocess.call(['bokeh serve --show ' + binPath + 'plot.py'], shell=True)
            self.plot_th = subprocess.Popen([command], shell=True)
            #print(command)
        
    
#################################
        def runSSimplePlot(self,signameFigure):
        
            self.is_plotting = True

    
            self.plot_th = threading.Thread(target=self.runThreadPlot(signameFigure))
            self.plot_th.start()
    
#################################
        def runStopPlot(self,signameFigure):
            
            self.plot_th.kill()
#             try:
#                 self.plot_th.terminate()
#             except (AttributeError, RuntimeError):
#                 pass
            
        #################################
        ##                             ##
        #################################
        def askOpenFile(self, fileIdentifier):
                filesNames[fileIdentifier] = tkFileDialog.askopenfilename()
                if not filesNames[fileIdentifier]:
                        filesNames[fileIdentifier] = "NONE"
                        return
                labelFileName[fileIdentifier]["text"] = filesNames[fileIdentifier]
                labelFileLines[fileIdentifier]["text"] = sum(1 for line in open(filesNames[fileIdentifier]))


        #################################
        ##                             ##
        #################################
        def createFileInformation(self, fileIdentifier, labelFrame):
                Tkinter.Button(labelFrame, text='Choose file', command=lambda: self.askOpenFile(fileIdentifier)).grid(row=0, column=0, rowspan=2, padx=5, pady=5, sticky='WE')

                Tkinter.Label(labelFrame, text="File:").grid(row=0, column=1, padx=5, pady=5, sticky='WE')
                labelFileName[fileIdentifier] = Tkinter.Label(labelFrame, text="NONE")
                labelFileName[fileIdentifier].grid(row=0, column=2, padx=10, sticky='WE')

                Tkinter.Label(labelFrame, text="Lines:").grid(row=1, column=1, padx=5, pady=5, sticky='WE')
                labelFileLines[fileIdentifier] = Tkinter.Label(labelFrame, text="---")
                labelFileLines[fileIdentifier].grid(row=1, column=2, padx=10, sticky='WE')


        #################################
        ##                             ##
        #################################
        def createField(self, frame, entry, label, button):
                nFields = len(entry)
                label[nFields] = Tkinter.Label(frame,text="Group "+str(nFields+1)+": ")
                label[nFields].grid(row=nFields,column=0,padx=5, pady=5)
                entry[nFields] = Tkinter.Entry(frame)
                entry[nFields].grid(row=nFields, column=1, padx=5, pady=5)
                button.grid(row=nFields+1)


        #################################
        ## Validation ##
        #################################
        def validateFloatKey(self, key, command):
                if (command == "1"):
                        return (key in ["0","1","2","3","4","5","6","7","8","9","."])
                else:
                        return True

        def validateIntegerKey(self, key, command):
                if (command == "1"):
                        return (key in ["0","1","2","3","4","5","6","7","8","9"])
                else:
                        return True

        def validateFiles(self, listOfFilesToValidate):
                for fileIdentifier in listOfFilesToValidate:
                        if filesNames[fileIdentifier] == "NONE":
                                tkMessageBox.showwarning("Error","You need to select all the files from this tab!")
                                return False
                return True

        def validateFilesOptional(self, listOfFilesToValidate):
                for fileIdentifier in listOfFilesToValidate:
                        if filesNames[fileIdentifier] == "NONE":
                                return False
                return True
                

root = Tkinter.Tk(className='Transcriptogramer')
platform = platform.system()

if(platform=='Linux'):
	try:
		img = Tkinter.PhotoImage(file='/usr/local/share/transcriptogram/logo.gif')
	except:
		img = Tkinter.PhotoImage(file='~/bin/logo.gif')
else:
	img = Tkinter.PhotoImage(file=(binPath + '/logo.gif'))

                

root.tk.call('wm', 'iconphoto', root._w, img)

Application(root).grid()

root.title("Transcriptogramer")
root.mainloop()
