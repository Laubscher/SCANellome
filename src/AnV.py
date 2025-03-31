#!/usr/bin/env python3

# Import modules - dependencies
import os                                 # system
import shutil
import time
from os import listdir
import mappy as mp                        # mapper
from statistics import median
import tkinter as tk                      # graphic interface
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter.messagebox import showinfo
import Dicodb                             # module that contain the database in a python dico
import base64
import img
from ttkthemes import ThemedTk
from tkinter import Menu
import webbrowser
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import glob
import pysam
from PIL import ImageTk, Image
from tkinter import Frame
from tkinter import BOTH
from tkinter import Scrollbar
from tkinter import VERTICAL
from tkinter import RIGHT
from tkinter import Y
from tkinter import Canvas
#from tkhtmlview import HTMLLabel

class ScrollableFrame:
    def __init__ (self,master,width,height,mousescroll=0):
        self.mousescroll = mousescroll
        self.master = master
        self.height = height
        self.width = width
        self.main_frame = Frame(self.master)
        self.main_frame.pack(fill=BOTH,expand=1)

        self.scrollbar = Scrollbar(self.main_frame, orient=VERTICAL)
        self.scrollbar.pack(side=RIGHT,fill=Y)

        self.canvas = Canvas(self.main_frame,yscrollcommand=self.scrollbar.set)
        self.canvas.pack(expand=True,fill=BOTH)

        self.scrollbar.config(command=self.canvas.yview)

        self.canvas.bind('<Configure>', lambda e: self.canvas.configure(scrollregion = self.canvas.bbox("all")))

        self.frame = Frame(self.canvas,width=self.width,height=self.height)
        self.frame.pack(expand=True,fill=BOTH)
        self.canvas.create_window((0,0), window=self.frame, anchor="nw")

        self.frame.bind("<Enter>", self.entered)
        self.frame.bind("<Leave>", self.left)

    def _on_mouse_wheel(self,event):
        self.canvas.yview_scroll(-1 * int((event.delta / 120)), "units")

    def entered(self,event):
        if self.mousescroll:
            self.canvas.bind_all("<MouseWheel>", self._on_mouse_wheel)
        
    def left(self,event):
        if self.mousescroll:
            self.canvas.unbind_all("<MouseWheel>")


# The software window
main1 = ThemedTk(theme="ubuntu", background=True, className="SCANellome v. 2.0.0")
main1.title('                                                                    SCANellome                                                                         v. 2.0.0')
main1.geometry("4000x1000")
obj = ScrollableFrame(main1,height=11300,width=4000 )
main = obj.frame


global pb

def start():   # start is a function that check or make the file structure for the whole application
  #Variables
  global fastq1Path
  global fastq2Path
  global fastq1
  global fastq2
  global fq1
  global fq2
  global sampleList
  global yAdd
  global pathLastDir
  global pairedLabel
  global text1Label
  global text2Label
  global minion
  global pathData
  global projectList

  pathData = os.path.expanduser("~/.AnV")                               # check if .AnV otherwise mk it
  if not os.path.exists(pathData):
    os.mkdir(pathData)

  #####################################
  #               img                 #
  #####################################

  if not os.path.exists(pathData + "/IMG"):
    os.mkdir(pathData + "/IMG")
    image = open(pathData + "/IMG/blue.png", 'wb')
    image.write(base64.b64decode((img.img)))
    image.close()

  try:
      icon = tk.PhotoImage(file = pathData + "/IMG/blue.png")
      main1.iconphoto(False, icon)
  except:
      pass

  #####################################
  #               tree                #
  #####################################

  if not os.path.exists(pathData + "/IMG/TREE/"):
    os.mkdir(pathData + "/IMG/TREE")
    image = open(pathData + "/IMG/TREE/GENcol.png", 'wb')
    image.write(base64.b64decode((img.gen)))
    image.close()

    image = open(pathData + "/IMG/TREE/ALPHA.png", 'wb')
    image.write(base64.b64decode((img.alpha)))
    image.close()

    image = open(pathData + "/IMG/TREE/BETA.png", 'wb')
    image.write(base64.b64decode((img.beta)))
    image.close()

    image = open(pathData + "/IMG/TREE/GAMMA.png", 'wb')
    image.write(base64.b64decode((img.gamma)))
    image.close()





  #####################################
  #               log                 #
  #####################################

  if not os.path.exists(pathData + "/LOG"):
    os.mkdir(pathData + "/LOG")
    log = open(pathData + "/LOG/log.txt", "w")
    firstSessionTime = time.time()
    log.write("Set up first session: " + str(firstSessionTime) + "\n")
    log.close()

  log = open(pathData + "/LOG/log.txt", "a")
  sessionTime = time.time()

  log.write("\nSession start: " + str(sessionTime) +"\n")

  #####################################

  #####################################
  #             USERDATA              #
  #####################################

  if not os.path.exists(pathData + "/USERDATA"):
    os.mkdir(pathData + "/USERDATA")
    os.mkdir(pathData + "/USERDATA/default")
    log.write("USERDATA set up" + "\n")

  pathLastDir = os.path.expanduser("~")
  projectList = listdir(pathData + "/USERDATA/")  # list of all "project" from past session - for each project one directory in USERDATA/
  #sampleUniq = set(listdir(pathData + "/USERDATA/"))

  log.write("List of project found in USERDATA/: " + str(projectList) + "\n")

  #####################################
  #               db                  #
  #####################################
  #make the fasta db from module Dicodb
  #####################################

  if not os.path.exists(pathData + "/DATABASE"):
    os.mkdir(pathData + "/DATABASE")
  if not os.path.exists(pathData + "/DATABASE/2024.2.fasta"):
    fasta=open(pathData + "/DATABASE/2024.2.fasta", "w")
    for entry in Dicodb.db:
      if entry != "none":
        fasta.write(">"+entry+"\n")
        fasta.write(Dicodb.db[entry][5]+"\n")
    fasta.close()
    log.write("DATABASE set up" + "\n")
  #####################################

  log.close()

  sampleList = []

def run():        #start the run -> mapping + analyse of all sample in the sample list, throws an error if mapping fail
    global pb
    global main
    log = open(pathData + "/LOG/log.txt", "a")
    log.write("Start mapping samples ... " + "\n")

    for s in sampleList:   #sample = [sampleName, fastq1Path, fastq2Path, yAdd + 20, minion.get()]
        print("")
        log.write("Mapping sample: " + str(s[0]) + "\n")
        pb = ttk.Progressbar(main, orient='horizontal', mode='determinate', length=280)  # Progress bar
        pb.place(x=350, y=s[3])
        main.update()
        pb['value'] += 1
        main.update()
        try:
          curr = time.time()
          if s[4] == 1:
              type="Nanopore"
              print("Nanopore")
          elif s[2] != "<empty>":
              type = "Paired"
              print("Paired")
          else:
              type="Single"
              print("single")
          mapping(s[1], pathData + "/DATABASE/2024.2.fasta", s[0], type, s[2])

        except:
          log.write("Error during mapping.. sample:" + str(s[0]) + "\n files may be corrupted!" + "\n")
          tk.messagebox.showinfo("Error..", "Error during mapping.. sample:" + str(s[0]) + "\n files may be corrupted!")
          sampleList.remove(s)
          os.rmdir(pathData + "/USERDATA/" + projectSelected + "/" + s[0])
          pb['value'] = 0

        else:
          pb['value'] += 1
          main.update()
          log.write("Mapping time sample: " + str(s[0]) + " " + str(time.time() - curr) + "\n")

    save_button = ttk.Button(
        main,
        text='Save CSV',
        command=file_save
    )
    save_button.place(x=150, y=195)

    save2_button = ttk.Button(
        main,
        text='Save FASTA',
        command=fasta_save
    )
    if makeConsensus.get() == 1:
      save2_button.place(x=290, y=195)
    reset_button.destroy()
    tk.messagebox.showinfo("Analysis completed", "Analysis completed" + "\n")

    log.close()

def map1Fastq(pathToFastq, genome_ref_covered, a):
      print("map1")
      print(pathToFastq)
      generator = mp.fastx_read(pathToFastq)
      for name, seq, qual in generator:
        for hit in a.map(seq):
            mgRef = hit.ctg.split(",")[0]                  #if not in dico make a new entry
            if mgRef in genome_ref_covered:
                pass
            else:
                genome_ref_covered[mgRef] = [[], 0,
                                             hit.ctg_len]  # list [0] -> depth for each pos (incrementation), [2] lenref, [1] count of reads (incrementation)
                for k in range(0, hit.ctg_len):
                    genome_ref_covered[mgRef][0].append(0) # populate depth with 0

            genome_ref_covered[mgRef][1] += 1
            for i in range(hit.r_st, hit.r_en):             # r_st  = ref start match, r_en -> end
                genome_ref_covered[mgRef][0][i] += 1

def map2sam(pathToFastq, genome_ref_covered, a, nameS):
    print("map4consensus")
    print(pathToFastq)
    generator = mp.fastx_read(pathToFastq)
    fichierSAM1 = open(pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/" + nameS + "1.sam", "w")
    fichierfastq2 = open(pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/" + nameS + "_2.fastq", "w")

    fichierHEADER = open(pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/" + nameS + ".sam", "w")
    #fichierHEADER.write("@HD	VN:1.0	SO:coordinate\n")
    reflist=[]

    phred = "B"*999999
    for name, seq, qual in generator:
        for hit in a.map(seq):
            if hit.ctg not in reflist:
              fichierHEADER.write("@SQ	SN:" + hit.ctg + "	LN:" + str(hit.ctg_len) + "\n")
              reflist.append(hit.ctg)
            if hit.strand == 1:    # orientation sequence +1

              fichierSAM1.write(name + "\t0\t" + hit.ctg + "\t" + str(hit.r_st+1) + "\t" + str(hit.mapq)+"\t" + str(hit.cigar_str) + "\t" + "*" + "\t" + "0" + "\t" + str( hit.blen) + "\t" + seq[hit.q_st:hit.q_en] + "\t" + phred[hit.q_st:hit.q_en] + "\n")
            if hit.strand == -1:

              # complement strand
              seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
              seq = seq.upper()
              seq = seq[::-1]
    fichierSAM1.close()
    fichierfastq2.close()
    fichierHEADER.close()

    fichierSAM = open(pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/" + nameS + ".sam", "a")  # open the file in append mode

    path2FastqTemp=pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/" + nameS + "_2.fastq"
    # reprocess for -1 misclassified
    generator = mp.fastx_read(path2FastqTemp)
    for name, seq, qual in generator:
        for hit in a.map(seq):
            print(hit.strand)
            if hit.strand == 1:

                fichierSAM.write( name + "\t0\t" + hit.ctg + "\t" + str(hit.r_st + 1) + "\t" + str(hit.mapq) + "\t" + str(hit.cigar_str) + "\t" + "*" + "\t" + "0" + "\t" + str(hit.blen) + "\t" + seq[hit.q_st:hit.q_en] + "\t" + phred[hit.q_st:hit.q_en] + "\n")


    fichierSAM.close()


    fichierSAM1=open(pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/" + nameS + "1.sam", "r")

    #header + append sam
    fichierSAM=open(pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/" + nameS + ".sam", "a")    # open the file in append mode

    for lane in fichierSAM1:
        fichierSAM.write(lane)
    fichierSAM1.close()
    fichierSAM.close()

    #ham + sam -> sam
    #variables for sam name

    PATH=pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/"

    SAM=PATH + nameS + ".sam"

    BAM=PATH + nameS + ".bam"

    FASTA = PATH + nameS + ".fasta"
    pysam.sort("-o", BAM, SAM)
    pysam.consensus("-o", FASTA, BAM)

    SAM1 =PATH + nameS + "1.sam"
    os.remove(SAM1)
    os.remove(BAM)
    os.remove(SAM)
    os.remove(path2FastqTemp)

def mapping(pathToFastq, db, nameS, type, pathToFastq2):
    global pb
    genome_ref_covered = dict()
    genome_ref_class = Dicodb.db                                    # key acc number; value list of info from the fasta header
    if type == "Nanopore":
        preset = "map-ont"     # for Oxford Nanopore read mapping
        a = mp.Aligner(db, preset=preset)
    elif type == "single":
        preset = "sr"           # single-end short reads
        a = mp.Aligner(db, preset=preset)
    else:
        a = mp.Aligner(db)
    print(a)
    pb['value'] += 2
    main.update()

    print(type)

    if type == "Paired":
      print("map")

      #fastq1
      map1Fastq(pathToFastq, genome_ref_covered, a)

      pb['value'] += 40
      main.update()

      #fastq2
      map1Fastq(pathToFastq2, genome_ref_covered, a)

    else:
      map1Fastq(pathToFastq, genome_ref_covered, a)

      pb['value'] += 40
      main.update()

    if makeConsensus.get() == 1:
    #slow mode get consensus
      map2sam(pathToFastq, genome_ref_covered, a, nameS)

    pb['value'] += 58
    main.update()

    #for each genome with match -> compute metrics

    for l in genome_ref_covered:
        genome_ref_covered[l].append(genome_ref_covered[l][2] - genome_ref_covered[l][0].count(0))  # 3 coverage (nt)
        genome_ref_covered[l].append(
            round((genome_ref_covered[l][2] - genome_ref_covered[l][0].count(0)) / genome_ref_covered[l][2] * 100,
                  2))  # 4 % cov
        genome_ref_covered[l].append(median(genome_ref_covered[l][0]))  # 5 cov depth

    pb['value'] += 2
    main.update()

    # *# ############################# #*#
    # *#                               #*#
    # *#   Get best hit at sp level    #*#
    # *#                               #*#
    # *# ############################# #*#

    listSpecies = set()   # set of sp name
    dictSpecies = dict()  # key sp ; value list of acc number.
    # populate the set and dict
    for l in genome_ref_covered:
        if genome_ref_covered[l][4] >= 50:

            listSpecies.add(genome_ref_class[l][2])
            if genome_ref_class[l][2] in dictSpecies:
                pass
            else:
                dictSpecies[genome_ref_class[l][2]] = set()
            dictSpecies[genome_ref_class[l][2]].add(l)

    pb['value'] += 2
    main.update()

    resultSp = set()
    for sp in listSpecies:  # for each sp keep the one with the best coverage
        bestCov = 49
        for acc in dictSpecies[sp]:
            if genome_ref_covered[acc][4] > bestCov:
                keep = acc
                bestCov = genome_ref_covered[acc][4]
        resultSp.add(keep)

    pb['value'] += 2
    main.update()

    fichierCSV = open(pathData + "/USERDATA/" + projectSelected + "/" + nameS + "/species.csv", "w")
    #Sample Name, ACC. NUMBER, Reads, ref_len, cov, %cov, depth (median), GENUS, GROUP, SPECIES, GENOTYPE, HOST
    for accResultSp in resultSp:
        fichierCSV.write(nameS + ", " + accResultSp + ", " + str(genome_ref_covered[accResultSp][1]) + ", " + str(genome_ref_covered[accResultSp][2]) + ", " + str(
            genome_ref_covered[accResultSp][3]) + ", " + str(genome_ref_covered[accResultSp][4]) + ", " + str(
            genome_ref_covered[accResultSp][5]) + ", " + genome_ref_class[accResultSp][0].split("=")[1]+ ", " + genome_ref_class[accResultSp][1].split("=")[1]+ ", " + genome_ref_class[accResultSp][2].split("=")[1] + ", " + genome_ref_class[accResultSp][3].split("=")[1] + ", " + genome_ref_class[accResultSp][4].split("=")[1] + "\n")
    fichierCSV.close()
    pb['value'] += 4
    main.update()

######################### batch ##################################################
def select_file_batch():
    global pathLastDir

    filetypes = (
        ('Gzipped fastq', '*.fastq.gz'),
        ('Fastq files', '*.fastq'),
        ('All files', '*.*')
    )
    filenames = fd.askopenfilenames(
        title='select fastq',
        initialdir=pathLastDir,
        filetypes=filetypes)

    for path in filenames:
        fastq1Path=str(path)
        if fastq1Path.rstrip().split(".")[-1] == "fastq":

            pass
        elif fastq1Path.rstrip().split(".fastq.")[-1] == "gz":
            pass

        else:
            pass
        add_sample_batch(fastq1Path)

    pathL = fastq1Path.split("/")           #for keep the path of the directory
    pathL.pop()
    pathLastDir = "/".join(pathL)

def add_sample_batch(fastq1Path):
    global yAdd
    global sampleList
    global sampleUniq

    isSample = True
    if illuminaPE.get() == 1:
        if fastq1Path.split("_R1_001")[-1] != ".fastq":
          if fastq1Path.split("_R1_001.fastq")[-1] != ".gz":
            isSample=False
          else:
            sampleName = fastq1Path.split("/")[-1].split("_R1_001.fastq")[0]
        else:
          sampleName = fastq1Path.split("/")[-1].split("_R1_001.fastq")[0]
        if fastq1Path.split("_R1_001")[-1] == ".fastq":
          fastq2Path = fastq1Path.split("_R1_001.fastq")[0] + "_R2_001.fastq"
        if fastq1Path.split("_R1_001.fastq")[-1] == ".gz":
            fastq2Path = fastq1Path.split("_R1_001.fastq")[0] + "_R2_001.fastq.gz"

    else :
      fastq2Path = "<empty>"
      if fastq1Path.split(".")[-1] != "fastq":
          if fastq1Path.split(".")[-1] != "gz":
              isSample = False
          else:
              sampleName = fastq1Path.split("/")[-1].split(".fastq")[0]

      else:
          sampleName = fastq1Path.split("/")[-1].split(".fastq")[0]


    if isSample:
        suffix = fastq1Path.split(".")[-1]
    # check if fastq

    # check if sample name is uniq
        if sampleName in sampleUniq:
          print("Sample name not uniq!!")
          if tk.messagebox.askokcancel("Can't add your sample..",'\n"' + sampleName + '"\n\nThe sample name already exist in this project.\n\n'+'Do you want to re analyze it?') == True:
              shutil.rmtree(pathData + "/USERDATA/" + projectSelected + "/" + sampleName)
              os.mkdir(pathData + "/USERDATA/" + projectSelected + "/" + sampleName)
              sample = [sampleName, fastq1Path, fastq2Path, yAdd + 20, minion.get()]
              sampleList.append(sample)
              labelListSample = ttk.Label(main, text=sampleName)
              yAdd += 20
              labelListSample.place(x=10, y=yAdd)
          else:
              print(False)



        elif suffix != "fastq":
          print("Sample suffix is not fastq!!")
          if suffix == "gz":
              os.mkdir(pathData + "/USERDATA/" + projectSelected + "/" + sampleName)
              sample = [sampleName, fastq1Path, fastq2Path, yAdd + 20, minion.get()]
              sampleList.append(sample)
              labelListSample = ttk.Label(main, text=sampleName)
              yAdd += 20
              labelListSample.place(x=10, y=yAdd)
              sampleUniq.add(sampleName)

          else:

            tk.messagebox.showinfo("Sorry can't add your sample..", '\n"' + suffix + '"\n\nThe file type is not fastq!')

        else:
          os.mkdir(pathData + "/USERDATA/" + projectSelected + "/" + sampleName)
          sample = [sampleName, fastq1Path, fastq2Path, yAdd + 20, minion.get()]
          sampleList.append(sample)
          labelListSample = ttk.Label(main, text=sampleName)
          yAdd += 20
          labelListSample.place(x=10, y=yAdd)
          sampleUniq.add(sampleName)

    print(sampleList)

def file_save():
    nameNewCsv = fd.asksaveasfile(mode='w', defaultextension=".csv")
    text2save="Sample Name, ACC. NUMBER, Reads, ref_len, cov, %cov, depth (median), GENUS, GROUP, SPECIES, GENOTYPE, HOST\n"
    for s in sampleList:
        try:
            fichier = open(pathData + "/USERDATA/" + projectSelected + "/" + s[0] + "/species.csv", "r")    # if error during mapping file will not existe
            for lane in fichier.read():
              text2save += lane
            fichier.close()
        except:
         pass
    nameNewCsv.write(text2save)
    nameNewCsv.close

def html_save():
    nameNewHtml = fd.asksaveasfile(mode='w', defaultextension=".html")
    text2save=""
    fichier = open(pathData + "/USERDATA/" + projectSelected + "/" + "/file.html", "r")   
    for lane in fichier.read():
              text2save += lane
    fichier.close()

    nameNewHtml.write(text2save)
    nameNewHtml.close()





def fasta_save():

    nameNewFasta = fd.asksaveasfile(mode='w', defaultextension=".FASTA")

    text2save=""
    for s in sampleList:
        print(s)
        fastaDico={}
        try:
            fasta=open(pathData + "/USERDATA/" + projectSelected + "/" + s[0] + "/" + s[0]+".fasta", "r")
            print(pathData + "/USERDATA/" + projectSelected + "/" + s[0] + "/" + s[0]+".fasta")
            for lane in fasta:
              if lane[0]==">":
                  ac=lane.rstrip().split(">")[1]
                  fastaDico[ac] = ""
              else:
                  fastaDico[ac]+=lane
            fasta.close()
        except:
            pass
        headerDico={}
        try:
            fichier = open(pathData + "/USERDATA/" + projectSelected + "/" + s[0] + "/species.csv", "r")    # if error during mapping file will not exist
            for lane in fichier:
              try:
                ac=lane.split(",")[1].split(" ")[1]

                headerDico[ac]=lane.split(",")[9].split(" ")[1]+"_"+lane.split(",")[0]+"_"+lane.split(",")[7].split(" ")[1]
              except:
                pass
            fichier.close()

        except:
           pass

        for ac in headerDico:

            try:
              text2save += ">" + str(headerDico[ac]) +"\n" + str(fastaDico[ac]) +"\n"
            except:
              pass
    nameNewFasta.write(text2save)




    nameNewFasta.close

def deleteA():
  if tk.messagebox.askokcancel("Warning..",'\n\nAre you sure you want to delete all data?\n\n') == True:
    global sampleUniq
    sampleUniq = []
    shutil.rmtree(pathData)
    start()
    for widgets in main.winfo_children():
        widgets.destroy()
    topButton()

def deleteP():
  if tk.messagebox.askokcancel("Warning..", '\n\nAre you sure you want to this project?\n\n') == True:
    global sampleUniq
    global projectSelected
    sampleUniq = []
    try:
        projectSelected = cb1.get()
    except:
        pass
    if projectSelected == "":
        print("no project to delete")
    else:
        print("delete" + projectSelected)
        shutil.rmtree(pathData + "/USERDATA/" + projectSelected)
    start()
    for widgets in main.winfo_children():
        widgets.destroy()
    topButton()

def resetA():
    global sampleList
    global sampleUniq
    for sample in sampleList:
        os.rmdir(pathData + "/USERDATA/" + projectSelected + "/" + str(sample[0]))
    sampleUniq = []
    start()
    for widgets in main.winfo_children():
        widgets.destroy()
    topButton()

def testCheck():
    for sample in dicoSampleVar:
        print(dicoSampleVar[sample].get())

def dataA():
    global projectList
    global dicoSampleVar
    global projectSelected

    try:
      projectSelected = cb1.get()
    except:
      pass

    start()
    for widgets in main.winfo_children():
        widgets.destroy()
    topButton()

    # test if project exist otherwise mk directory
    if projectSelected == "" :
        projectSelected="default"
    if not os.path.exists(pathData + "/USERDATA/" + projectSelected):
        tk.messagebox.showinfo("Error..", "Error " + projectSelected + "do not exist!")

    text0Label = ttk.Label(main, text="Selected project: ", foreground="gray")
    text05Label = ttk.Label(main, text=projectSelected, foreground="darkorange1")
    text05Label.place(x=510, y=5)  #
    text0Label.place(x=350, y=5)  #
    # file name in a list, then read once all the file to make a set of virus name in all the project
    # then reread all the file to make a dictionary where all key are sampleID and values are
    # an other dictionary with keys as virusNames and values as nb of reads

    sampleInProject = listdir(pathData + "/USERDATA/" + projectSelected)  # list of all sample in the project

    dicoSampleCheck = {}
    dicoSampleVar = {}

    for sample in sampleInProject:
      if sample!="file.html":
        dicoSampleVar[sample] = tk.IntVar(main, 0)               #variable for checkbox one for each sample in a dico
        dicoSampleVar[sample].set(1)                             #we don't know the key in advance so loop in a dico
        dicoSampleCheck[sample] = ttk.Checkbutton(               #same for the check button it self
        main,
        text=str(sample),
        variable=dicoSampleVar[sample],
        command=testCheck
        )
        print(str(sample))

    y1=90
    for i in dicoSampleCheck:
        y1+=20
        dicoSampleCheck[i].place(x=30, y=y1)


    grid_button = ttk.Button(
        main,
        text='Display selected samples',
        command=grid
    )
    grid_button.place(x=30, y=60)

def grid():
    global projectList
    global dicoSampleVar

    fileList0 = glob.glob(pathData + "/USERDATA/" + projectSelected +"/*/species.csv")

    for sample in dicoSampleVar:
        if dicoSampleVar[sample].get() == 0:
            fileList0.remove(pathData + "/USERDATA/" + projectSelected + "/" + sample + "/species.csv")

    sampleDico = {}

    genusNameSet = set()

    sampleNameSet = set()

    fileList=[]

    for file in fileList0:
        print(file)
        f = open(file, "r")
        isRSLTS=False              #if any output for this sample
        for line in f:
              genusNameSet.add(line.split(",")[7])
              sampleNameSet.add(line.split(",")[0])
              isRSLTS=line.split(",")[0]
        print(isRSLTS)
        if isRSLTS:
          print(isRSLTS, file)
          fileList.append(file)
        f.close()

    # empty files are excluded from the list of samples We can include them must retrieve the name and put none everywhere


    if fileList==[]:
        tk.messagebox.showinfo("No detection..", "None Anellovirus detected in any selected sample.")
        print("No anellovirus detected")
    else:
      sampleNameList = list(sampleNameSet)
      genusNameList = sorted(list(genusNameSet))

    # virusNameSet.remove('Virus')

      fig = make_subplots(len(genusNameList), 1)  # make one subplot for each genus
      n = 0  # counter n-eme subplot

      clb = {}  # dictionary for colorbar gestion
      #get nb virus by genus in
      listLenVirusByGenus=[]
      for genus in genusNameList:
        n += 1
        virusNameSet = set()
        for file in fileList:
            f = open(file, "r")
            for line in f:
                if genus == line.split(",")[7]:
                    virusNameSet.add(line.split(",")[9])
            f.close()
        virusNameList = list(virusNameSet)
        listLenVirusByGenus.append(len(virusNameList))  #to get the max number of virus in a grid for later set the height of the graph
        for file in fileList:
            virusNameDico = {}

            # getsample name should be unique in one file
            sampleNameSetInFile = set()
            f = open(file, "r")
            for line in f:
                sampleNameSetInFile.add(line.split(",")[0])
            f.close()
            sampleName = list(sampleNameSetInFile)[0]
            print(sampleName)

            # populate dictionary entry with each virus found in the project and value 0
            for virusName in virusNameSet:
                virusNameDico[virusName] = 0

            sampleDico[sampleName] = virusNameDico
            f = open(file, "r")
            for line in f:
                if line.split(",")[9] in sampleDico[sampleName]:
                    sampleDico[sampleName][line.split(",")[9]] += int(
                        line.split(",")[5].split(".")[0])
            f.close()


        data = []
        for virus in virusNameList:
            valueList = []
            for sample in sampleNameList:
                if sampleDico[sample][virus] > 0:
                  valueList.append(sampleDico[sample][virus])
                else:
                  valueList.append("none")                              # 0-> none in heatmap hoverongaps = False for no defaults values
            data.append(valueList)
        print(data)

        fig.add_trace(go.Heatmap(z=data, name=str(genus),colorscale= [[0, 'whitesmoke'], [0.33, '#fad146'], [0.66, 'orange'], [1, 'crimson']],
                                 hovertemplate='Coverage: %{z} %' + '<br>Virus: %{y}' + '<br>Sample: %{x}',
                                 y=virusNameList,
                                 x=sampleNameList,
                                 zmin=0,
                                 zmax=100,


                                 ), n, 1)

        fig.data[0].update(zmin=50, zmid=51, zmax=100 )

        fig.update_coloraxes(showscale=False)


      w= 200 + (len(sampleNameList) * 100)
      h= 200 + max(listLenVirusByGenus) * len(genusNameList) * 25

      fig.update_layout(

          autosize=False,
          width=w,
          height=h,
          margin=dict(
              l=50,
              r=50,
              b=100,
              t=100,
              pad=4
          ),
          paper_bgcolor="lightgrey",
      )

      fig.update_traces(showscale=False)
      fig.write_html(pathData + "/USERDATA/" + projectSelected +"/file.html")
      new = 2  # open in a new tab, if possible
      url=pathData + "/USERDATA/" + projectSelected +"/file.html"
     
      webbrowser.get().open("file:"+url,new=new)

      save_button3 = ttk.Button(
        main,
        text='Save HTML',
        command=html_save
      )
      save_button3.place(x=300, y=60)


      #my_label = HTMLLabel(main, html="<a href="+ url + ">click here if browser does not open</a>",state="normal", font=("Calibri", 3))
      #my_label.place(x=500, y=60)
      #fig.show()




def default():

    global projectSelected
    global projectList
    global cb1

    for widgets in main.winfo_children():
        widgets.destroy()
    topButton()
    start()

    text3Label = ttk.Label(main, text="Enter a project name:")

    text3Label.place(x=10, y=45)          # Add a sample

    cb1 = ttk.Combobox(main, values=projectList)
    cb1.place(x=215, y=45)

    projectSelected = cb1.get().split("'")[-1]
    print(projectSelected)

    enter_button = ttk.Button(   #button enter-like
        main,
        text='⏎',
        command=analyse_batch
    )

    enter_button.place(x=380, y=40)           #analyze
    tree=Image.open(pathData + "/IMG/TREE/GENcol.png")
    img = ImageTk.PhotoImage(tree)
    label = ttk.Label(main, image=img)
    label.image = img
    label.place(x=1, y=75)

    text4Label = ttk.Label(main, text="Tree of Primate Anelloviruses: \n\nClick on Alpha- Beta- or Gamma- \ntorquevirus to display sub-trees.\n\nNon-human hosts are indicated\nin parenthesis\n\nNon-ICTV approved (new) species\nbranches are colored in red.")

    text4Label.place(x=5, y=80)          # Add a legend




#button a placer sur l'image pour afficher les sous-arbres

    beta_button = ttk.Button(   
        main,
        text='Betatorquevirus',
        command=placeBeta
    )

    beta_button.place(x=485, y=280)           


    alpha_button = ttk.Button(   
        main,
        text='Alphatorquevirus',
        command=placeAlpha
    )

    alpha_button.place(x=485, y=895)    

    gamma_button = ttk.Button(   
        main,
        text='Gammatorquevirus',
        command=placeGamma
    )

    gamma_button.place(x=485, y=805)    

def placeAlpha():
    alpha=Image.open(pathData + "/IMG/TREE/ALPHA.png")
    imgAlpha = ImageTk.PhotoImage(alpha)
    labelA = ttk.Label(main, image=imgAlpha)
    labelA.image = imgAlpha
    labelA.place(x=700, y=75)

def placeBeta():
    beta=Image.open(pathData + "/IMG/TREE/BETA.png")
    imgBeta = ImageTk.PhotoImage(beta)
    labelB = ttk.Label(main, image=imgBeta)
    labelB.image = imgBeta
    labelB.place(x=700, y=75)

def placeGamma():
    gamma=Image.open(pathData + "/IMG/TREE/GAMMA.png")
    imgGamma = ImageTk.PhotoImage(gamma)
    labelG = ttk.Label(main, image=imgGamma)
    labelG.image = imgGamma
    labelG.place(x=700, y=75)


def analyse_batch():

    global fastq1Path
    global sampleList
    global projectSelected
    global sampleUniq
    global yAdd
    global pathLastDir
    global pairedLabel
    global text1Label
    global text2Label
    global minion
    global illuminaSE
    global illuminaPE
    global reset_button
    global run_button
    global makeConsensus
    global pathData
    try:
      projectSelected = cb1.get()
    except:
      pass  #if not combobox window
    # test if project exist otherwise mk directory
    if projectSelected == "":
        projectSelected="default"
    if not os.path.exists(pathData + "/USERDATA/" + projectSelected):
        os.mkdir(pathData + "/USERDATA/" + projectSelected)
        log = open(pathData + "/LOG/log.txt", "a")
        log.write("Project: USERDATA/" + projectSelected + "added" + "\n")
        log.close()

    for widgets in main.winfo_children():
        widgets.destroy()
    topButton()

    minion = tk.IntVar(main, 0)  # 1 if Nanopore data 0 otherwise         # The check mark box
    minion.set(0)

    illuminaSE = tk.IntVar(main, 1)  #if illumina Single End 1 otherwise 0
    illuminaSE.set(1)

    illuminaPE = tk.IntVar(main, 0)  # Paired End
    illuminaPE.set(0)

    makeConsensus = tk.IntVar(main, 0)  # if check make consensus
    makeConsensus.set(0)

    yAdd = 250    # height

    sampleUniq = set(listdir(pathData + "/USERDATA/" + projectSelected + "/"))

    open_button = ttk.Button(
        main,
        text='Add Fastq File(s)',
        command=select_file_batch
    )

    run_button = ttk.Button(
        main,
        text='Run',
        command=run
    )

    # checkbox

    minion_check = ttk.Checkbutton(
        main,
        text="Oxford Nanopore",
        variable=minion,
        command=excludeM
    )

    illuminaSE_check = ttk.Checkbutton(
        main,
        text="Illumina single end",
        variable=illuminaSE,
        command=excludeS
    )

    illuminaPE_check = ttk.Checkbutton(
        main,
        text="Illumina paired end",
        variable=illuminaPE,
        command=excludeP
    )

    consensus_check = ttk.Checkbutton(
        main,
        text="genarate FASTA consensus [slow mode]",
        variable=makeConsensus,
    )

    reset_button = ttk.Button(
        main,
        text='Reset',
        command=resetA
    )

    text0Label = ttk.Label(main, text="Selected project: ", foreground="gray")
    text05Label = ttk.Label(main, text=projectSelected, foreground="darkorange1")
    text1Label = ttk.Label(main, text="Select samples to add:")
    text2Label = ttk.Label(main, text="Added sample: ")
    text05Label.place(x=510, y=5)
    text0Label.place(x=350, y=5)
    text1Label.place(x=10, y=45)         # Add a sample
    open_button.place(x=10, y=70)
    run_button.place(x=10, y=195)
    text2Label.place(x=10, y=240)
    minion_check.place(x=10, y=110)       # Oxford Nanopore
    illuminaSE_check.place(x=10, y=130)
    illuminaPE_check.place(x=10, y=150)
    consensus_check.place(x=10, y=170)

    reset_button.place(x=150, y=195)            #reset

def excludeM():
    if minion.get() == 1:
        illuminaSE.set(0)
        illuminaPE.set(0)

def excludeS():
    if illuminaSE.get() == 1:
        minion.set(0)
        illuminaPE.set(0)
def excludeP():
    if illuminaPE.get() == 1:
        minion.set(0)
        illuminaSE.set(0)

#######################################################################################################################

# button
def topButton():

  # mini menubar

  file = tk.Menubutton(main, text="File", fg="gray", activebackground="darkorange1")

  file.menu = Menu(file, tearoff=0, fg="gray", activebackground="darkorange1")
  file["menu"] = file.menu
  file.place(x=0, y=0)

  file.menu.add_command(label='Project selection', command=default)
  file.menu.add_command(label='Analysis', command=analyse_batch)
  file.menu.add_command(label='Exit', command=main.destroy)

  file.menu.insert_separator(2)

  edit = tk.Menubutton(main, text="Edit", fg="gray", activebackground="darkorange1")
  edit.menu = Menu(edit, tearoff=0, fg="gray", activebackground="darkorange1")
  edit["menu"] = edit.menu
  edit.place(x=45, y=0)

  edit.menu.add_command(label='Delete project', command=deleteP)

  edit.menu.add_command(label='Delete all data', command=deleteA)

  view = tk.Menubutton(main, text="View", fg="gray", activebackground="darkorange1")

  view.menu = Menu(view, tearoff=0, fg="gray", activebackground="darkorange1")
  view["menu"] = view.menu
  view.place(x=90, y=0)
  view.menu.add_command(label='Data', command=dataA)

# run the application
default()

main.mainloop()
