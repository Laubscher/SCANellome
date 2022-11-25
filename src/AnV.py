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
#from tkinter import *
import Dicodb                             # module that contain the database in a python dico
import base64
import img
from ttkthemes import ThemedTk

# The software window
#main = tk.Tk()
main = ThemedTk(theme="ubuntu", background=True)

main.title('                                                                      AnV                                                                         v. 0.0.2')
main.geometry("800x500")

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
      main.iconphoto(False, icon)
  except:
      pass

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
    fasta=open(pathData + "/DATABASE/Anello.fasta", "w")
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
              print("P-aired")
          else:
              type="Single"
              print("single")
          mapping(s[1], pathData + "/DATABASE/Anello.fasta", s[0], type, s[2])


        except:
          log.write("Error during mapping.. sample:" + str(s[0]) + "\n files may be corrupted!" + "\n")
          tk.messagebox.showinfo("Error..", "Error during mapping.. sample:" + str(s[0]) + "\n files may be corrupted!")
          pb['value'] = 0

        else:
          pb['value'] += 1
          main.update()
          log.write("Mapping time sample: " + str(s[0]) + " " + str(time.time() - curr) + "\n")

    save_button = ttk.Button(
        main,
        text='Save',
        command=file_save
    )
    save_button.place(x=290, y=195)
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
            if genome_ref_covered[acc][4] > bestCov:  # tie ?
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
        ('Fastq files', '*.fastq'),
        ('All files', '*.*')
    )
    filenames = fd.askopenfilenames(
        title='select fastq',
        initialdir=pathLastDir,
        filetypes=filetypes)

    for path in filenames:
        fastq1Path=str(path)
        if fastq1Path.rstrip().split(".")[-1] == "fastq":           #delete this ?
            pass
        else:
            pass # add a window error here
        add_sample_batch(fastq1Path)

    pathL = fastq1Path.split("/")           #for remember the path of the directory
    pathL.pop()
    pathLastDir = "/".join(pathL)

def add_sample_batch(fastq1Path):
    global yAdd
    global sampleList
    global sampleUniq

    isSample = True
    if illuminaPE.get() == 1:
        if fastq1Path.split("_R1_001")[-1] != ".fastq":
          isSample=False
        else:
          fastq2Path = fastq1Path.split("_R1_001.fastq")[0] + "_R2_001.fastq"
          sampleName = fastq1Path.split("/")[-1].split("_R1_001.fastq")[0]
    else :
        fastq2Path = "<empty>"  # Todo add it if R2 found (and /or check box)
        sampleName = fastq1Path.split("/")[-1].split(".fastq")[0]
    if isSample:
        suffix = fastq1Path.split(".")[-1]
    # check if fastq

    # check if sample name is uniq
        if sampleName in sampleUniq:
          print("Sample name not uniq!!")
          tk.messagebox.showinfo("Sorry can't add your sample..", '\n"' + sampleName + '"\n\nThe sample name already exist in this project.\n')

        elif suffix != "fastq":
          print("Sample suffix is not fastq!!")
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



'''############################ one by one ###############################################################################

def select_file1():
    global fastq1Path
    global pathLastDir
    filetypes = (
        ('Fastq files', '*.fastq'),
        ('All files', '*.*')
    )
    filename = fd.askopenfilename(
        title='Open a file',
        initialdir=pathLastDir,
        filetypes=filetypes)

    fastq1Path = str(filename)
    print(fastq1Path)

    fastq1.set(fastq1Path)
    pathL = fastq1Path.split("/")
    pathL.pop()
    pathLastDir = "/".join(pathL)

    if fastq1Path.rstrip().split(".")[-1] == "fastq":
        fq1.config(background="darkorange1", foreground="white")
    else:
        fq1.config(background="yellow", foreground="black")

def select_file2():
    global fastq2Path
    global pathLastDir
    filetypes = (
        ('Fastq files', '*.fastq'),
        ('All files', '*.*')
    )

    filename = fd.askopenfilename(
        title='Open a file',
        initialdir=pathLastDir,
        filetypes=filetypes)

    fastq2Path = str(filename)

    fastq2.set(fastq2Path)

    pathL = fastq2Path.split("/")
    pathL.pop()
    pathLastDir = "/".join(pathL)

    if fastq2Path.rstrip().split(".")[-1] == "fastq":
        fq2.config(background="orange", foreground="white")
    else:
        fq2.config(background="yellow", foreground="black")

def add_sample():
    global yAdd
    sampleName = fastq1Path.split("/")[-1].split(".fastq")[0]
    # check if sample name is uniq
    suffix  = fastq1Path.split(".")[-1]
    # check if fastq

    if sampleName in sampleUniq:
        print("Sample name not uniq!!")
        tk.messagebox.showinfo("Sorry can't add your sample..", "The sample name is not uniq!")

    elif suffix != "fastq":
        print("Sample suffix is not fastq!!")
        tk.messagebox.showinfo("Sorry can't add your sample..", "The file type is not fastq!")

    else:
        os.mkdir(pathData + "/USERDATA/" + projectSelected + "/" + sampleName)
        sample = [sampleName, fastq1Path, fastq2Path, yAdd + 20, minion.get()]
        sampleList.append(sample)
        labelListSample = ttk.Label(main, text=sampleName)
        yAdd += 20
        labelListSample.place(x=10, y=yAdd)
        sampleUniq.add(sampleName)
    print(minion.get())
    fastq1.set("<empty-mandatory>")
    fq1.config(background="gray", foreground="black")
    fastq2.set("<empty>")
    fq2.config(background="yellow", foreground="black")

#######################################################################################################################'''


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

def deleteA():
    #TODO add a confirmation window
    global sampleUniq
    sampleUniq = []
    shutil.rmtree(pathData)
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

def dataA():
    global projectList
    #global sampleUniq # list of all directory one by sample in USERDATA
    start()
    for widgets in main.winfo_children():
        widgets.destroy()
    topButton()

    projectSelected = cb1.get()

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

    cb1 = ttk.Combobox(main, values=projectList)#, width=7)
    cb1.place(x=215, y=45)

    projectSelected = cb1.get().split("'")[-1]
    print(projectSelected)

    enter_button = ttk.Button(   #button enter-like
        main,
        text='⏎',
        command=analyse_batch
    )
    enter_button.place(x=380, y=40)           #analyze

'''def analyse():

    global fastq1Path
    global fastq2Path
    global fastq1
    global fastq2
    global fq1
    global fq2
    global sampleList
    global projectSelected
    global sampleUniq
    global yAdd
    global pathLastDir
    global pairedLabel
    global text1Label
    global text2Label
    global minion
    global pathData

    projectSelected = cb1.get()
    # test if project exist otherwise mk directory

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

    yAdd = 200

    sampleUniq = set(listdir(pathData + "/USERDATA/" + projectSelected + "/"))

    open_button_fq1 = ttk.Button(
        main,
        text='Select Fastq File 1',
        command=select_file1
    )

    open_button_fq2 = ttk.Button(
        main,
        text='Select Fastq File 2',
        command=select_file2
    )

    add_sample_button = ttk.Button(
        main,
        text='Add',
        command=add_sample
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
        variable=minion
    )
    fastq1Path = "<empty-mandatory>"
    fastq2Path = "<empty>"

    fastq1 = tk.StringVar()
    fastq2 = tk.StringVar()

    fastq1.set(fastq1Path)
    fastq2.set(fastq2Path)

    fq1 = ttk.Label(main, textvariable=fastq1, background="gray", foreground="black", relief="ridge")
    fq2 = ttk.Label(main, textvariable=fastq2, background="yellow", foreground="black", relief="ridge")

    pairedLabel = ttk.Label(main, text="(if paired)")

    text1Label = ttk.Label(main, text="Add a sample:")
    text2Label = ttk.Label(main, text="Sample list:")

    text1Label.place(x=10, y=45)          # Add a sample
    open_button_fq1.place(x=10, y=70)     #
    open_button_fq2.place(x=10, y=105)    #
    fq1.place(x=330, y=75)                #
    fq2.place(x=330, y=110)               #
    pairedLabel.place(x=220, y=110)       #
    add_sample_button.place(x=10, y=155)  #
    run_button.place(x=150, y=155)        #
    text2Label.place(x=10, y=200)         #
    minion_check.place(x=330, y=50)       # Oxford Nanopore

######################################################################################################################'''

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

    global pathData

    projectSelected = cb1.get()
    # test if project exist otherwise mk directory

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

    text1Label = ttk.Label(main, text="Select samples to add:")
    text2Label = ttk.Label(main, text="Sample list: ")

    text1Label.place(x=10, y=45)         # Add a sample
    open_button.place(x=10, y=140)       #
    run_button.place(x=10, y=195)        #
    text2Label.place(x=10, y=240)        #
    minion_check.place(x=10, y=70)       # Oxford Nanopore
    illuminaSE_check.place(x=10, y=90)
    illuminaPE_check.place(x=10, y=110)


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
  delete_button = ttk.Button(
    main,
    text='Delete all data',
    command=deleteA
  )

  reset_button = ttk.Button(
    main,
    text='Reset',
    command=resetA
  )

  analyse_button = ttk.Button(
    main,
    text='Analysis',
    command=analyse_batch
  )

  data_button = ttk.Button(
    main,
    text='Data',
    command=dataA
  )

  project_button = ttk.Button(
    main,
    text='Project selection',
    command=default
  )

  delete_button.place(x=625, y=-5)           #delete all data
  analyse_button.place(x=-5, y=-5)           #analyze
  reset_button.place(x=125, y=-5)            #reset
  data_button.place(x=500, y=-5)             #data
  project_button.place(x=260, y=-5)          #project selection


# run the application

#topButton()
#start()
default()

main.mainloop()
