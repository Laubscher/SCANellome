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

from ttkthemes import ThemedTk

# The software window

#main = tk.Tk()
main = ThemedTk(theme="ubuntu", background=True)

main.title('AnV')
main.geometry("800x500")

global pb

def start():
  #Variables
  global fastq1Path
  global fastq2Path
  global fastq1
  global fastq2
  global fq1
  global fq2
  global sampleList
  global sampleUniq
  global yAdd
  global pathLastDir
  global pairedLabel
  global text1Label
  global text2Label
  global minion
  global pathData

  pathData = os.path.expanduser("~/.AnV")                               # check if .AnV otherwise mk it
  if not os.path.exists(pathData):
    os.mkdir(pathData)

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
    log.write("USERDATA set up" + "\n")

  pathLastDir = os.path.expanduser("~")

  sampleUniq = set(listdir(pathData + "/USERDATA/"))     # list of all "samples" from past session - for each sample one directory in USERDATA/

  log.write("List of samples found in USERDATA/: " + str(sampleUniq) + "\n")

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

    showinfo(
        title='Selected File',
        message=filename
    )

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

    showinfo(
        title='Selected File',
        message=filename
    )
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
        os.mkdir(pathData + "/USERDATA/" + sampleName)
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

def run():
    global pb
    global main
    log = open(pathData + "/LOG/log.txt", "a")
    log.write("Start mapping samples ... " + "\n")

    for s in sampleList:
        log.write("Mapping sample: " + str(s[0]) + "\n")
        pb = ttk.Progressbar(main, orient='horizontal', mode='determinate', length=280)  # Progress bar
        pb.place(x=350, y=s[3])
        main.update()
        pb['value'] += 1
        main.update()
        try:
          curr = time.time()
          mapping(s[1], pathData + "/DATABASE/Anello.fasta", s[0])

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
    save_button.place(x=290, y=155)
    tk.messagebox.showinfo("Analysis completed", "Analysis completed" + "\n")
    log.close()

def mapping(pathToFastq, db, nameS, type="single"):
    global pb

    genome_ref_covered = dict()
    genome_ref_class = Dicodb.db                                    # key acc number; value list of info from the fasta header
    preset = None
    if type == "Nanopore":
        preset = "map-ont"     # for Oxford Nanopore read mapping
    elif type == "single":
        preset = "sr"           # single-end short reads
    a = mp.Aligner(db, preset=preset)

    pb['value'] += 1
    main.update()

    for name, seq, qual in mp.fastx_read(pathToFastq):
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
    pb['value'] += 88
    main.update()

    for l in genome_ref_covered:
        genome_ref_covered[l].append(genome_ref_covered[l][2] - genome_ref_covered[l][0].count(0))  # 3 cov
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

    fichierCSV = open(pathData + "/USERDATA/" + nameS + "/species.csv", "w")
    #Sample Name, ACC. NUMBER, Reads, ref_len, cov, %cov, depth (median), GENUS, GROUP, SPECIES, GENOTYPE, HOST
    for accResultSp in resultSp:
        fichierCSV.write(nameS + ", " + accResultSp + ", " + str(genome_ref_covered[accResultSp][1]) + ", " + str(genome_ref_covered[accResultSp][2]) + ", " + str(
            genome_ref_covered[accResultSp][3]) + ", " + str(genome_ref_covered[accResultSp][4]) + ", " + str(
            genome_ref_covered[accResultSp][5]) + ", " + genome_ref_class[accResultSp][0].split("=")[1]+ ", " + genome_ref_class[accResultSp][1].split("=")[1]+ ", " + genome_ref_class[accResultSp][2].split("=")[1] + ", " + genome_ref_class[accResultSp][3].split("=")[1] + ", " + genome_ref_class[accResultSp][4].split("=")[1] + "\n")
    fichierCSV.close()
    pb['value'] += 4
    main.update()

    return ("done")

def file_save():
    nameNewCsv = fd.asksaveasfile(mode='w', defaultextension=".csv")
    text2save="Sample Name, ACC. NUMBER, Reads, ref_len, cov, %cov, depth (median), GENUS, GROUP, SPECIES, GENOTYPE, HOST\n"
    for s in sampleList:
        try:
            fichier = open(pathData + "/USERDATA/" + s[0] + "/species.csv", "r")    # if error during mapping file will not existe
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
        os.rmdir(pathData + "/USERDATA/" + str(sample[0]))
    sampleUniq = []
    start()
    for widgets in main.winfo_children():
        widgets.destroy()
    topButton()

def analyse():

    global fastq1Path
    global fastq2Path
    global fastq1
    global fastq2
    global fq1
    global fq2
    global sampleList
    global sampleUniq
    global yAdd
    global pathLastDir
    global pairedLabel
    global text1Label
    global text2Label
    global minion
    global pathData
    global minion

    minion = tk.IntVar(main, 0)  # 1 if Nanopore data 0 otherwise         # The check mark box
    minion.set(0)

    yAdd = 200
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
    fq1.place(x=330, y=75)
    fq2.place(x=330, y=110)               #
    pairedLabel.place(x=220, y=110)       #
    add_sample_button.place(x=10, y=155)  #
    run_button.place(x=150, y=155)        #
    text2Label.place(x=10, y=200)         #
    minion_check.place(x=330, y=50)       # Oxford Nanopore


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
    command=analyse
  )

  delete_button.place(x=625, y=-5)           #delete all data
  analyse_button.place(x=-5, y=-5)           #analyze
  reset_button.place(x=125, y=-5)             #reset

# run the application

topButton()
start()
main.mainloop()
