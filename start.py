import os
import time
from os import listdir
import mappy as mp
from statistics import median
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter.messagebox import showinfo
from tkinter import *
new = tk.Tk()
new.title('Anelloviruses')
new.geometry("800x500")
global workingRep
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
global pb
global minion

minion = tk.IntVar(new, 0)  # 1 if Nanopore data 0 otherwise
minion.set(0)
pathLastDir = "/home/florian/PycharmProjects/interface/"

sampleUniq = set(listdir("USERDATA/")) #list of all "samples" from past session - for each sample one directory in USERDATA/
print(sampleUniq)

workingRep = "/home/florian/PycharmProjects/interface/working"

yAdd = 150
sampleList = []

fastq1Path = "<empty-mandatory>"
fastq2Path = "<empty>"

fastq1 = StringVar()
fastq2 = StringVar()

fastq1.set(fastq1Path)
fastq2.set(fastq2Path)

fq1 = Label(new, textvariable=fastq1, bg="red", fg="black", relief=RIDGE)
fq2 = Label(new, textvariable=fastq2, bg="yellow", fg="black", relief=RIDGE)
pairedLabel = Label(new, text="(if paired)")
text1Label = Label(new, text="Add a sample:")
text2Label = Label(new, text="Sample list:")


def select_file1():
    global fastq1Path
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

    if fastq1Path.rstrip().split(".")[-1] == "fastq":
        fq1.config(bg="green", fg="white")
    else:
        fq1.config(bg="yellow", fg="black")

def select_file2():
    global fastq2Path
    filetypes = (
        ('Fastq files', '*.fastq'),
        ('All files', '*.*')
    )

    filename = fd.askopenfilename(
        title='Open a file',
        initialdir='/',
        filetypes=filetypes)

    showinfo(
        title='Selected File',
        message=filename
    )
    fastq2Path = str(filename)

    fastq2.set(fastq2Path)

    if fastq2Path.rstrip().split(".")[-1] == "fastq":
        fq2.config(bg="green", fg="white")
    else:
        fq2.config(bg="yellow", fg="black")

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
        os.mkdir("USERDATA/" + sampleName)
        sample = [sampleName, fastq1Path, fastq2Path, yAdd + 20, minion.get()]
        sampleList.append(sample)
        labelListSample = Label(new, text=sampleName)
        yAdd += 20
        labelListSample.place(x=5, y=yAdd)
        sampleUniq.add(sampleName)
    print(minion.get())
    fastq1.set("<empty-mandatory>")
    fq1.config(bg="red", fg="black")
    fastq2.set("<empty>")
    fq2.config(bg="yellow", fg="black")

def run():
    global pb
    print(sampleList)

    for s in sampleList:
        pb = ttk.Progressbar(new, orient='horizontal', mode='determinate', length=280)  # Progress bar
        pb.place(x=350, y=s[3])
        new.update()
        curr = time.time()
        pb['value'] += 1
        new.update()
        mapping(s[1], "testdb/Anellovirus_2022.0.fasta", s[0])
        pb['value'] += 1
        new.update()

        print(time.time() - curr)
    save_button.place(x=200, y=120)


def mapping(pathToFastq, db, nameS, type="single"):
    global pb

    genome_ref_covered = dict()
    genome_ref_class = dict()  # key acc number; value list of info from the fasta header
    preset = None
    if type == "Nanopore":
        preset = "map-ont"     # for Oxford Nanopore read mapping
    elif type== "single":
        preset ="sr"           # single-end short reads
    a = mp.Aligner(db, preset=preset)

    pb['value'] += 1
    new.update()

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
                # info classification in a separated dict.
                genome_ref_class[mgRef] = [hit.ctg.split(",")[1], hit.ctg.split(",")[2], hit.ctg.split(",")[3],
                                           hit.ctg.split(",")[4], hit.ctg.split(",")[5]]
                                                 #TODO keep only acc Number in fasta db make a correspendence table

            genome_ref_covered[mgRef][1] += 1
            for i in range(hit.r_st, hit.r_en):             # r_st  = ref start match, r_en -> end
                genome_ref_covered[mgRef][0][i] += 1
    pb['value'] += 88
    new.update()

    for l in genome_ref_covered:
        genome_ref_covered[l].append(genome_ref_covered[l][2] - genome_ref_covered[l][0].count(0))  # 3 cov
        genome_ref_covered[l].append(
            round((genome_ref_covered[l][2] - genome_ref_covered[l][0].count(0)) / genome_ref_covered[l][2] * 100,
                  2))  # 4 % cov
        genome_ref_covered[l].append(median(genome_ref_covered[l][0]))  # 5 cov depth


    pb['value'] += 2
    new.update()

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
    new.update()

    resultSp = set()
    for sp in listSpecies:  # for each sp keep the one with the best coverage
        bestCov = 49
        for acc in dictSpecies[sp]:
            if genome_ref_covered[acc][4] > bestCov:  # tie ?
                keep = acc
                bestCov = genome_ref_covered[acc][4]
        resultSp.add(keep)

    pb['value'] += 2
    new.update()

    fichierCSV = open("USERDATA/" + nameS + "/species.csv", "w")
    #Sample Name, ACC. NUMBER, Reads, ref_len, cov, %cov, depth (median), GENUS, GROUP, SPECIES, GENOTYPE, HOST
    for accResultSp in resultSp:
        fichierCSV.write(nameS + ", " + accResultSp + ", " + str(genome_ref_covered[accResultSp][1]) + ", " + str(genome_ref_covered[accResultSp][2]) + ", " + str(
            genome_ref_covered[accResultSp][3]) + ", " + str(genome_ref_covered[accResultSp][4]) + ", " + str(
            genome_ref_covered[accResultSp][5]) + ", " + genome_ref_class[accResultSp][0].split("=")[1]+ ", " + genome_ref_class[accResultSp][1].split("=")[1]+ ", " + genome_ref_class[accResultSp][2].split("=")[1] + ", " + genome_ref_class[accResultSp][3].split("=")[1] + ", " + genome_ref_class[accResultSp][4].split("=")[1] + "\n")
    fichierCSV.close()
    pb['value'] += 4
    new.update()

    return ("done")

def file_save():
    nameNewCsv = fd.asksaveasfile(mode='w',defaultextension=".csv")
    text2save="Sample Name, ACC. NUMBER, Reads, ref_len, cov, %cov, depth (median), GENUS, GROUP, SPECIES, GENOTYPE, HOST\n"
    for s in sampleList:
        fichier = open("USERDATA/" + s[0] + "/species.csv", "r")
        for lane in fichier.read():

            text2save += lane
        fichier.close()

    nameNewCsv.write(text2save)
    nameNewCsv.close

# button
open_button_fq1 = ttk.Button(
    new,
    text='Select Fastq File 1',
    command=select_file1
)

open_button_fq2 = ttk.Button(
    new,
    text='Select Fastq File 2',
    command=select_file2
)

add_sample_button = ttk.Button(
    new,
    text='Add',
    command=add_sample
)

run_button = ttk.Button(
    new,
    text='Run',
    command=run
)

save_button = ttk.Button(
    new,
    text='Save',
    command=file_save
)

#checkbox

minion_check = ttk.Checkbutton(
    new,
    text="Oxford Nanopore",
    variable=minion
    )



text1Label.place(x=5, y=15)
open_button_fq1.place(x=5, y=40)
fq1.place(x=230, y=45)
open_button_fq2.place(x=5, y=70)
fq2.place(x=230, y=75)
pairedLabel.place(x=150, y=75)
add_sample_button.place(x=5, y=120)
run_button.place(x=100, y=120)
text2Label.place(x=5, y=145)
minion_check.place(x=230, y=15)

# run the application

new.mainloop()
