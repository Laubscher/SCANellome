#!/usr/bin/env python3
import plotly.express as px

import glob

# file name in a list, then read once all the file to make a set of virus name in all the project # ? devide by genus
# then reread all the file to make a dictionary where all key are sampleID and values are
# an other dictionary with keys as virusNames and values as nb of reads


fileList = glob.glob("*/species.csv")

sampleDico = {}

virusNameSet = set()
sampleNameSet = set()

for file in fileList:
    f = open(file, "r")
    for line in f:
        virusNameSet.add(line.split(",")[9])
        sampleNameSet.add(line.split(",")[0])
    f.close()

# virusNameSet.remove('Virus')

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
            sampleDico[sampleName][line.split(",")[9]] += int(line.split(",")[5].split(".")[0])  # attention arrondir
    f.close()

sampleNameList = list(sampleNameSet)
virusNameList = list(virusNameSet)

data = []
for sample in sampleNameList:
    valueList = []
    for virus in virusNameList:
        valueList.append(sampleDico[sample][virus])
    data.append(valueList)

fig = px.imshow(data,
                labels=dict(x="Viruses", y="Samples", color="% Coverage"),
                x=virusNameList,
                y=sampleNameList
                )

fig.update_xaxes(side="top")
fig.show()