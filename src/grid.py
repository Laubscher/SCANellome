#!/usr/bin/env python3

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import glob

# file name in a list, then read once all the file to make a set of virus name in all the project # ? devide by genus
# then reread all the file to make a dictionary where all key are sampleID and values are
# an other dictionary with keys as virusNames and values as nb of reads

fileList = glob.glob("*/species.csv")

sampleDico = {}

genusNameSet = set()

sampleNameSet = set()

for file in fileList:
    f = open(file, "r")
    for line in f:
        genusNameSet.add(line.split(",")[7])
        sampleNameSet.add(line.split(",")[0])
    f.close()

sampleNameList = list(sampleNameSet)
genusNameList = sorted(list(genusNameSet))

# virusNameSet.remove('Virus')


fig = make_subplots(len(genusNameList),1)   # make one subplot for each genus
n=0  # counter n-ème subplot

for genus in genusNameList:
  n+=1
  virusNameSet = set()
  for file in fileList:
        f = open(file, "r")
        for line in f:
            if genus == line.split(",")[7]:
              virusNameSet.add(line.split(",")[9])
        f.close()
  virusNameList = list(virusNameSet)
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


# mettre none à la place de 0 et dans le heatmap hoverongaps = False


  data = []
  for virus in virusNameList:
    valueList = []
    for sample in sampleNameList:
        valueList.append(sampleDico[sample][virus])
    data.append(valueList)
  print(data)



  fig.add_trace( go.Heatmap(z=data,name=str(genus),
                hovertemplate='Coverage: %{z} %'+'<br>Virus: %{y}'+'<br>Sample: %{x}',
                y=virusNameList,
                x=sampleNameList,
                ), n,1)




fig.write_html("file.html")
#fig.update_xaxes(side="top")
fig.show()
# labels=dict(x="Viruses", y="Samples", z="% Coverage"),