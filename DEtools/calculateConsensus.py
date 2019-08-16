#!/usr/bin/env python3
import os
import sys
import json

##python consensusTable_new.py jobID /dbs/conDE/upload outputFile

#jobid media_root outputfile
def tableConsensus(input_folder):
    outputTable = os.path.join(input_folder,"consensus.tsv")
    jobDir = os.path.join(input_folder,"consensus")
    methodsFiles = next(os.walk(jobDir))[2]  
    
    consensus = {}
    for method in methodsFiles:
        methodFile = open(jobDir+"/"+method)
        methodName = method.split(".")[0]

        header = methodFile.readline().strip().split("\t")
        nameIndex = header.index("name")+1
        FoldChangeIndex = header.index("FoldChange")+1

        for line in methodFile:
            line = line.strip().split("\t")
            name = line[nameIndex]
            FoldChange = float(line[FoldChangeIndex])
            
            try:
                data = consensus[name]
            except:
                data = [0,0,[],[]]

            countOver = data[0]
            countInfra = data[1]
            methodsOver = data[2]
            methodsInfra = data[3]

            if FoldChange>1:
                countOver = countOver+1
                methodsOver.append(methodName)
            else:
                countInfra = countInfra + 1
                methodsInfra.append(methodName)
            data = [countOver,countInfra,methodsOver,methodsInfra]
            consensus[name] = data

    headerOutfile = "name\t#Over\t#Infra\t#Consensus\tExpression\tconsensusMethods\n"

    outFile = open(outputTable,'w')
    outFile.write(headerOutfile)

    for element in consensus:
        countOver = str(consensus[element][0])
        countInfra = str(consensus[element][1])
        if countOver>countInfra:
            label = "Over"
            methodsConsensus = ",".join(consensus[element][2])
            countConsensus = str(countOver)
        elif countInfra>countOver:
            label = "Infra"
            methodsConsensus = ",".join(consensus[element][3])
            countConsensus = str(countInfra)
        else:
            label = "Both"
            methodsConsensus = "-"
            countConsensus = str(countInfra)
        toWrite = element+"\t"+countOver+"\t"+countInfra+"\t"+countConsensus+"\t"+label+"\t"+methodsConsensus+"\n"
        outFile.write(toWrite)

    outFile.close()

# with open(sys.argv[1], "r") as jf:
with open(sys.argv[1], "r") as jf:
    input_dict = json.load(jf)
# input_dict = json.load()

tableConsensus(input_dict.get("folder"))