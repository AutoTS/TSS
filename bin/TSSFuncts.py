#!/bin/env python3
#small comment
#Written By Taylor Nielson

import sys, os,re, shutil, numpy as np, time, subprocess

#We need to add a crest_selectivity default
def default():
	inputs = {}
	iF = open(os.path.expanduser("~/TSS/bin/.default"), "r")
	line = iF.readline()
	while line:
		if "spin" in line:
			inputs["spin"] = line.split(':')[1]
		elif "charge" in line:
			inputs["charge"] = line.split(':')[1]
		elif "basis" in line:
			inputs["basis"] = line.split(':')[1]
		elif "method" in line:
			inputs["method"] = line.split(':')[1]
		elif "batch" in line:
			inputs["batch"] = line.split(':')[1]
		elif "metal" in line:
                	inputs["mbasis"] = line.split(':')[1]
		elif "denfit" in line:
			inputs["denfit"] = line.split(':')[1]
		line = iF.readline()
	iF.close()
	inputs["opt"] = "modred"
	return inputs


#Need to include a crest_selectivity parameter
def parseInput(inputFile, inputs):
	name = inputFile.split('.')[0]
	extension = inputFile.split('.')[1]
	if extension == "in":
		pass
	else:
		sys.exit("input file wasn't a .in file")
	iF = open(inputFile, "r")
	coords = []
	line = iF.readline()
	while line:
		#Required - Default in file
		if "spin" in line:
			inputs["spin"] = line.split(':')[1]
		 #Required - Default in file
		elif "charge" in line:
			inputs["charge"] = line.split(':')[1]
		 #Required - Default in file
		elif "basis" in line:
			inputs["basis"] = line.split(':')[1]
		 #Required - Default in file
		elif "method" in line:
			inputs["method"] = line.split(':')[1]
		 #Required - Gaussian has a default
		elif "temperature" in line:
			inputs["temperature"] = line.split(':')[1]
		#Not Required, Gaussian Defaults to Gas
		elif "solvent" in line:
			inputs["solvent"] = line.split(':')[1]
		 #Required - Default in file
		elif "batch" in line:
			inputs["batch"] = line.split(':')[1]
		 #Required - Default in file
		elif "metal" in line:
			inputs["mbasis"] = line.split(':')[1]
		#Required, default not in file
		elif "library" in line:
			inputs["library"] = line.split(':')[1]
		#Not Required, Gaussian Defaults to gas
		elif "solvent_model" in line:
			inputs["solvent_model"] = line.split(':')[1]
		 #Required - Default in file
		elif "denfit" in line:
                	inputs["denfit"] = line.split(':')[1]
		elif "memory" in line:
			inputs["mem"] = line.split(':')[1]
		elif "time" in line:
			inputs["time"] = line.split(':')[1]
		#Not Required
		elif "subtract" in line:
			parseChanges(line, "subtract", inputs)
		#Not Required
		elif "substitute" in line:
			parseChanges(line, "substitute", inputs)
		#Not Required
		elif "add" in line:
			parseChanges(line, "add", inputs)
		elif "difficulty" in line:
			inputs["difficulty"] = line.split(':')[1]
		elif "conformation_leniency" in line:
			inputs["leniency"] = line.split(':')[1]
		line = iF.readline()
	iF.close()
	coords = getCoords(inputs)
	return inputs, coords


def getCoords(inputs):
	coords = []
	iF = open(os.path.expanduser("~/TSS/libs/base_templates/" + inputs["library"].strip()), "r")
	line = iF.readline()
	line = iF.readline()
	line = iF.readline()
	while line:
		coords.append(line)
		line = iF.readline()
	return coords

def parseChanges(line, option, inputs):
	split_line = line.split(':')[1]
	split_line = split_line[1:-2]
	split_again = split_line.split(',')
	if option == "subtract":
		inputs["subtract"] = [val for val in split_again]
	elif option == "substitute":
		inputs["substitute"] = [val for val in split_again]
	elif option == "add":
		inputs["add"] = [val for val in split_again]


def buildCom(inputs, coords, f_name):
	oF = open(f_name, 'w')
	oF.write("%nprocshared=1\n")
	oF.write("%mem="+ str(int(inputs["mem"])-1) + "GB\n")
	oF.write("#opt=" + inputs["opt"].strip() + " freq=noraman " + inputs["method"].strip() + "/" + "genecp"  + " integral = ultrafine")
	oF.write("\n\n")
	oF.write("TSS")
	oF.write("\n\n")
	oF.write(inputs["charge"].strip() + " " + inputs["spin"])
	for coord in coords:
		oF.write(coord)
	writeFreezes(oF, coords, inputs)
	oF.write("\n")
	oF.write("\n")
	oF.write("\n")
	oF.write("\n")
	oF.close()


def writeFreezes(outFile,coords, inputs):
	Freezes = []
	libFile = open(os.path.expanduser("~/TSS/libs/base_templates/" + inputs["library"].strip()), "r")
	libFile.readline()
	freeze_line = libFile.readline()
	freeze_line = freeze_line[2:]
	freeze_line = freeze_line.split(';')
	freeze_line.pop(-1)
	for i in range(0, len(freeze_line)):
        	val = freeze_line[i].split('-')
        	Freezes.append(val[0])
        	Freezes.append(val[1])
	outFile.write("\n")
	for i in range(0,len(Freezes)-1,2):
        	outFile.write("B " + str(int(Freezes[i]) + 1) + " " + str(int(Freezes[i+1]) + 1) + " F\n")
	writeGenecp(libFile, outFile, coords, inputs)	

def writeGenecp(libFile, outFile,coords, inputs):
	metals, non_metals = getAtomTypes(coords)
	outFile.write("\n")
	for val in metals:
        	outFile.write(val + " ")
	outFile.write("0\n")
	outFile.write(inputs["mbasis"].strip() + "\n****\n")
	for val in non_metals:
        	outFile.write(val + " ")
	outFile.write("0\n" + inputs["basis"].strip() + "\n****\n\n")
	for val in metals:
        	outFile.write(val + " ")
	outFile.write("0\n")
	outFile.write(inputs["mbasis"].strip())
		

def getAtomTypes(coords):
	metals = set()
	non_metals = set() 
	for coord in coords:
        	if "Fe" in coord.split()[0]:
                	metals.add("Fe")
        	else:
                	non_metals.add(coord[0].strip())
	return metals, non_metals
	


def buildInputs(library, library_location):
	if library is True:
		buildLibraryInputs(library_location)
	else:
		buildGuessedInputs()


#Copies the base_input as defined in the input file and modifies it to build a new xyz file
#Currently only does subractions
def buildLibraryInputs(lib_location):
	shutil.copy(os.path.expanduser("~/TSS/libs/base_templates/" + inputs["library"].strip()), "temp.xyz")
	tempF = open("temp.xyz", "r")
	finalF = open("crest_conformers.xyz", "w")
	atomNumber = int(tempF.readline())
	atomNumber -= len(inputs["subtract"])
	#need to include adding and substituting
	tempF.readline()	
	finalF.write(str(atomNumber) + "\n\n")
	for i in range(1, atomNumber + 1):
		line = tempF.readline()
		if str(i) in inputs["subtract"]:
			pass
		#elif i in inputs["substitute"]:
		#elif i in inputs["add"]:
		else:
			finalF.write(line)


def modredCrest(crest_file, inputs):
	os.chdir("modred")
	iF = open("../" + crest_file, 'r')
	num_structures = 1
	line = iF.readline()
	while line:
		line = iF.readline()
		line = iF.readline()
		coords = []
		while line:
			if line.strip()[0].isalpha():
				coords.append(line)
				line = iF.readline()
			else:
				break
		oF_name = "conf" + str(num_structures) + ".com"
		buildCom(inputs, coords, oF_name)
		num_structures += 1
	os.chdir("../")
# def modredRangeCreation():
	#This will take the current modreds and create multiple ones with differing frozen bond lengths

def logtoxyz(f_name):
	inFile = open(f_name, 'r')
	iF = inFile.readlines()
	myLine = 0
	for i, line in enumerate(iF):
		if 'Standard orientation' in line:
			myLine = i
	coords = []
	done = False
	i = myLine + 5
	myRegex = r'\s*\d*\s*(\d*)\s*\d*\s*(.*\s*.*\s*.*)'
	while not done:
		if '--' in iF[i]:
			break
		l = re.findall(myRegex, iF[i], flags=0)	
		line = str(l[0][0]) + '\t' + str(l[0][1])
		coords.append(line)
		i += 1
	inFile.close()
	return coords

def gaussianProcesses():
	commands = []
	switched = []
	optType = []	
	processes = []
	file_names = []
	args = sys.argv
	allDone = False
	os.chdir("modred")
	for file in os.listdir(os.getcwd()):
        	file_names.append(str(file.split('.')[0]))
        	commands.append(['/apps/gaussian16/B.01/AVX2/g16/g16', file])
        	optType.append("modred")
        	switched.append(0)
	for com in commands:
        	processes.append(subprocess.Popen(com))
	while not allDone:
        	allDone = True
        	time.sleep(10)
        	i = 0
        	drawStatus(file_names,processes, optType, switched)
        	for p in processes:
                	if p.poll() is None:
                        	allDone = False
                	else:
                        	if (not switched[i]):
                                	hasNeg = checkNegVib(file_names[i] + ".log")
                                	if hasNeg:
                                        	allDone = False
                                        	coords = logtoxyz(file_names[i] + ".log")
                                        	os.chdir('../gaussianTS/')
                                        	buildCom(inputs, coords, file_names[i] + ".com")
                                        	processes[i] = subprocess.Popen(['/apps/gaussian16/B.01/AVX2/g16/g16', (file_names[i] + ".com")])
                                        	optType[i] = "TS Calc"
                                       		os.chdir('../modred')
                                        	switched[i] = 1
                                	else:
                                        	optType[i] = "killed"
                	i += 1
	drawStatus(file_names,processes, optType, switched)
	print("all done")

def makeDirectories():
	os.mkdir("modred")
	os.mkdir("gaussianTS")

def drawStatus(file_names, processes, optType, switched):
	os.chdir('../')
	statusFile = open(".status", "w")
	i = 0
	for p in processes:
                        if p.poll() is None:
                                statusFile.write(file_names[i] + "          ------> Running " + optType[i] + "\n\n")
                        else:
                                #Check for negative frequencyi
                        	if(not switched[i]):
                                	if optType[i] is "killed":
                                        	statusFile.write(file_names[i] + "          ------> Killed - No Negative Vibration at end of Modred\n\n")
                                	else:
                                        	statusFile.write(file_names[i] + "          ------> Transitioning from modred to TS Calc\n\n")
                        	else:
                                	statusFile(file_names[i] + "          ------> Done\n\n")
                        i += 1 
	statusFile.close()
	os.chdir('modred')

def checkNegVib(inFile):
	iF = open(inFile, "r")
	line = iF.readline()
	while line:
        	if "Frequencies --" in line:
                	Freq = float(line.split()[2])
                	if Freq < 0:
                        	return float(line.split() [3]) > 0
        	line = iF.readline()
	return False	



