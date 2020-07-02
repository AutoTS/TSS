#!/bin/env python3
#small comment
#Written By Taylor Nielson

import sys, os, shutil, numpy as np, time, subprocess

inputs = {}
coords = []

#We need to add a crest_selectivity default
def default():
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
		elif "m-basis" in line:
			inputs["m-basis"] = line.split(':')[1]
		elif "m-method" in line:
			inputs["m-method"] = line.split(':')[1]
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
		elif "m-basis" in line:
			inputs["m-basis"] = line.split(':')[1]
		 #Required - Default in file
		elif "m-method" in line:
			inputs["m-method"] = line.split(':')[1]	
		#Required, default not in file
		elif "library" in line:
			inputs["library"] = line.split(':')[1]
		#Not Required, Gaussian Defaults to gas
		elif "solvent_model" in line:
			inputs["solvent_model"] = line.split(':')[1]
		 #Required - Default in file
		elif "denfit" in line:
			inputs["denfit"] = line.split(':')[1]
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
	coords = getCoords()
	return inputs, coords


def getCoords():
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
	#need to error check and to check whether the input has the desired info or if the defaults are going to be used
	#
	#if inputs["mult"] == "":
 	#	inputs["mult"] = 2
	#
	oF.write("#opt=" + inputs["opt"].strip() + " freq=noraman " + inputs["method"].strip() + "/" + inputs["basis"].strip() + " integral = ultrafine")
	oF.write("\n\n")
	oF.write("TSS")
	oF.write("\n\n")
	oF.write(inputs["charge"].strip() + " " + inputs["spin"])
	for coord in coords:
		oF.write(coord)
	oF.write("\n")
	oF.write("\n")
	oF.write("\n")
	oF.write("\n")
	oF.close()

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
	finalF = open("out.xyz", "w")
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

def gaussianProcesses():
	commands = []
	processes = []
	file_names = []
	args = sys.argv
	allDone = False
	os.chdir("modred")
	for file in os.listdir(os.getcwd()):
        	file_names.append(str(file))
        	commands.append(['/apps/gaussian16/B.01/AVX2/g16/g16', file])
	for com in commands:
        	processes.append(subprocess.Popen(com))
	while not allDone:
        	allDone = True
        	time.sleep(30)
        	i = 0
        	drawStatus(file_names,processes)
        	for p in processes:
                	if p.poll() is None:
                        	allDone = False
				#Run Check to make sure it is still TS
                	else:
				#Check for negative frequency
                        	yellow = 2
				#coords = logtoxyz()
                        	#buildcom(inputs, coords, f_name)
                        	
                        	#IF negative frequency run TS
                	i += 1
	print("all done")

def makeDirectories():
	os.mkdir("modred")
	os.mkdir("gaussianTS")

def drawStatus(file_names, processes):
	os.system('clear')
	i = 0
	for p in processes:
                        if p.poll() is None:
                                allDone = False
                                #Run Check to make sure it is still TS
                                print(file_names[i] + "          ------> Running Modred\n\n")
                        else:
                                #Check for negative frequency
                                #IF negative frequency run TS
                                print(file_names[i] + "          ------> Done\n\n")
                        i += 1 
