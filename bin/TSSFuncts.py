#!/bin/env python3
#Written By Taylor Nielson

import sys, os,re, shutil, numpy as np, time, subprocess


#This list is currently incomplete But it is a mapping of metals and nonmetals to their atomic number
metals_lib = {"26":"FE","51":"SB","75":"RE","3":"LI","4":"BE","11":"NA","12":"MG","78":"PT"}
non_metals_lib = {"7":"N","15":"P","1":"H","6":"C","16":"S","8":"O","9":"F","10":"NE","17":"CL","18":"AR","34":"SE","35":"BR","36":"KR","53":"I","54":"XE","86":"RN"}

#Loads defaults from the default file into the inputs variable
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


#Parses the input file and overrides the defaults with values found in the input file
def parseInput(inputFile, inputs,xyz_file):
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
	#	if "spin" in line:
	#		inputs["spin"] = line.split(':')[1]
		 #Required - Default in file
	#	elif "charge" in line:
	#		inputs["charge"] = line.split(':')[1]
		 #Required - Default in file
		if "basis" in line:
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
#		elif "batch" in line:
#			inputs["batch"] = line.split(':')[1]
		 #Required - Default in file
		elif "metal" in line:
			inputs["mbasis"] = line.split(':')[1]
		#Required, default not in file
#		elif "library" in line:
#			inputs["library"] = line.split(':')[1]
		#Not Required, Gaussian Defaults to gas
		elif "solvent_model" in line:
			inputs["solvent_model"] = line.split(':')[1]
		 #Required - Default in file
		elif "denfit" in line:
                	inputs["denfit"] = line.split(':')[1]
		elif "memory" in line:
			inputs["mem"] = line.split(':')[1]
		elif "memper" in line:
			inputs["memper"] = line.split(':')[1]
		elif "confs" in line:
			inputs["confs"] = line.split(':')[1]
		elif "extra" in line:
			inputs["extra"] = line.split(':')[1]
#		elif "time" in line:
#			inputs["time"] = line.split(':')[1]
#		#Not Required
#		elif "subtract" in line:
#			parseChanges(line, "subtract", inputs)
#		#Not Required
#		elif "substitute" in line:
#			parseChanges(line, "substitute", inputs)
#		#Not Required
#		elif "add" in line:
#			parseChanges(line, "add", inputs)
		elif "difficulty" in line:
			inputs["difficulty"] = line.split(':')[1]
		elif "conformational_leniency" in line:
			inputs["leniency"] = line.split(':')[1]
		line = iF.readline()
	iF.close()
	#coords = getCoords(inputs,xyz_file)
	return inputs# coords

#Gets the coordinates from the xyz file and puts them into a list
def getCoords(inputs,xyz_file):
	coords = []
	#iF = open(os.path.expanduser("~/TSS/libs/base_templates/" + inputs["library"].strip()), "r")
	iF = open(xyz_file,"r")
	line = iF.readline()
	line = iF.readline()
	parseFreezes(line)
	line = iF.readline()
	while line:
		coords.append(line)
		line = iF.readline()
	return coords
	
		
#Builds a com file using values from input file and defaults for any not specified
def buildCom(inputs, coords, f_name):
	oF = open(f_name, 'w')
	oF.write("%nprocshared=1\n")
	oF.write("%mem="+ str(int(inputs["memper"])//int(inputs["confs"])) + "GB\n")
	oF.write("#opt=" + inputs["opt"].strip() + " freq=noraman " + inputs["method"].strip() + "/" + "genecp"  + " integral=ultrafine" + " " + inputs["extra"].strip())
	oF.write("\n\n")
	oF.write("TSS")
	oF.write("\n\n")
	oF.write(inputs["charge"].strip() + " " + inputs["spin"] + "\n")
	for coord in coords:
        	coord = coord.strip()
        	if coord[0].isalpha():
                        oF.write(coord + "\n")
        	else:
        		coord = coord.split()
        		if coord[0] in metals_lib:
        			string = " "
        			coord[0] = metals_lib[str(coord[0])]
        			oF.write(string.join(coord)+ "\n")
        		else:
        			string = " "
        			coord[0] = non_metals_lib[str(coord[0])]
        			oF.write(string.join(coord)+ "\n")
	writeFreezes(oF, coords, inputs)
	oF.write("\n")
	oF.write("\n")
	oF.write("\n")
	oF.write("\n")
	oF.close()

#Writes the freezes in the modred file
def writeFreezes(outFile,coords, inputs):
	#Freezes = []
	#libFile = open(os.path.expanduser("~/TSS/libs/base_templates/" + inputs["library"].strip()), "r")
	#libFile.readline()
	#freeze_line = libFile.readline()
	#freeze_line = freeze_line[2:]
	#freeze_line = freeze_line.split(';')
	#freeze_line.pop(-1)
	#if "modred" in inputs["opt"]:
        #	for i in range(0, len(freeze_line)):
        #        	val = freeze_line[i].split('-')
        #        	Freezes.append(val[0])
        #        	Freezes.append(val[1])
        #	outFile.write("\n")
	#for i in range(0,len(Freezes)-1,2):
        #	outFile.write("B " + str(int(Freezes[i]) + 1) + " " + str(int(Freezes[i+1]) + 1) + " F\n")
	outFile.write("\n")
	if "modred" in inputs["opt"]:
		for bond in inputs["bonds"]:
			outFile.write("B " + str(bond[0]) + " " + str(bond[1]) + " F\n")
	writeGenecp(outFile, coords, inputs)	


#Writes the genecp section in the modred and TS search file
def writeGenecp(outFile,coords, inputs):
	metals, non_metals = getAtomTypes(coords)
	if "modred" in inputs["opt"]:
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
		


#Gets a list of what types of atoms are in the file to be used in the write genecp section
def getAtomTypes(coords):
	metals = set()
	non_metals = set() 
	for coord in coords:
        	coord = coord.strip()
        	if coord.split()[0].upper() in metals_lib.values():
                	metals.add(coord.split()[0])
        	elif coord.split()[0] in metals_lib:
                        metals.add(metals_lib[str(coord.split()[0])])
        	else:
                	if coord[0].isalpha():
                        	non_metals.add(coord.split()[0])
                	else:
                        	non_metals.add(non_metals_lib[str(coord.split()[0])])
	return metals, non_metals
	

#Copies the base_input as defined in the input file and modifies it to build a new xyz file
#Currently only does subractions
#Currently not being used
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

#Uses the results of the modred to make com files for the modred section
def modredCrest(crest_file, inputs):
	energies = []
	acceptable_energy = True
	os.chdir("modred")
	coords_list = []
	names_list = []
	num_structures = 0
	iF = open("../"+crest_file,"r")
	line = iF.readline()
	while line:
		line = iF.readline()
		energies.append(float(line)*627.51)
		line = iF.readline()
		coords = []
		while line:
			if line.strip()[0].isalpha():
				coords.append(line)
				line = iF.readline()
			else:
				break
		#for i in range(0,len(energies)-1):
		#	if energies[-1] -energies[i] < 0.10:
		#		acceptable_energy = False
		#		energies.pop()
		#		break
		if acceptable_energy:
			coords_list.append(coords)
			names_list.append("conf" + str(num_structures) + ".com")
			num_structures +=1
		acceptable_energy = True
	inputs["numconfs"] = num_structures
	#for i in range(0,len(coords_list)):
#		buildCom(inputs, coords_list[i], names_list[i])
	try:
		for i in range(len(coords_list)-1, len(coords_list)-int(inputs["confs"])-1,-1):
			buildCom(inputs, coords_list[i], names_list[i])
	except:
		print("okay\n")
	os.chdir("../")
# def modredRangeCreation():
	#This will take the current modreds and create multiple ones with differing frozen bond lengths


#Takes the log of the modred and gets the resulting xyz coords
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
	z = 0
	return coords


#Runs the gaussian jobs and monitors them. If the modred finishes but doesn't have an imaginary vibration it is killed otherwise it runs the transition state search
def gaussianProcesses(inputs):
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
        	time.sleep(300)
        	i = 0
        	drawStatus(file_names, processes, optType, switched)
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
                                        	inputs["numconfs"] = findAliveProcesses(processes)
                                        	processes[i] = subprocess.Popen(['/apps/gaussian16/B.01/AVX2/g16/g16', (file_names[i] + ".com")])
                                        	optType[i] = "TS Calc"
                                       		os.chdir('../modred')
                                        	switched[i] = 1
                                	else:
                                        	optType[i] = "killed"
                	i += 1
	drawStatus(file_names, processes, optType, switched)
	os.chdir("../gaussianTS")
	files = [f for f in os.listdir('.') if f.split('.')[1] is "log"]
	for f in files:
        	CheckPassFail(f)
	print("all done")

#Helper function to get the number of processes still alive
def findAliveProcesses(processes):
	num = 0
	for p in processes:
		if p.poll() is None:
			num +=1
	if num == 1:
		return 2
	return num

#Helper function to make all the directories
def makeDirectories():
	os.mkdir("modred")
	os.mkdir("gaussianTS")
	os.mkdir("crest")
	os.mkdir("completed")

#Helper function to make and update the status file to see the status of the jobs
def drawStatus(file_names, processes, optType, switched):
	os.chdir('../')
	statusFile = open("status", "w")
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
                                	statusFile.write(file_names[i] + "          ------> Done\n\n")
                        i += 1 
	statusFile.close()
	os.chdir('modred')

#Helper function to check the files for imaginary vibrations
def checkNegVib(inFile):
	iF = open(inFile, "r")
	line = iF.readline()
	while line:
        	if "Frequencies --" in line:
                	Freq = float(line.split()[2])
                	if Freq < 0:
                        	return float(line.split() [3]) > 0
        	line = iF.readline()
	iF.close()
	return False
#Helper function to check to see if the file finished without errors	
def checkCompleted(inFile):
	iF = open(inFile, "r")
	for line in iF:
		pass
	last_line = line
	iF.close()
	if "Normal termination" in last_line:
		return true
	return false
	
#Helper function to check to see if the file has imaginary vibration and completed normally
def CheckPassFail(inputFile):
	pass1 = checkNegVib(inputFile)
	if pass1:
		pass2 = checkCompleted(inputFile)
		if pass2:
			copyfile(inputFile, "../completed/" + inputFile)


#Helper function to make the output for each file
def outputFunc(f_name):
	thval = ""
	iF = open(f_name, 'r')
	line = iF.readline()
	while line:
		if "Zero-point correction" in line:
			no_of_lines = 7
			lines = line
			for i in range(no_of_lines):
				lines+=iF.readline()
			thval = lines
		line = iF.readline()
	iF.close()
	return thval
        
#Creates the completed output file
def finalOutput():
	os.chdir("../completed")
	divider = "-" * 50 + '\n'
	results = []
	oF = open("autots.out", 'w')
	oF.write("\nAuto TS Output File\n")
	for file in os.listdir("."):
		results.append(outputFunc(file))
	for result in results:
		oF.write(divider)
		oF.write(result)
	oF.close()

def runCrest(xyz_file, leniency, inputs):
	coords = []
	header = ''
	bonds = []
	libFile = open(xyz_file, "r")
	header = libFile.readline()
	bonds_line = libFile.readline()[2:]
	bond_strings = bonds_line.split(';')
	bond_strings.pop()
	charge_multiplicity = bond_strings.pop()
	charge_multiplicity = charge_multiplicity[2:]
	charge_multiplicity = charge_multiplicity.split(',')
	inputs["charge"]=charge_multiplicity[0]
	inputs["spin"]=charge_multiplicity[1]
	#with open(xyz_file, "r") as coorid_file:
	header += "\n"
	bond_strings = bond_strings[0].split(',')
	if bond_strings[0] == '':
		peanut = 5
	else:
		for bond in bond_strings:
			atoms = bond.split('-')
			bonds.append([str(int(atoms[0]) + 1), str(int(atoms[1]) + 1)])
	coords = libFile.readlines()
	inputs["bonds"] = bonds
	os.chdir("crest")
	with open("cinp", "w") as constraint_file:
		constraint_file.write("$constrain\n")
		for atoms in bonds:
			constraint_file.write("\tdistance: " + atoms[0] + ", " + atoms[1] + ", auto\n")
		constraint_file.write("$end\n")
	with open("coords.xyz", "w") as crest_coords:
		crest_coords.write(header)
		crest_coords.writelines(coords)
	run_crest = subprocess.Popen([os.path.expanduser("~/crest"), "coords.xyz", "-cinp", "cinp", "-ewin", leniency])
	run_crest.wait()##				    ^^^ replace later with global path to crest.exe
	shutil.copy("crest_conformers.xyz", "../crest_conformers.xyz")
	os.chdir("../")
	tempFile = open("crest_conformers.xyz","a")
	origFile = open(xyz_file,"r")
	line = origFile.readline()
	while line:
		if "F:" in line:
			tempFile.write("0\n")
		else:
			tempFile.write(line)
		line = origFile.readline()
	tempFile.close()
	origFile.close()
