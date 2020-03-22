#! /usr/bin/env python3

#####################################################################
#  PrimedRPA: RPA Primer and Probe Set Finder                       #
#  Higgins M et al. Submitted. 2018                                 #
#                                                                   #
#  Dependencies:                                                    #
#     Python 3.7
#     Glob 0.6
#     Pandas 0.20.3                                                 #
#     Sys 3.6.3                                                     #
#     Bio 1.70                                                      #
#####################################################################


# Install neccessary python libraries
import os
import sys
import math
import glob
import time
import subprocess
import itertools
import datetime
import random
import numpy
import pandas as pd
from collections import Counter

# Wrapper for clustal omega
def RunningClustalo1(ClustalOPath,
					ListOfFileNames,
					overwriteOutput=True):

	for FileName in ListOfFileNames:
		OutputName = FileName.replace(".fasta",'_Aligned.fasta')
		command = "{0} -i {1} -o {2} --outfmt=fasta".format(ClustalOPath,FileName, OutputName)
		result = subprocess.call([command], stdout=subprocess.PIPE, shell=True,)


# Function to generate reverse complement sequence
def getComplement(seq,
				 reverse=False,
				 rule='N2N'):

	seqComp = ""
	for base in seq.upper():
		base = base.upper()
		if base == "A":
			seqComp += "T"
		elif base == "T":
			seqComp += "A"
		elif base == "C":
			seqComp += "G"
		elif base == "G":
			seqComp += "C"
		elif base == "N":
			if rule == 'N2-':
				seqComp += '-'
			elif rule == 'N2N':
				seqComp += "N"
		elif base == "-":
			seqComp += "N"
	if reverse:
		return(seqComp[::-1])
	else:
		return(seqComp)


def CheckingAlignedOutputFile(Targetsequencefile,
							  UserConservedDecimal,
							  DesiredPrimerLength,
							  ProbeRequired,
							  DesiredProbeLength,
							  AmpliconSize,
							  BindingSitesSeed):


	fastadict = {}
	with open(Targetsequencefile) as file_one:
		for line in file_one:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				active_sequence_name = line[1:]
				if active_sequence_name not in fastadict:
					fastadict[active_sequence_name] = ''
				continue
			sequence = line
			fastadict[active_sequence_name] += sequence.upper()


	FirstSeq = True
	FastaSeqLength = 0
	for fastaseq in fastadict.values():
		if FirstSeq == True:
			FastaSeqLength = len(fastaseq)
			FirstSeq = False
		if len(fastaseq) != FastaSeqLength:
			sys.exit('ERROR: Alignment file error length of all sequences is not equal')


	AlignedDF = pd.DataFrame()

	for seqindex in list(range(FastaSeqLength)):

		TempNucleotides = []

		for faseq in fastadict.values():
			TempNucleotides.append(faseq[seqindex])


		if (len(TempNucleotides) == 2 and
			'-' in TempNucleotides):
			MostCommonN = '-'
			NRepresentation = -100 #Harsh weighting against splits

		else:
			BugFix = Counter(TempNucleotides)
			MostCommonN = max(TempNucleotides, key=BugFix.get)

			if MostCommonN == '-':
				NRepresentation = -100 #Harsh weighting against splits
			else:
				NRepresentation = TempNucleotides.count(MostCommonN)/len(TempNucleotides)


		AlignedDF = AlignedDF.append({'Index_Pos':seqindex,'Nucleotide':MostCommonN,'Abundance':NRepresentation},ignore_index=True)

	TargetSiteLengths = [DesiredPrimerLength]
	if ProbeRequired in ['Exo','Nfo']:
		TargetSiteLengths.append(DesiredProbeLength)

	PrimerStartSites = []
	ProbeStartSites = []

	RoundOneIter = True
	for tsl in TargetSiteLengths:
		for i in list(range(FastaSeqLength-tsl)):
			MeanHomologyScore = AlignedDF.loc[i:i+tsl,'Abundance'].mean()
			if MeanHomologyScore >= UserConservedDecimal:
				if RoundOneIter == True:
					PrimerStartSites.append(i)
				else:
					ProbeStartSites.append(i)
		RoundOneIter = False

	PrimerProbeBindingSiteSets = []

	AllUnfilteredCombos = []

	for pss in PrimerStartSites:

		if len(ProbeStartSites) != 0:
			PrSSSubsets= [prss for prss in ProbeStartSites if prss >= pss+DesiredPrimerLength if prss<= pss+(AmpliconSize-DesiredProbeLength)]
			for TempPrSSS in PrSSSubsets:
				RPSSubsets = [rpss for rpss in PrimerStartSites if rpss >= TempPrSSS+DesiredProbeLength if rpss <= pss+(AmpliconSize-DesiredPrimerLength)]
				AllUnfilteredCombos += list(itertools.product([pss],[TempPrSSS],RPSSubsets))

		else:
			RPSSubsets = [rpss for rpss in PrimerStartSites if rpss >= pss+DesiredPrimerLength if rpss <= pss+(AmpliconSize-DesiredPrimerLength)]
			AllUnfilteredCombos += list(itertools.product([pss],RPSSubsets))

	print('{}: Binding Site Combinations Identified'.format(str(len(AllUnfilteredCombos))))

	#Subset AllUnfilteredCombos if necessary to maintain computational efficiency
	if BindingSitesSeed <= len(AllUnfilteredCombos):
		AllUnfilteredCombos = random.sample(AllUnfilteredCombos,BindingSitesSeed)
		print('{}: Binding Site Combinations Remaining after Subsetting'.format(str(len(AllUnfilteredCombos))))

	for ufcombo in AllUnfilteredCombos:


		# No Probe necessary
		if len(ufcombo) == 2:
			FPBindingSite = ''.join(AlignedDF.loc[ufcombo[0]:ufcombo[0]+DesiredPrimerLength,'Nucleotide'].tolist())
			RPBindingSite = ''.join(AlignedDF.loc[ufcombo[1]:ufcombo[1]+DesiredPrimerLength,'Nucleotide'].tolist())
			PPPriorSet = [FPBindingSite,RPBindingSite]

		# With probe
		else:
			FPBindingSite = ''.join(AlignedDF.loc[ufcombo[0]:ufcombo[0]+DesiredPrimerLength,'Nucleotide'].tolist())
			RPBindingSite = ''.join(AlignedDF.loc[ufcombo[2]:ufcombo[2]+DesiredPrimerLength,'Nucleotide'].tolist())
			ProbeBindingSite = ''.join(AlignedDF.loc[ufcombo[1]:ufcombo[1]+DesiredProbeLength,'Nucleotide'].tolist())
			PPPriorSet = [FPBindingSite,ProbeBindingSite,RPBindingSite]


		if '-' not in ''.join(PPPriorSet):

			PrimerProbeBindingSiteSets.append(PPPriorSet)

	if len(PrimerProbeBindingSiteSets) != 0:
		print('Associated Binding Site Sequences Successfully Extracted')
	return PrimerProbeBindingSiteSets


def RunFluroProbeAnalysis(ProbeBindingSeq,
						  minIndexPosition,
						  maxIndexPosition):

	ProbeValidPass = False
	basenumber = 0
	proberegionlength = len(ProbeBindingSeq)
	for base in ProbeBindingSeq:
		if (basenumber + 4) < proberegionlength:
			basenumber += 1
			if minIndexPosition < basenumber < maxIndexPosition and base == "T":
				if (basenumber + 2) < proberegionlength:
					if ProbeBindingSeq[(basenumber+2)] == "T" or ProbeBindingSeq[(basenumber+3)] == "T":
						# State that probe seq passed in forward sense.
						ProbeValidPass = True
						# New Break For Loop
						break
	return ProbeValidPass

def CompareSequenceForMatches(firstseq,
							  secondseq):

	MaximumScore = 0
	MaximumLength = 0
	MinimumLength = 0

	if len(secondseq) > len(firstseq):
		MaximumLength = len(secondseq)
		MinimumLength = len(firstseq)

	else:
		MaximumLength = len(firstseq)
		MinimumLength = len(secondseq)

	# Adjusted sequences
	fistseqadj = firstseq + '?' * int(MaximumLength-len(firstseq))
	secondseqadj = secondseq + '?' * int(MaximumLength-len(secondseq))

	ForwardShiftfistseqadj = fistseqadj
	ReverseShiftfistseqadj = fistseqadj

	FirstCompTest = True

	TestCounter = 0

	# Length of fistseqadj and secondseqadj should be equal
	for i in list(range(MaximumLength)):

		FSMScore = 0
		RSMScore = 0

		# Ensure we do a first test with no changes
		if TestCounter != 0:

			# Forward Shift FS
			ForwardShiftfistseqadj = '?' + ForwardShiftfistseqadj[:-1]
			# Reverse Shift FS
			ReverseShiftfistseqadj = ReverseShiftfistseqadj[1:] + '?'


		TestCounter +=1

		for it in list(range(MaximumLength)):

			if (ForwardShiftfistseqadj[it] == secondseqadj[it] and
				ForwardShiftfistseqadj[it] != '?' and
				secondseqadj[it] != '?'):

				FSMScore +=1

			if (ReverseShiftfistseqadj[it] == secondseqadj[it] and
				ReverseShiftfistseqadj[it] != '?' and
				secondseqadj[it] != '?'):

				RSMScore +=1


		if FSMScore > RSMScore:
			if FSMScore > MaximumScore:
				MaximumScore = FSMScore

		else:
			if RSMScore > MaximumScore:
				MaximumScore = RSMScore


	PercentageMatch = MaximumScore / MinimumLength

	return PercentageMatch


def BlastnBackgroundCheck(PrimerProbeSetList,
						  BlastnPath,
						  BlastnDBName,
						  RunReferenceName):


	MaxBackgroundScoreBindingScore = 0

	for seq in PrimerProbeSetList:

		#Create temp fasta file
		tempfastainput = open('./{}_Blastn_Input.fa'.format(RunReferenceName),'w')
		tempfastainput.write('>Temp_Blastn_Fasta\n{}\n'.format(seq))
		tempfastainput.close()

		#Triggure Blastn command
		blastncommandrun = '{0} -task "blastn-short"  -query  {1}_Blastn_Input.fa -db {2}/{2} -out {1}_Blastn_Output.csv -outfmt "10 sseqid pident qstart qend" '.format(BlastnPath,RunReferenceName,BlastnDBName)
		subprocess.call([blastncommandrun],shell=True)

		# Try to read dataframe (may be empty if no alignments found)
		try:

			blastnoutdf = pd.read_csv('{}_Blastn_Output.csv'.format(RunReferenceName),header=None)

			blastnoutdf['Temp_Score'] = (blastnoutdf[3] - blastnoutdf[2])*(blastnoutdf[1]/100)

			QuerySpecificPercentageId = []

			for TPSI in blastnoutdf['Temp_Score'].tolist():

				AdjustedPIScore = TPSI/len(seq)

				QuerySpecificPercentageId.append(AdjustedPIScore)


			MaximumPercentageMatch = max(QuerySpecificPercentageId) * 100


		# If dataframe empty e.g. no alignments at all
		except pd.errors.EmptyDataError:
			MaximumPercentageMatch = 0

		if MaximumPercentageMatch >= MaxBackgroundScoreBindingScore:
			MaxBackgroundScoreBindingScore = MaximumPercentageMatch

		# Delete temporary files
		os.remove('{}_Blastn_Input.fa'.format(RunReferenceName))
		os.remove('{}_Blastn_Output.csv'.format(RunReferenceName))



	return MaxBackgroundScoreBindingScore


# Basic Command Line User Interface
print('\n\n')
print('-------------------------------------------')
print('----------------PrimedRPA------------------')
print('-----Finding RPA Primer and Probe Sets-----')
print('-------------Higgins M et al.--------------')
print('-------------------------------------------\n\n')

# Basic Help Functions
if '-h' in sys.argv or '--help' in sys.argv or len(sys.argv) == 1:
	print('Usage: PrimedRPA <Path to Parameters File>\n')
	print('Example: PrimedRPA ./Validation/Validation_1_PrimedRPA_Parameters.txt\n\n')
	print('Additional information can be found at https://github.com/MatthewHiggins2017/PrimedRPA')
	print('If used, please cite:\n\tHiggins M et al. Submitted. 2018')
	sys.exit()

parametersFile = sys.argv[1]

if parametersFile in glob.glob(parametersFile):
#if os.path.isfile(parametersFile):
	paraFile = open(parametersFile,"r")
	iu = []
	for line in paraFile.readlines():
		if ">" in line:
			n = line.strip('\n')
			h = n.strip('>')
			iu.append(h)
	u = iu[1:]

	RunReferenceName = str(u[0])
	InputFileType = str(u[2])
	Targetsequencefile = str(u[1])
	UserConservedDecimal = int(u[3])/100
	DesiredPrimerLength = int(u[4])
	ProbeRequired = str(u[5])
	DesiredProbeLength = int(u[6])
	AmpliconSize = int(u[7])
	NumberOfContinRepeates = int(u[8])
	MinimumGCcontentFilter = int(u[9])
	MaximumGCcontentFilter = int(u[10])
	PrimerProbeSelfComplementaryBinding = int(u[11])
	Lengthofregionswhichcantoleratesecondarystructure = int(u[12])
	BackgroundCheckRequired = str(u[13])
	BackgroundMatchCutOffPercentage = int(u[14])
	BackgroundGenomeFile = str(u[15])
	BindingSitesSeed = int(u[16])
	PrimedSetsSeed = int(u[17])


	# Error checks contents of parameters file
	if '{}_PrimedRPA_Primer_Binding_Sites.csv'.format(RunReferenceName) in glob.glob('{}_PrimedRPA_Primer_Binding_Sites.csv'.format(RunReferenceName)):

		print('Previously Created Primer Binding Sites Identified')

	elif Targetsequencefile not in glob.glob(Targetsequencefile):
	#elif os.path.isfile(Targetsequencefile) == False:

		print("Error: Input fasta not found. Please ensure that the following path is correct: '{}'.\nPrimedRPA Terminating".format(Targetsequencefile))
		sys.exit()


	else:

		fileInfo = os.stat(Targetsequencefile)
		if fileInfo.st_size > 10000000:
			InputSizeCheck = input("We advise the input file to be less than 10MB due to delays in primer generation.\nDo you want to continue [y/n]:")
			if InputSizeCheck.lower() in ['n','no']:
				print('PrimedRPA Terminating')
				sys.exit()

	if AmpliconSize <= (DesiredPrimerLength*2 + DesiredProbeLength):
		print("Error: Defined amplicon size is less than combined primer and probe lengths.\nPrimedRPA Terminating")
		sys.exit()

	if BackgroundGenomeFile not in glob.glob(BackgroundGenomeFile) and BackgroundCheckRequired.lower() == 'yes':
		print("Error: Background fasta not found.Please ensure that the following path is correct: '{}'.\nPrimedRPA Terminating".format(BackgroundGenomeFile))
		sys.exit()


	# UPDATE: As part of Stage Two to make sure everything is reflected correctly
	print(('\n'+
		'Received Parameters:\n\n'+
		'Run Name:{RunReferenceName}\n' +
		'Input File Type: {InputFileType}\n' +
		'Input Fasta: {Targetsequencefile}\n'+
		'Homology of Conserved Target DNA: {UserConservedDecimal}%\n'+
		'Desired Primer Length: {DesiredPrimerLength:,}\n'+
		'Desired Probe Length: {DesiredProbeLength:,}\n'+
		'Desired Amplicon Length: {AmpliconSize:,}\n'+
		'Minimum GC Content: {MinimumGCcontentFilter}%\n'+
		'Maximum GC Content: {MaximumGCcontentFilter}%\n'+
		'Self-binding cut-off threshold (%): {PrimerProbeSelfComplementaryBinding:,}\n'+
		'Seconday Structure cut-off threshold (%): {Lengthofregionswhichcantoleratesecondarystructure:,}\n'+
		'Performing Background Binding Check: {BackgroundCheckRequired}\n'+
		'Background Binding cut-off threshold (%) {BackgroundMatchCutOffPercentage}\n' +
		'Background Check File: {BackgroundGenomeFile}\n').format(**globals()))
else:
	print("Error: Parameters file not found. Please ensure that the following path is correct '{}'\nPrimedRPA Terminating.".format(parametersFile))
	sys.exit()


#Define tool dependancy paths
InstalledSourceScriptPath = os.path.realpath(__file__)
BlastnPath = InstalledSourceScriptPath.replace('PrimedRPA.py','Tool_Dependancies/blastn')
BlastnDBCreationPath = InstalledSourceScriptPath.replace('PrimedRPA.py','Tool_Dependancies/makeblastdb')
ClustalOPath = InstalledSourceScriptPath.replace('PrimedRPA.py','Tool_Dependancies/clustalo')



# Check if File Housing Primer Binding Sites Already Exists For Reference Run
if '{}_PrimedRPA_Primer_Binding_Sites.csv'.format(RunReferenceName) not in glob.glob('{}_PrimedRPA_Primer_Binding_Sites.csv'.format(RunReferenceName)):
#if os.path.isfile('{}_PrimedRPA_Primer_Binding_Sites.csv'.format(RunReferenceName)) == False:

	# If necessary perform alignment on the input fasta file
	if InputFileType == 'MS':

		# Check if MS alignment already Exists in Working Directory

		if Targetsequencefile.replace('.fasta','_Aligned.fasta') not in glob.glob(Targetsequencefile.replace('.fasta','_Aligned.fasta')):
		#if os.path.isfile(Targetsequencefile.split('/')[-1].replace('.fasta','_Aligned.fasta')) == False:
			RunningClustalo1(ClustalOPath,[Targetsequencefile], overwriteOutput=False)

		Targetsequencefile = Targetsequencefile.replace('.fasta','_Aligned.fasta')


	PrimerProbeBindingSiteSets = CheckingAlignedOutputFile(Targetsequencefile,
														  UserConservedDecimal,
														  DesiredPrimerLength,
														  ProbeRequired,
														  DesiredProbeLength,
														  AmpliconSize,
														  BindingSitesSeed)

	if len(PrimerProbeBindingSiteSets) != 0:
		# Create CSV to store temp primers so the first stage does not have to be rerun.
		tempprimerdf = pd.DataFrame(PrimerProbeBindingSiteSets)

		# Provide user with option to save primers as csv
		CSVSaveInput = 'No'
		UserDefinedAnswer = False

		timeout = time.time() + 60
		tempprimerdfsize = str(sys.getsizeof(tempprimerdf)/1000000000)
		while UserDefinedAnswer == False:

			if time.time() > timeout:
				break
			else:
				CSVSaveInput = input('Save possible primer binding sites. Estimated storage space useage {} GB. [Yes/No]'.format(tempprimerdfsize))
				UserDefinedAnswer = True


		if CSVSaveInput.lower() in ['y','yes']:
			tempprimerdf.to_csv('{}_PrimedRPA_Primer_Binding_Sites.csv'.format(RunReferenceName),
																		   index=None,
																		   header=None)

# If CSV housing primer binding sites exist, extract sets accordingly.
else:
	tempdf = pd.read_csv('{}_PrimedRPA_Primer_Binding_Sites.csv'.format(RunReferenceName),
						header=None)
	PrimerProbeBindingSiteSets = tempdf.values.tolist()


# Check if background binding check is needed and if so generate the Blastn Database
if BackgroundCheckRequired.lower() in ["y","yes"]:

	BlastnDBName =  '{}_Blastn_DB_PrimedRPA'.format(BackgroundGenomeFile.split('/')[-1].split('.')[0])

	if os.path.exists(BlastnDBName) == False:

		os.makedirs(BlastnDBName)

		makeblastdbcommand = '{0} -in {1} -dbtype nucl -parse_seqids -out {2}/{2}'.format(BlastnDBCreationPath,
																					 BackgroundGenomeFile,
																					 BlastnDBName)

		print('\nCreating Blastn Database:')
		subprocess.call([makeblastdbcommand],shell=True)
		print('\n')


# Prepare for filtering
FilterSSets = []

RepeatRegionList = []
for RROI in ['N','A','G','C','T']:
	RROIString = RROI * NumberOfContinRepeates
	RepeatRegionList.append(RROIString)

GCFilteringPrimerPass = 0
RRFilteringPP = 0
ProbeFilteringPP = 0
SSandSBFilterPP = 0
BackgroundFilterPP = 0
SetCounter = 0
FilteredSetPassed  = 0
ContinueSearch = 'Yes'

MaxBackgroundBS = []

# Begin subsetting the PrimerProbe set based on 1000 seed defined seed
random.shuffle(PrimerProbeBindingSiteSets)

ListSplitsNeeded = math.ceil(len(PrimerProbeBindingSiteSets)/PrimedSetsSeed)

for PPBSSTranche in [PrimerProbeBindingSiteSets[ili::ListSplitsNeeded] for ili in range(ListSplitsNeeded)]:

	if ContinueSearch.lower() in ['no','n']:
		break

	# Begin for filtering
	for PPSet in list(PPBSSTranche):

		SetCounter += 1

		CombinedSeq = ''.join(PPSet)

		SN = 0
		for base in CombinedSeq:
			if base == "G" or base == "C":
				SN = SN + 1

		# New GC Filter + N Limit (Filter Step)
		if (SN/len(CombinedSeq) >= MinimumGCcontentFilter/100 and
			SN/len(CombinedSeq) <= MaximumGCcontentFilter/100 and
			CombinedSeq.count('N') < 4):

			# Repeat region filter
			RepeatRegionPresent = False
			for bsseq in PPSet:
				if any(rss in bsseq for rss in RepeatRegionList) == True:
					RepeatRegionPresent = True
					break

			GCFilteringPrimerPass +=1

			# Check if successfully passed repeat region filtering (Filter Step)
			if RepeatRegionPresent == False:
				RRFilteringPP +=1

				# Check if probe required
				if ProbeRequired in ['Exo','Nfo']:

					ProbeSuccess = False

					# Add in commands here to generate nfo or exo probe as required.
					ProbeBindingSeq = PPSet[1]

					if ProbeRequired == 'Exo':
						minIndexPosition = int(DesiredProbeLength*0.45)
						maxIndexPosition = int((DesiredProbeLength*0.75))


						ProbeSuccessfulOutput = RunFluroProbeAnalysis(ProbeBindingSeq,
																	  minIndexPosition,
																	  maxIndexPosition)

						# If unsuccessful try reverse complement of the probe
						if ProbeSuccessfulOutput == False:
							RCProbeSuccessfulOutput = RunFluroProbeAnalysis(getComplement(ProbeBindingSeq,reverse=True),
																			minIndexPosition,
																			maxIndexPosition)

							# If probe successful when reverse complementing
							if RCProbeSuccessfulOutput == True:

								# Update the PPSet variable with reverse complement of probe sequence
								PPSet[1] = getComplement(ProbeBindingSeq)

								ProbeSuccess = True

						# If successful as forward sense (first attempt)
						else:
							ProbeSuccess = True

					if ProbeRequired =='Nfo':

						# Add in necessary check steps / Antibody attachment

						ProbeSuccess = True


				# Check if probe not required or if required it is successful (Filter Step)
				if ProbeRequired.lower() in ['no','n'] or ProbeSuccess == True:

					ProbeFilteringPP +=1

					# Reverse complement the reverse primer to obtain actual primer sequence.
					PPSet[-1] = getComplement(PPSet[-1],reverse=True)

					# Filter sets which are not self-binding. ### TO DO AND IMPROVE ###
					SelfBindingFilterPass = False

					MaximumSBPercentageMatch = 0
					MaximumSSPercentageMatch = 0


					# Loop through index of seqs in set
					for porpi in list(range(len(PPSet))):

						PorpiTempPMCValue = CompareSequenceForMatches(PPSet[porpi],getComplement(PPSet[porpi]))

						# Get list of other sequences in set.
						OtherSequences = [oss for oss in PPSet if oss != PPSet[porpi]]

						# Loop through other sequences
						for osst in OtherSequences:

							ReverseTempPMCValue = CompareSequenceForMatches(PPSet[porpi],getComplement(osst))
							ReverseCompTempPMCValue = CompareSequenceForMatches(PPSet[porpi],getComplement(osst,reverse=True))

							if MaximumSBPercentageMatch < max([PorpiTempPMCValue,ReverseTempPMCValue,ReverseCompTempPMCValue]):
								MaximumSBPercentageMatch = max([PorpiTempPMCValue,ReverseTempPMCValue,ReverseCompTempPMCValue])

						#Check Seq for SS ability:
						SSPTempValue = CompareSequenceForMatches(PPSet[porpi],getComplement(PPSet[porpi],reverse=True))

						if  MaximumSSPercentageMatch < SSPTempValue:
							MaximumSSPercentageMatch = SSPTempValue


					# Filter on percentage of matches to self-bind between members of primer set or form secondary structure.
					if (MaximumSBPercentageMatch <= PrimerProbeSelfComplementaryBinding/100 and
						MaximumSSPercentageMatch <= Lengthofregionswhichcantoleratesecondarystructure/100):

						SSandSBFilterPP += 1

						# Filter on background binding check if it is required.
						if BackgroundCheckRequired.lower() in ['y',"yes"]:

							BackgroundBindingMaxScore = BlastnBackgroundCheck(PPSet,BlastnPath,BlastnDBName,RunReferenceName)
							if BackgroundBindingMaxScore <= BackgroundMatchCutOffPercentage:

								BackgroundFilterPP += 1

								# If set passess the filter append to final list
								FilterSSets.append(PPSet)
								FilteredSetPassed +=1

								MaxBackgroundBS.append(BackgroundBindingMaxScore)

						# If no background binding check is required:
						else:

							# If set passess all filtering stages add set to list.
							FilterSSets.append(PPSet)

							FilteredSetPassed +=1

		print('Sets Total : Sets Assessed : Sets Passed: [{}:{}:{}]\r'.format(len(PrimerProbeBindingSiteSets),SetCounter,FilteredSetPassed), end="")

	# Ask if user would like to continue if primers found but no all possibilities explored.
	if (FilteredSetPassed != 0 and
		SetCounter != len(PrimerProbeBindingSiteSets)):

		ContinueSearch = input('{} Primer Sets Found, would you like to continue searching for more? [Yes/No] : '.format(FilteredSetPassed))
		if ContinueSearch.lower() in ['n','no']:
			break

# Print output of filtering stages.
print('\n\nPrimers Passed GCFiltering: [{}]\nPrimers Passed Repeat Region Filtering [{}]\nPrimers Passed Probe Filtering [{}]\nPrimers Passed Seconday Structure and Self-Binding Filtering [{}]'.format(GCFilteringPrimerPass,RRFilteringPP,ProbeFilteringPP,SSandSBFilterPP))


if BackgroundCheckRequired.lower() == ['y','yes']:
	print('Primers Passed Background Filtering [{}]\n'.format(BackgroundFilterPP))


# Report if no primers were found
if len(FilterSSets) == 0:
	print('No primers passed the filtering stage')

	sys.exit()


# Export primers if filtering successful
else:

	ColumnsNamesTemp = ['Forward_Primer','{}_Probe'.format(ProbeRequired),'Reverse_Primer']
	if ProbeRequired.lower() in ['no','n']:
		ColumnsNamesTemp.remove('{}_Probe'.format(ProbeRequired))

	OutputDataFrame = pd.DataFrame(FilterSSets,columns=ColumnsNamesTemp)

	# Add necessary max Blastn score
	if len(MaxBackgroundBS) != 0:
		OutputDataFrame['Max_Blastn_Percentage_Match'] = MaxBackgroundBS

	OutputDataFrame.to_csv('{}_PrimedRPA_Output.csv'.format(RunReferenceName),index=None)
	print('\nPrimedRPA Run {} Complete!'.format(RunReferenceName))
	print('\nIf you found this program useful please reference:\n\nM. Higgins et al., “PrimedRPA: primer design for recombinase polymerase amplification assays,”\nBioinformatics, vol. 35, no. 4, pp. 682–684, Feb. 2019.\n\n')
	sys.exit()
