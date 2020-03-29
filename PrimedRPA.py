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
from multiprocessing import Pool
import argparse

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



## Background Binding Check = TO IMPROVE
def BlastnBackgroundCheck(seq,AllParameter):


	MaxBackgroundScoreBindingScore = 0
	MaxScoreBackSeq = ''


	#Create temp fasta file
	tempfastainput = open('./{}/{}_Blastn_Input.fa'.format(AllParameter.PrimerBlastnOutput,seq),'w')
	tempfastainput.write('>Temp_Blastn_Fasta\n{}\n'.format(seq))
	tempfastainput.close()

	#Triggure Blastn command
	blastncommandrun = '{0} -task "blastn-short"  -query  {3}/{1}_Blastn_Input.fa -db {2}/{2} -out {3}/{1}_Blastn_Output.csv -outfmt "10 sseqid pident qstart qend sstart send evalue gaps" '.format(AllParameter.BlastnPath,
																																											 				   seq,
																																											 			  	   AllParameter.BlastnDBName,
																																											 			  	   AllParameter.PrimerBlastnOutput)
	subprocess.call([blastncommandrun],shell=True)

	# Try to read dataframe (may be empty if no alignments found)
	try:

		blastnoutdf = pd.read_csv('{}/{}_Blastn_Output.csv'.format(AllParameter.PrimerBlastnOutput,seq),header=None)
		blastnoutdf.columns=['Background_SeqID','Percentage Identity','qStart','qEnd','BackSeq_Start','BackSeq_End','Evalue','Number Of Gaps']

		blastnoutdf['Cross Reactivity Score'] = ((((blastnoutdf['qEnd']+1) - blastnoutdf['qStart'])*(blastnoutdf['Percentage Identity']/100))/len(seq))*100

		blastnoutdf = blastnoutdf.sort_values(by=['Cross Reactivity Score'],ascending=False)


		blastnoutdf.to_csv('{}/{}_Blastn_Output.csv'.format(AllParameter.PrimerBlastnOutput,seq),index=None)

		IndexOfInterest = blastnoutdf.index.tolist()[0]

		MaximumPercentageMatch = blastnoutdf.loc[IndexOfInterest,'Cross Reactivity Score']
		MaxHomologyBackgroundSeq = blastnoutdf.iloc[IndexOfInterest,0]

	# If dataframe empty e.g. no alignments at all
	except pd.errors.EmptyDataError:
		MaximumPercentageMatch = 0

	if MaximumPercentageMatch > MaxBackgroundScoreBindingScore:
		MaxBackgroundScoreBindingScore = MaximumPercentageMatch
		MaxScoreBackSeq = MaxHomologyBackgroundSeq

	return (MaxBackgroundScoreBindingScore, MaxScoreBackSeq)


# Creates Exo Fluorescent Probe
def RunFluroProbeAnalysis(ProbeBindingSeq):

	minIndexPosition = int(len(ProbeBindingSeq)*0.45)
	maxIndexPosition = int(len(ProbeBindingSeq)*0.75)

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


# Secondary Structure Filter
def SSIdentification(SeqOne, SeqTwo, threshold=4):

	# Nucleotide Dict
	NucDict = {'A':'T',
			   'T':'A',
			   'C':'G',
			   'G':'C'
				}

	# Maximum Binding sites
	MaxBindingSites = 0
	MaxBindingString = ''
	MaxBindingPercentage = 0

	# Max length of string according to shift.
	MaxLength = len(SeqOne) + len(SeqTwo) - 1

	# Add white space to end of sequences
	SeqOneNorm = SeqOne + ' '*(MaxLength-len(SeqOne))
	SeqTwoNorm = SeqTwo + ' '*(MaxLength-len(SeqTwo))
	SeqTwoReversed = SeqTwo[::-1] + ' '*(MaxLength-len(SeqTwo))

	# Loop through var pairs whereby the first element will be shifted
	for VarPair in [(SeqTwoNorm,SeqOneNorm),
					(SeqOneNorm,SeqTwoNorm),
					(SeqOneNorm,SeqTwoReversed),
					(SeqTwoReversed,SeqOneNorm)]:

		# Loop through possible shift positions
		for i in list(range(VarPair[0].count(' ')+1)):

			if i == 0:
				DynamicSeq = VarPair[0]
			else:
				DynamicSeq = ' '*i + VarPair[0][:-i]

			DyanimicSeqList = list(DynamicSeq)
			FixedSeqList = list(VarPair[1])

			# This shall house the syntax to represent binding
			PossibleBindingString = ''

			for StringPos in list(range(len(DyanimicSeqList))):
				if DyanimicSeqList[StringPos] in list(NucDict.keys()):
					if NucDict[DyanimicSeqList[StringPos]] == FixedSeqList[StringPos]:
						PossibleBindingString += '|'
					else:
						PossibleBindingString += '-'
				else:
					PossibleBindingString += '-'

			NumberOfBindingMatches = PossibleBindingString.count('|')
			CompleteBindingString = DynamicSeq + '\n' + PossibleBindingString + '\n' + VarPair[1]


			if MaxBindingSites < NumberOfBindingMatches:
				MaxBindingSites = NumberOfBindingMatches
				MaxBindingString = CompleteBindingString
				MaxBindingPercentage = (MaxBindingSites/min([len(SeqOne),len(SeqTwo)]))*100


	#SSPossibilitiesDF = SSPossibilitiesDF.sort_values(by=['Complementary_Bases'], ascending=False)
	return (MaxBindingPercentage,MaxBindingString)


# Creating Alignment DF Multithread Function
def CreatingInputHomologyDF(fastadict,FastaIndexList):
	AlignedDF = pd.DataFrame()

	for seqindex in FastaIndexList:
		TempNucleotides = []
		for faseq in fastadict.values():
			TempNucleotides.append(faseq[seqindex])

		# Remove NoN-specific Nucleotides = IUPAC designated the symbols for nucleotides
		TempNucleotides = [TN for TN in TempNucleotides if TN not in ['W','S','M',
																	  'K','R','Y',
																      'N']]
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


	return AlignedDF


def IndentifyingAndFilteringOligos(AllParameter,
								   AlignedDF,
								   PossibleStartIndexes):


	OligoDF = pd.DataFrame()
	TargetSiteLengths = [AllParameter.PrimerLength]
	if AllParameter.ProbeRequired in ['EXO','NFO']:
		TargetSiteLengths.append(AllParameter.ProbeLength)

	# Loop through primer or probe
	for TSLP in TargetSiteLengths:
		if TSLP == AllParameter.PrimerLength:
			OligoType = 'Primer'
		else:
			OligoType = 'Probe'

		# Add in ability to handle specified range or primer or probe values.
		if '-' in TSLP:
			TSLRangePrior = [int(TSLI) for TSLI in TSLP.split('-')]
			TSLList = list(range(TSLRangePrior[0],TSLRangePrior[1]+1))

		else:
			TSLList = [int(TSLP)]

		# Loop through possible oligo lengths.
		for TSL in TSLList:

			# Loop through possible start sites
			for i in PossibleStartIndexes:

				# Make Sure dont go too far out of dataframe
				if i+TSL<len(AlignedDF)-TSL:

					DFSubSet = AlignedDF[(AlignedDF['Index_Pos']>=i)&(AlignedDF['Index_Pos']<=(i+TSL-1))]
					MeanHomologyScore = DFSubSet['Abundance'].mean()


					if MeanHomologyScore > AllParameter.IdentityThreshold:
						NucleotideSeq = ''.join(DFSubSet['Nucleotide'].tolist())
						NucleotideSeq = NucleotideSeq.upper()

						#Check Length of Oligo
						if len(NucleotideSeq) == TSL:

							# Assess GC Content
							GCPercentage = ((NucleotideSeq.count('G') + NucleotideSeq.count('C'))/len(NucleotideSeq))*100
							if (GCPercentage < AllParameter.MaxGC and GCPercentage > AllParameter.MinGC):

								# Assess Repeat Nucleotide
								NoRepeat = True
								for RROI in ['N','A','G','C','T']:
									RROIString = RROI * AllParameter.NucleotideRepeatLimit
									if RROIString in NucleotideSeq:
										NoRepeat = False

								if NoRepeat == True:

									# Run Secondary Structure Check / Self Dimerisation
									SDSSFilterPass=True
									MaxBindingSites, MaxBindingString = SSIdentification(NucleotideSeq,NucleotideSeq)
									if MaxBindingSites > AllParameter.DimerisationThresh:
										SDSSFilterPass=False


									if SDSSFilterPass == True:
										# Define row dictionary if passed all above stages
										RowDict = {'Oligo_Sequence':NucleotideSeq,
												   'Oligo_Type':OligoType,
												   'Binding_Site_Start_Index':i,
												   'Oligo_Length':TSL,
												   'Conservation_Score':MeanHomologyScore,
												   'GC_Content': GCPercentage,
												   'Self Dimerisation / Secondary Structure Percentage':MaxBindingSites,
												   'Self Dimerisation / Secondary Structure String':MaxBindingString}


										# Assess if Specific Probe Check is Needed.
										ProbePass = True
										if OligoType == 'Probe':

											# Run Specific Exo Probe Checks
											if AllParameter.ProbeRequired == 'EXO':
												ProbePass = RunFluroProbeAnalysis(NucleotideSeq)



										if ProbePass == True:

											# Add Run Blastn Check
											BlastnPass = True
											#If Background Check Neccessary Look To Change Blastn pass
											if AllParameter.BackgroundCheck.upper() != 'NO':

												# Run Blastn Check And If Passess Write Out Set
												MaxBackgroundScoreBindingScore, MaxScoreBackSeq  = BlastnBackgroundCheck(NucleotideSeq, AllParameter)


												# Check to see if Max Binding Score Less Than Threshold
												if MaxBackgroundScoreBindingScore > AllParameter.CrossReactivityThresh:
													BlastnPass = False
													# Remove Background If Neccessary
													os.remove('{}/{}_Blastn_Input.fa'.format(AllParameter.PrimerBlastnOutput, NucleotideSeq))
													os.remove('{}/{}_Blastn_Output.csv'.format(AllParameter.PrimerBlastnOutput, NucleotideSeq))

												else:
													RowDict['Max Background Binding Score'] = MaxBackgroundScoreBindingScore
													RowDict['Max Background Binding SeqID'] = MaxScoreBackSeq


											# If it passed Background Filtering
											if BlastnPass == True:
												OligoDF = OligoDF.append(RowDict,ignore_index=True)

	# Return Primer/Probe DataFrame
	return OligoDF


def ComboIdentifyier(PrimerSS,ReversePrimerSS,ProbeSS,AllParameter,PassedOligos,FPrimerL,RPrimerL,ProbeL):

	PassedSetsDataFrame = pd.DataFrame()
	CombosList = []

	for FP in PrimerSS:

		# If Probe Required
		if AllParameter.ProbeRequired in ['EXO','NFO']:
			Probes = [prss for prss in ProbeSS if (prss >= (FP+FPrimerL) and
												   prss<= (FP+AllParameter.AmpliconSizeLimit-RPrimerL-ProbeL))]

			for Probe in Probes:
				ReversePrimers = [rpss for rpss in ReversePrimerSS if (rpss >= (Probe+ProbeL) and
																rpss <= (FP+AllParameter.AmpliconSizeLimit-RPrimerL))]

				CombosList += list(itertools.product([FP],ReversePrimers,[Probe]))

		# Probe Not Required
		else:
			ReversePrimers = [rpss for rpss in ReversePrimerSS if (rpss >= (FP+FPrimerL) and
															rpss <= (FP+AllParameter.AmpliconSizeLimit-RPrimerL))]

			CombosList += list(itertools.product([FP],ReversePrimers))


	# Extract Valid Primer-Probe Combos
	for Combo in CombosList:
		# Identify Indexes
		FPIndex = PassedOligos[(PassedOligos['Oligo_Type']=='Primer')&(PassedOligos['Oligo_Length']==FPrimerL)&(PassedOligos['Binding_Site_Start_Index']==Combo[0])].index.tolist()[0]
		RPIndex = PassedOligos[(PassedOligos['Oligo_Type']=='Primer')&(PassedOligos['Oligo_Length']==RPrimerL)&(PassedOligos['Binding_Site_Start_Index']==Combo[1])].index.tolist()[0]


		# Extract Combo Primer Sequences - (Reverse Complement RP)
		ComboSeqList = [PassedOligos.loc[FPIndex,'Oligo_Sequence'],getComplement(PassedOligos.loc[RPIndex,'Oligo_Sequence'],True)]

		# Add Probe Sequence If Necessary
		if AllParameter.ProbeRequired in ['EXO','NFO']:
			ProbeIndex = PassedOligos[(PassedOligos['Oligo_Type']=='Probe')&(PassedOligos['Oligo_Length']==ProbeL)&(PassedOligos['Binding_Site_Start_Index']==Combo[2])].index.tolist()[0]
			ComboSeqList.append(PassedOligos.loc[ProbeIndex,'Oligo_Sequence'])


		#Run Dimerisation Check And If Passess Continue
		DimerisationPass = True
		MaxComboSSScore = 0
		MaxComboSSString = ''
		for SSSCombo in list(itertools.combinations(ComboSeqList,r=2)):
			SSMaxBindingSites, SSMaxBindingString = SSIdentification(SSSCombo[0],SSSCombo[1])

			if SSMaxBindingSites > MaxComboSSScore:
				MaxComboSSScore = SSMaxBindingSites
				MaxComboSSString = SSMaxBindingString

		if MaxComboSSScore > AllParameter.DimerisationThresh:
			DimerisationPass = False

		if DimerisationPass == True:

			BlastnPass = True ### REMOVE INDEX AT LATER DATE


			# Add Everything To DF If Passed All Parameters
			if BlastnPass == True:

				PassedComboRow = {'Forward Primer (FP)':ComboSeqList[0],
								  'FP GC%':  PassedOligos.loc[FPIndex,'GC_Content'],
								  'FP Binding Start Site':  PassedOligos.loc[FPIndex,'Binding_Site_Start_Index'],
								  'Reverse Primer (RP)': ComboSeqList[1],
								  'RP GC%': PassedOligos.loc[RPIndex,'GC_Content'],
								  'RP Binding Start Site': PassedOligos.loc[RPIndex,'Binding_Site_Start_Index'],
								  'Amplicon Size': PassedOligos.loc[RPIndex,'Binding_Site_Start_Index'] + RPrimerL - PassedOligos.loc[FPIndex,'Binding_Site_Start_Index'],
								  'Max Secondary Structure / Dimerisation Percentage Score':SSMaxBindingSites,
								  'Max Secondary Structure / Dimerisation String':MaxComboSSString,
								  'Forward Primer Length':FPrimerL,
								  'Reverse Primer Length':RPrimerL
								  }

				MaxBackgroundScoresIndexes = [FPIndex,RPIndex]

				if AllParameter.ProbeRequired in ['EXO','NFO']:
					PassedComboRow['Probe (P)'] = ComboSeqList[2]
					PassedComboRow['Probe GC%'] = PassedOligos.loc[ProbeIndex,'GC_Content']
					PassedComboRow['Probe Binding Start Site'] = PassedOligos.loc[ProbeIndex,'Binding_Site_Start_Index']
					PassedComboRow['Probe Length'] = ProbeL
					MaxBackgroundScoresIndexes.append(ProbeIndex)


				if AllParameter.BackgroundCheck != 'NO':

					MaxBackgroundScoreBindingScore = 0
					MaxScoreBackSeq = ''
					for SucIndex in MaxBackgroundScoresIndexes:
						if PassedOligos.loc[SucIndex,'Max Background Binding Score'] > MaxBackgroundScoreBindingScore:
							MaxBackgroundScoreBindingScore = PassedOligos.loc[SucIndex,'Max Background Binding Score']
							MaxScoreBackSeq = PassedOligos.loc[SucIndex,'Max Background Binding SeqID']


					PassedComboRow['Max Background Binding Score'] = MaxBackgroundScoreBindingScore
					PassedComboRow['Max Background Binding Seq'] = MaxScoreBackSeq


				# Add Passed Row To DataFrame
				PassedSetsDataFrame = PassedSetsDataFrame.append(PassedComboRow,ignore_index=True)

	return PassedSetsDataFrame

# Main Function To Generate Primer - Probe Combos
def CheckingAlignedOutputFile(AllParameter):


	# Check If Binding Sites Already Exist
	if AllParameter.PriorBindingSite != 'NO':

		print('Reading In Oligo Binding Sites')
		PassedOligos = pd.read_csv('{}'.format(AllParameter.PriorBindingSite))

		# Filter on Parameters Passed In
		PassedOligos = PassedOligos[(PassedOligos['GC_Content']>=AllParameter.MinGC)&
									(PassedOligos['GC_Content']<=AllParameter.MaxGC)&
									(PassedOligos['Self Dimerisation / Secondary Structure Percentage']<=AllParameter.DimerisationThresh)&
									(PassedOligos['Conservation_Score']>=AllParameter.IdentityThreshold)]


		# Filter on Primer + Probe Ranges - (Tidy Up At Later Date)
		PPLengthDict = {}
		PrimerRange=False
		ProbeRange = False
		for LL in [AllParameter.PrimerLength,AllParameter.ProbeLength]:
			if LL == AllParameter.PrimerLength:
				LLType = 'Primer'
			else:
				LLType = 'Probe'

			if '-' in LL:
				LLRangePrior = [int(LLLI.strip()) for LLLI in LL.split('-')]
				PPLengthDict['{}_Max'.format(LLType)] = max(LLRangePrior)
				PPLengthDict['{}_Min'.format(LLType)] = min(LLRangePrior)
				if LLType == 'Primer':
					PrimerRange = True
				else:
					ProbeRange = True
			else:
				PPLengthDict['{}_Len'.format(LLType)] =  int(LL)


		if PrimerRange == True:
			PrimerPassedSubet = PassedOligos[(PassedOligos['Oligo_Type']=='Primer') &
											 (PassedOligos['Oligo_Length']>=PPLengthDict['Primer_Min']) &
											 (PassedOligos['Oligo_Length']<=PPLengthDict['Primer_Max']) ]
		else:
			PrimerPassedSubet = PassedOligos[(PassedOligos['Oligo_Type']=='Primer') &
											 (PassedOligos['Oligo_Length']==PPLengthDict['Primer_Len'])]


		if ProbeRange == True:
			ProbePassedSubet = PassedOligos[(PassedOligos['Oligo_Type']=='Probe') &
											 (PassedOligos['Oligo_Length']>=PPLengthDict['Probe_Min']) &
											 (PassedOligos['Oligo_Length']<=PPLengthDict['Probe_Max']) ]


		else:
			ProbePassedSubet = PassedOligos[(PassedOligos['Oligo_Type']=='Probe') &
											(PassedOligos['Oligo_Length']==PPLengthDict['Probe_Len'])]


		PassedOligos = pd.concat([PrimerPassedSubet,ProbePassedSubet]).reset_index().drop(['index'],axis=1)



		# If Background Check is Present filter
		if AllParameter.BackgroundCheck == 'YES':
			PassedOligos = PassedOligos[PassedOligos['Max Background Binding Score']<=AllParameter.CrossReactivityThresh]



		if (len(PrimerPassedSubet)==0 or len(ProbePassedSubet)==0 or len(PassedOligos) == 0):
			print('No Oligos Passed Filtering')
			sys.exit()


	#Create Binding DataFrame If It Doesnt Exist
	else:

		# Load In User Specified
		if AllParameter.PriorAlign != 'NO':
			print('Reading Alignment Summary')
			MergedAlignedDF = pd.read_csv('{}'.format(AllParameter.PriorAlign))
			HDFLST = list(range(len(MergedAlignedDF)))
			HomoDFInputIndexBlocks = [HDFLST[i:i + 1000] for i in range(0, len(HDFLST), 1000)]


		# Create Alignment DF if it doesnt exist
		else:
			print('Generating Alignment Summary')
			# Read in the aligned fasta sequences
			fastadict = {}
			with open(AllParameter.InputFile) as file_one:
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

			# Extract Alignment Length and QC
			FirstSeq = True
			FastaSeqLength = 0
			for fastaseq in fastadict.values():
				if FirstSeq == True:
					FastaSeqLength = len(fastaseq)
					FirstSeq = False
				if len(fastaseq) != FastaSeqLength:
					sys.exit('ERROR: Alignment file error length of all sequences is not equal')


			# Generate Neccessary Alignment Summary DF
			HDFLST = list(range(FastaSeqLength))
			HomoDFInputIndexBlocks = [HDFLST[i:i + 1000] for i in range(0, len(HDFLST), 1000)]
			MTDFOVI = list(zip([fastadict]*len(HomoDFInputIndexBlocks),HomoDFInputIndexBlocks))
			with Pool(processes=AllParameter.Threads) as pool:
				AlignedDFMultiThreadOupt = pool.starmap(CreatingInputHomologyDF,MTDFOVI)
			MergedAlignedDF = pd.concat(AlignedDFMultiThreadOupt).reset_index().drop(['index'],axis=1)
			MergedAlignedDF = MergedAlignedDF.sort_values(by=['Index_Pos'])
			MergedAlignedDF.to_csv('{}_Alignment_Summary.csv'.format(AllParameter.RunID),index=None)


		print('Generating Primer/Probe Binding Site DataFrame')
		# Generating Primer/Probe Binding Site DataFrame
		PrimerProbeCheckParallelInput = list(zip([AllParameter]*len(HomoDFInputIndexBlocks),
												 [MergedAlignedDF]*len(HomoDFInputIndexBlocks),
												 HomoDFInputIndexBlocks))
		with Pool(processes=AllParameter.Threads) as pool:
			PotentialPrimerProbeOut = pool.starmap(IndentifyingAndFilteringOligos,PrimerProbeCheckParallelInput)


		PassedOligos =  pd.concat(PotentialPrimerProbeOut).reset_index().drop(['index'],axis=1)
		if len(PassedOligos) == 0:
			print('No Oligos Passed Filtering')
			sys.exit()
		PassedOligos.to_csv('{}_PrimedRPA_Oligo_Binding_Sites.csv'.format(AllParameter.RunID),index=None)







	# Identify Primer-Probe Sets
	print('Identifying Valid Primer-Probe Combinations')
	FinalOutputDF = pd.DataFrame()

	# Determine all probe lengths available
	if AllParameter.ProbeRequired in ['EXO','NFO']:
		PossibleProbeLengths = PassedOligos[PassedOligos['Oligo_Type']=='Probe'].loc[:,'Oligo_Length'].unique().tolist()

	else:
		PossibleProbeLengths = [0]

	# Loop through possible FP  lengths
	for SFPL in PassedOligos[PassedOligos['Oligo_Type']=='Primer'].loc[:,'Oligo_Length'].unique():

		# Identify all possible FP starting sites
		PossibleForwardPrimerSS = PassedOligos[(PassedOligos['Oligo_Type']=='Primer')&(PassedOligos['Oligo_Length']==SFPL)].loc[:,'Binding_Site_Start_Index'].tolist()

		# Loop through all possible RP lengths
		for SRPL in PassedOligos[PassedOligos['Oligo_Type']=='Primer'].loc[:,'Oligo_Length'].unique():

			# Identify all possible RP starting sites
			PossibleReversePrimerSS = PassedOligos[(PassedOligos['Oligo_Type']=='Primer')&(PassedOligos['Oligo_Length']==SRPL)].loc[:,'Binding_Site_Start_Index'].tolist()


			# Split FP sites into tranches based on number of threads available.
			random.shuffle(PossibleForwardPrimerSS)
			PrimerStartSiteTranchesOne = [PossibleForwardPrimerSS[i:i + 2] for i in range(0, len(PossibleForwardPrimerSS), 2)]
			PrimerStartSiteTranchesTwo = [PrimerStartSiteTranchesOne[i:i + AllParameter.Threads] for i in range(0, len(PrimerStartSiteTranchesOne), AllParameter.Threads)]

			# Loop through possible probe lengths.
			for PPL in PossibleProbeLengths:

				if PPL == 0:
					PossibleProbeSS = []

				else:
					PossibleProbeSS = PassedOligos[(PassedOligos['Oligo_Type']=='Probe')&(PassedOligos['Oligo_Length']==PPL)].loc[:,'Binding_Site_Start_Index'].tolist()

				for PSST in PrimerStartSiteTranchesTwo:
					if len(FinalOutputDF) < AllParameter.MaxSets: ## Check this threshold setting.
						ComboSearchInput = list(zip(PSST,
												  [PossibleReversePrimerSS]*len(PSST),
												  [PossibleProbeSS]*len(PSST),
												  [AllParameter]*len(PSST),
												  [PassedOligos]*len(PSST),
												  [SFPL]*len(PSST),
												  [SRPL]*len(PSST),
												  [PPL]*len(PSST)))



						with Pool(processes=AllParameter.Threads) as pool:
							SuccessfulSets = pool.starmap(ComboIdentifyier,ComboSearchInput)

						TempOutputDF = pd.concat(SuccessfulSets)
						FinalOutputDF = pd.concat([FinalOutputDF,TempOutputDF])
						if len(FinalOutputDF) != 0:
							FinalOutputDF.to_csv('{}_Output_Sets.csv'.format(AllParameter.RunID),index=None)


	print('PrimedRPA Complete')
	print('If helpful, please cite:\n\nHiggins M et al. Submitted. 2018\nDOI: 10.1093/bioinformatics/bty701')


########################################################################################################################
########################################################################################################################
########################################################################################################################


# Basic Command Line User Interface
print('\n\n')
print('-------------------------------------------')
print('----------------PrimedRPA------------------')
print('-----Finding RPA Primer and Probe Sets-----')
print('-------------Higgins M et al.--------------')
print('-------------------------------------------\n\n')

# Utilise Either Parameters File Or Command Line Interface
if 'PrimedRPA_Parameters.txt' in sys.argv[1]:

	# Import and Extract Parameters
	parametersFile = sys.argv[1]
	try:
		paraFile = open(parametersFile,"r")
		iu = []
		for line in paraFile.readlines():
			if ">" in line:
				n = line.strip('\n')
				h = n.strip('>')
				iu.append(h)
		u = iu[1:]

		class AllParameter:
			RunID = str(u[0])
			PriorAlign = str(u[1])
			PriorBindingSite = str(u[2])
			InputFile = str(u[3])
			InputFileType = str(u[4])
			IdentityThreshold = int(u[5])/100
			PrimerLength = u[6].strip()
			ProbeRequired = str(u[7]).upper()
			ProbeLength = u[8].strip()
			AmpliconSizeLimit = int(u[9])
			NucleotideRepeatLimit = int(u[10])
			MinGC = int(u[11])
			MaxGC = int(u[12])
			DimerisationThresh = int(u[13])
			BackgroundCheck = str(u[14])
			CrossReactivityThresh = int(u[15])
			MaxSets = int(u[16])
			Threads = int(u[17])


	except:
		print('Parameters File Could Not Be Opened\nCheck File Path.')
		sys.exit()

else:

	parser = argparse.ArgumentParser()
	parser.add_argument('--RunID', help='Desired Run ID', required=True)
	parser.add_argument('--PriorAlign', help='Optional: Path to Prior Binding File',default='NO')
	parser.add_argument('--PriorBindingSite', help='Optional: Path to Prior Binding File',default='NO')
	parser.add_argument('--InputFile', help='Path to Input File',default='NO')
	parser.add_argument('--InputFileType', help='Options [SS,MS,MAS]')
	parser.add_argument('--IdentityThreshold', help='Desired Identity Threshold',default=0.99)
	parser.add_argument('--PrimerLength', help='Desired Primer Length',type=str,default='30')
	parser.add_argument('--ProbeRequired', help='Options [NO,EXO,NFO]',type=str,default='NO')
	parser.add_argument('--ProbeLength', help='Desired Probe Length',default='50')
	parser.add_argument('--AmpliconSizeLimit', help='Amplicon Size Limit',default=500)
	parser.add_argument('--NucleotideRepeatLimit', help='Nucleotide Repeat Limit (i.e 5 = AAAAA)',default=5)
	parser.add_argument('--MinGC', help='Minimum GC Content',default=30)
	parser.add_argument('--MaxGC', help='Maximum GC Content',default=70)
	parser.add_argument('--DimerisationThresh', help='Percentage Dimerisation Tolerated',default=40)
	parser.add_argument('--BackgroundCheck', help='Options [NO, Path To Background Fasta File]',default='NO')
	parser.add_argument('--CrossReactivityThresh', help='Max Cross Reactivity Threshold',default=60)
	parser.add_argument('--MaxSets', help='Default Set To 100',default=100)
	parser.add_argument('--Threads', help='Default Set To 1',default=1)
	AllParameter = parser.parse_args()



# Check Files Defined Exist
FileLocationsCheckList = [AllParameter.PriorAlign,
						  AllParameter.PriorBindingSite,
						  AllParameter.BackgroundCheck,
						  AllParameter.InputFile]

for FL in FileLocationsCheckList:
	if FL != 'NO':
		if FL not in glob.glob(FL):
			print('File Not Found: {}'.format(FL))
			sys.exit()



#Define tool dependancy paths
InstalledSourceScriptPath = os.path.realpath(__file__)
AllParameter.BlastnPath = InstalledSourceScriptPath.replace('PrimedRPA.py','Tool_Dependancies/blastn')
AllParameter.BlastnDBCreationPath = InstalledSourceScriptPath.replace('PrimedRPA.py','Tool_Dependancies/makeblastdb')
AllParameter.ClustalOPath = InstalledSourceScriptPath.replace('PrimedRPA.py','Tool_Dependancies/clustalo')


# Create Blastn Database if neccessary
if AllParameter.BackgroundCheck.upper() != 'NO':

	AllParameter.BlastnDBName =  '{}_Blastn_DB_PrimedRPA'.format(AllParameter.BackgroundCheck.split('/')[-1].split('.')[0])
	if os.path.exists(AllParameter.BlastnDBName) == False:
		os.makedirs(AllParameter.BlastnDBName)
		makeblastdbcommand = '{0} -in {1} -dbtype nucl -parse_seqids -out {2}/{2}'.format(AllParameter.BlastnDBCreationPath,
																						  AllParameter.BackgroundCheck,
																						  AllParameter.BlastnDBName)

		print('\nCreating Blastn Database:')
		subprocess.call([makeblastdbcommand],shell=True)
		print('\n')

	# Create directory to house blastn outputs
	AllParameter.PrimerBlastnOutput = '{}/{}'.format(AllParameter.BlastnDBName,AllParameter.RunID)
	if os.path.exists(AllParameter.PrimerBlastnOutput) == False:
		os.mkdir(AllParameter.PrimerBlastnOutput)


# Check If Alignment Is Necessary
if AllParameter.InputFileType == 'MS':

	# Check if MS Alignment already Exists in Working Directory
	if AllParameter.InputFile.replace('.fasta','_Aligned.fasta') not in glob.glob(AllParameter.InputFile.replace('.fasta','_Aligned.fasta')):
		RunningClustalo1(AllParameter.ClustalOPath,[AllParameter.InputFile], overwriteOutput=False)

	AllParameter.InputFile = AllParameter.InputFile.replace('.fasta','_Aligned.fasta')


# Run Primer, Probe Design process
CheckingAlignedOutputFile(AllParameter)
