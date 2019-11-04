#! /usr/bin/env python3

import subprocess

PackList = ['os','shutil','sys','pandas',
            'math','glob','time','itertools',
            'datetime','random','numpy']


# Checking Package Installation and if not present install
for DepPack in PackList:
    try:
        __import__(DepPack)
    except ImportError as e:
        print(e)
        print('Installing Python Package: {}'.format(DepPack))
        subprocess.run('pip install {}'.format(DepPack),shell=True)
        pass

# Install Necessary Packages
import os
import shutil
import sys
import pandas as pd

OperatingSystem = sys.platform


# Stage One = Sorting Out Program Dependencies Folder
if OperatingSystem in ["linux","linux2"]:
    print('Linux Operating System Detected')

    if os.path.isdir("./Tool_Dependancies_Mac") == True:
        shutil.rmtree('./Tool_Dependancies_Mac')
        #os.rename('Tool_Dependancies_Linux','Tool_Dependancies')
        subprocess.run('mv Tool_Dependancies_Linux Tool_Dependancies',shell=True)
    AliasFile = 'bashrc'

elif OperatingSystem == "darwin":
    print('Darwin Operating System Detected')

    if os.path.isdir("./Tool_Dependancies_Linux") == True:
        shutil.rmtree('./Tool_Dependancies_Linux')
        subprocess.run('mv Tool_Dependancies_Mac Tool_Dependancies',shell=True)
    #os.rename('Tool_Dependancies_Mac','Tool_Dependancies')
    AliasFile = 'bash_profile'

elif OperatingSystem == "win32" or 'win64':
    print('Windows Operating System Detected\n\
    This is not supported.\n')
    sys.exit()

for dp in ['clustalo','blastn','makeblastdb']:
    subprocess.run('chmod u+x ./Tool_Dependancies/{}'.format(dp),shell=True)


# Stage Two = File Presence
NecessaryFilePaths = ['PrimedRPA.py',
                      'PrimedRPA_Parameters.txt',
                      './Tool_Dependancies/blastn',
                      './Tool_Dependancies/makeblastdb',
                      './Tool_Dependancies/clustalo']

for FilePath in NecessaryFilePaths:
    if os.path.isfile(FilePath) == False:
        print('{}: Not Found\nInstallation Exit'.format(FilePath))
        sys.exit()


# Stage Three = Edit bashrc
PrimedRPAPath = os.path.dirname(os.path.realpath(__file__))

subprocess.run('cp ~/.{0} ~/.{0}_Prior_PrimedRPA'.format(AliasFile),shell=True)
subprocess.run('echo alias PrimedRPA="{0}/PrimedRPA.py" >> ~/.{1}'.format(PrimedRPAPath,AliasFile),shell=True)

# Stage Four Triggure basic Validation.
subprocess.run('source ~/.{}'.format(AliasFile),shell=True)

ValidationCommand = './PrimedRPA.py ./Validation/Validation_1_PrimedRPA_Parameters.txt'

try:
    subprocess.run('chmod +x PrimedRPA.py',shell=True)
    subprocess.run(ValidationCommand,shell=True,check=True)

except subprocess.CalledProcessError as e:
    print(e)
    print('PrimedRPA Validation Run 1 Fail\nInstallation Exit')
    sys.exit()


ValOneDF = pd.read_csv('Validation_1_PrimedRPA_Output.csv')
ValOneDF = ValOneDF.sort_values(by=['Max_Blastn_Percentage_Match'],ascending=True)


ExpectedResult = pd.DataFrame([['GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTC', 'TCCATGCATTTGGTATTTCGTCTGGCGGCTGTGCACGCGATAGCATTGCGAGACGCTGG', 'TACTTCAAAGACAGATACTGCGACATAGGGTGCTCCGGCTC', 17.073170731707318],
                                ['ATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCT', 'TCCATGCATTTGGTATTTCGTCTGGCGGCTGTGCACGCGATAGCATTGCGAGACGCTGG', 'TACTTCAAAGACAGATACTGCGACATAGGGTGCTCCGGCTC', 19.51219512195122],
                                ['GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTC', 'CTCCATGCATTTGGTATTTCGTCTGGCGGCTGTGCACGCGATAGCATTGCGAGACGCTG', 'ACTTCAAAGACAGATACTGCGACATAGGGTGCTCCGGCTCC', 17.073170731707318],
                                ['GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTC', 'CTCCATGCATTTGGTATTTCGTCTGGCGGCTGTGCACGCGATAGCATTGCGAGACGCTG', 'TACTTCAAAGACAGATACTGCGACATAGGGTGCTCCGGCTC', 17.073170731707318]],
                                columns=["Forward_Primer", "Nfo_Probe", "Reverse_Primer", "Max_Blastn_Percentage_Match"])





counter = 0
for i in ValOneDF.values.tolist():
    if i in ExpectedResult.values.tolist():
        counter+=1

if counter == len(ExpectedResult):
    print('Validation 1 Passed\n')

else:
    print('Validation 1 Failed\nInstallation Exit')
    sys.exit()

FilesToRemove = ['Validation_1_PrimedRPA_Output.csv']
for i in FilesToRemove:
    if os.path.exists(i):
        os.remove(i)

shutil.rmtree('./Validation_1_Background_Blastn_DB_PrimedRPA')

print('Installation Successfully Complete')
