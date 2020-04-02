#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

PackList = ['os','shutil','sys','pandas',
            'math','glob','time','itertools',
            'datetime','random','numpy','multiprocess']


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
PrimedRPAPath = os.path.dirname(os.path.realpath(__file__))


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
    AliasFile = 'bash_profile'

elif OperatingSystem == "win32" or 'win64':
    print('Windows Operating System Detected\n\
    This is not supported.\n')
    sys.exit()

for dp in ['clustalo','blastn','makeblastdb']:
    subprocess.run('chmod u+x ./Tool_Dependancies/{}'.format(dp),shell=True)


# Install And Establish Samtools - Commands Shown Below
SamtoolsInstall = 'cd {0}/Tool_Dependancies/ ;wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 ;tar -xf samtools-1.10.tar.bz2; cd ./samtools-1.10 ; ./configure --prefix=/{0}/Tool_Dependancies/; make; make install; chmod u+x {0}/Tool_Dependancies/bin/samtools; cd {0}'.format(PrimedRPAPath)
subprocess.run(SamtoolsInstall,shell=True)



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
ValidationCommand = './PrimedRPA.py ./PrimedRPA_Parameters.txt'

try:
    subprocess.run('chmod +x PrimedRPA.py',shell=True)
    subprocess.run(ValidationCommand,shell=True,check=True)
    print('\n\nValidation Run Successful')

except subprocess.CalledProcessError as e:
    print(e)
    print('PrimedRPA Validation Run Fail\nInstallation Exit')
    sys.exit()



FilesToRemove = ['Installation_Validation_Run_Output_Sets.csv',
                 'Installation_Validation_Run_Alignment_Summary.csv',
                 'Installation_Validation_Run_PrimedRPA_Oligo_Binding_Sites.csv']

for i in FilesToRemove:
    if os.path.exists(i):
        os.remove(i)

shutil.rmtree('./Validation_1_Background_Blastn_DB_PrimedRPA')

print('Installation Successfully Complete')
