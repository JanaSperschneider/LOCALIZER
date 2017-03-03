-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
LOCALIZER: subcellular localization prediction of plant and effector proteins in the plant cell
Copyright (C) 2016-2017 Jana Sperschneider	
Contact: jana.sperschneider@csiro.au
-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
Installation instructions for LOCALIZER 1.0
-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
LOCALIZER relies on two tools, the EMBOSS software and the WEKA 3.6 software. These have been shipped 
with the EffectorP 1.0.tar.gz file, but they need to be installed by the user. 

1) Extract the LOCALIZER 1.0 archive:

-----------------------------------------
tar xvf LOCALIZER_1.0.tar.gz
cd LOCALIZER_1.0 
-----------------------------------------

2) Install EMBOSS

-----------------------------------------
cd Scripts
tar xvf emboss-latest.tar.gz 
cd EMBOSS-6.5.7/
./configure
make
cd ../
-----------------------------------------

3) Install WEKA: simply unzip the file weka-3-6-12.zip

-----------------------------------------
unzip weka-3-6-12.zip
-----------------------------------------

3) Run LOCALIZER

To test that LOCALIZER is working, type the following command in the working directory LOCALIZER_1.0/Scripts

-----------------------------------------
python LOCALIZER.py -e -i Effector_Testing.fasta
-----------------------------------------

Note that LOCALIZER runs under Python 2.x, not under Python 3.x.
Note also that NLStradamus for predicting NLSs needs PERL to be installed on your computer.
You also need BioPython on your computer for running ProtParam (http://biopython.org/DIST/docs/install/Installation.html)

If you are having troube installing EMBOSS, please see here for help: http://emboss.sourceforge.net/download/
If you are having troube installing WEKA, please see here for help: http://www.cs.waikato.ac.nz/~ml/weka/index.html

