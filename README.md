#### What is LOCALIZER?

LOCALIZER is a machine learning method for predicting the subcellular localization of both plant proteins and pathogen effectors in the plant cell. It can currently predict localization to chloroplasts and mitochondria using transit peptide prediction and to nuclei using a collection of nuclear localization signals (NLSs). 
 
For plant protein localization prediction, submit full-length sequences and run it in 'plant mode'.

For effector protein localization prediction, submit full-length sequences and run it in 'effector mode'. 
It is recommended to use tools such as SignalP or Phobius	to predict first if a protein is likely to be secreted and to obtain the mature sequences. Alternatively, provide full sequences and let LOCALIZER delete the first 20 aas as the putative signal peptide region.

Do not submit short sequence fragments to LOCALIZER, it expects the full protein sequences. 

You can submit your proteins of interest to the webserver at http://localizer.csiro.au/ or install it locally.
All training and evaluation data can be found [here](http://localizer.csiro.au/data.html).

#### Installing LOCALIZER
LOCALIZER has been written in Python and uses pepstats from the EMBOSS software and the WEKA 3.6 software. It also requires that you have Perl and BioPython installed. **LOCALIZER from version 1.0.4 inclusive uses Python 3.** 

To get LOCALIZER to work on your local machine, you need to install the EMBOSS and WEKA softwares from source. Both are already provided in the LOCALIZER distribution to ensure that compatible versions are used. 

1. Make sure LOCALIZER has the permission to execute. Then unpack LOCALIZER in your desired location
```
tar xvf LOCALIZER_1.0.4.tar.gz
chmod -R 755 LOCALIZER_1.0.4/
cd LOCALIZER_1.0.4 
```

2. For the EMBOSS installation, you need to switch to the Scripts directory and unpack, configure and make. Alternatively, if you are on a computer cluster and EMBOSS is already installed, you can change the variable PEPSTATS_PATH in the LOCALIZER.py script to the EMBOSS location on the machine you are using.
```
cd Scripts
tar xvf emboss-latest.tar.gz
cd EMBOSS-6.5.7/
./configure
make
cd ../ 
```

3. For WEKA, you need to simply unzip the file weka-3-6-12.zip
```
unzip weka-3-6-12.zip
```
If you are having troube installing EMBOSS, please see [here](http://emboss.sourceforge.net/download/) for help.
If you are having troube installing WEKA, please see [here](https://www.cs.waikato.ac.nz/~ml/weka/index.html) for help. 

4. Test if LOCALIZER is working
```
python LOCALIZER.py -e -i Effector_Testing.fasta
```

5. Problems?

If you are getting an error message like 'ImportError: No module named Bio', you need to install BioPython on your computer. See [here](https://biopython.org/wiki/Download) for help. For example, you can try and run:
```
pip install biopython
```

Note also that you need PERL to be installed on your computer for running NLStradamus. 

#### Running LOCALIZER
Run this to get a feel for the output format:
```
python LOCALIZER.py -e -i Effector_Testing.fasta

# -----------------
# LOCALIZER 1.0.4 Predictions (-e mode)
# -----------------
Identifier      Chloroplast             Mitochondria            Nucleus
CRN15           -                       -                       Y (KRKR)
Ecp2            -                       -                       -
AVR-Pii         -                       -                       -
ToxA            Y (0.877 | 62-130)      -                       -
--------------------------------------
--------------------------------------
# Proteins analyzed: 4 from file: Effector_Testing.fasta

# Number of proteins with cTP: 1 (25.0%)
# Number of proteins with cTP & possible mTP: 0 (0.0%)
# Number of proteins with cTP & NLS: 0 (0.0%)
# Number of proteins with cTP & possible mTP & NLS: 0 (0.0%)
# Number of proteins with mTP: 0 (0.0%)
# Number of proteins with mTP & possible cTP: 0 (0.0%)
# Number of proteins with mTP & NLS: 0 (0.0%)
# Number of proteins with mTP & possible cTP & NLS: 0 (0.0%)
# Number of proteins with NLS and no transit peptides: 1 (25.0%)
--------------------------------------
--------------------------------------
# Summary statistics

# Number of proteins with chloroplast localization (cTP, cTP & possible mTP, cTP & NLS, cTP & possible mTP & NLS): 1 (25.0%)
# Number of proteins with mitochondrial localization (mTP, mTP & possible cTP, mTP & NLS, mTP & possible cTP & NLS): 0 (0.0%)
# Number of proteins with nuclear localization and no transit peptides: 1 (25.0%)
# Number of proteins with nuclear localization and with transit peptides: 0 (0.0%)
--------------------------------------
--------------------------------------

```
LOCALIZER will return the output as shown in the example above. First, a summary table will be shown which shows the predictions (chloroplast, mitochondria or nucleus) for each submitted protein. If a transit peptide is predicted, the start and end positions in the submitted sequences are shown, alongside the probability. In this example, ToxA has a predicted chloroplast transit peptide with probability 0.885 at position 62-130 in its sequence. LOCALIZER does not return a probability for nucleus localization, because it uses a simple NLS search. In this example, LOCALIZER found a NLS in CRN15, i.e. the sequence KRKR.

In the summary statistic, we count LOCALIZER predictions that are 'chloroplast', 'chloroplast and possible mitochondrial', 'chloroplast and nucleus' and 'chloroplast & possible mitochondrial and nucleus' as chloroplast predictions (same strategy for mitochondrial predictions). A protein that carries a predicted transit peptide with an additional predicted NLS might have experimental evidence only for one of those locations due to the technical hurdles of recognizing dual targeting and should thus not necessarily be counted as a false positive prediction. However, in the LOCALIZER paper, a protein was counted as a nucleus prediction only if it has the category 'nucleus' to avoid assigning a protein to multiple predictions in the evaluation and this is what we recommend. 


#### Citation for LOCALIZER:

Sperschneider, J., Catanzariti, A., DeBoer, K. et al. LOCALIZER: subcellular localization prediction of both plant and effector proteins in the plant cell. Sci Rep 7, 44598 (2017) doi:10.1038/srep44598
