#### What is LOCALIZER?

LOCALIZER is a method for predicting the subcellular localization of both plant proteins and pathogen effectors in the plant cell. It can currently predict localization to chloroplasts and mitochondria using transit peptide prediction and to nuclei using a collection of nuclear localization signals (NLSs). 
 
For plant protein localization prediction, submit full-length sequences and run it in 'plant mode'.

For effector protein localization prediction, submit full-length sequences and run it in 'effector mode'. 
It is recommended to use tools such as SignalP or Phobius	to predict first if a protein is likely to be secreted and to obtain the mature sequences. Alternatively, provide full sequences and let LOCALIZER delete the first 20 aas as the signal peptide region.

Do not submit short sequence fragments to LOCALIZER, it expects the full sequence. 
 
#### Citation for LOCALIZER:

Sperschneider, J., Catanzariti, A., DeBoer, K. et al. LOCALIZER: subcellular localization prediction of both plant and effector proteins in the plant cell. Sci Rep 7, 44598 (2017) doi:10.1038/srep44598

#### Running LOCALIZER

You can submit your proteins of interest to the webserver at http://localizer.csiro.au/.

Alternatively, you can install LOCALIZER on your machine to run it locally. 
For detailed installation instructions see here: http://localizer.csiro.au/software.html

For help on how to interpret the output format, see http://localizer.csiro.au/output.html
