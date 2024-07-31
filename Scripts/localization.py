#!/usr/bin/env python3
"""
    LOCALIZER: subcellular localization prediction of plant  
    and effector proteins in the plant cell
    Copyright (C) 2016-2017 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 3 of the License, or     
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Please also see the CSIRO Disclaimer provided with LOCALIZER (LICENCE.txt).

    Contact: jana.sperschneider@csiro.au
"""
import os
import sys
import functions
import subprocess
import errno
import uuid
import shutil
import re

import parameters
from operator import itemgetter

#from Bio import SeqIO
#from Bio.SeqUtils import ProtParam
# -----------------------------------------------------------------------------------------------------------
def GRAVY(sequence):
    # The GRAVY value for a peptide or protein is calculated as the sum of 
    # hydropathy values of all the amino acids, 
    # divided by the number of residues in the sequence. 

    gravy = 0.0
    for aa in sequence:
        if aa.upper() in parameters.GRAVY_DIC:
            gravy += parameters.GRAVY_DIC[aa.upper()]

    gravy = gravy/len(sequence)

    return gravy
# -----------------------------------------------------------------------------------------------------------
def bipart_NLS(sequence):
    '''
  * Bipartite nuclear localization signal profile *
  http://prosite.expasy.org/cgi-bin/prosite/get-prodoc-entry?PDOC00015
  (1) Two adjacent basic amino acids (Arg or Lys).
  (2) A spacer region of any 8-12 residues.
  (3) At least three basic residues (Arg or Lys) in the five positions
      after the spacer region.
    '''
    motif = []
    spacer_region = [8,9,10,11,12]

    for spacer in spacer_region:
        for position in range(0, len(sequence) - (7 + spacer) + 1, 1):
            # Two adjacent basic amino acids (Arg or Lys)
    	    if sequence[position] == 'R' or sequence[position] == 'K':
                if sequence[position + 1] == 'R' or sequence[position + 1] == 'K':
                    sub_seq = sequence[position+spacer+2:position+spacer+7]
                    if sub_seq.count('R') + sub_seq.count('K') >= 3.0:
                        motif.append(sequence[position:position+spacer+7])           

    return motif
# -----------------------------------------------------------------------------------------------------------
def predict_NLS(seq):
    """ Function: predict_NLS()

        Purpose:  Given a protein sequence, search for a list of NLS motifs
                  using regular expressions and exact matches.
              
        Input:    Protein sequence.
    
        Return:   List of NLS motifs that were found in the sequence.
    """ 
    NLS_motifs = []

    for motif in parameters.NLS_MOTIFS:
        match = re.search(motif, seq)
        if match:	
            NLS_motifs.append(match.group())

    return NLS_motifs
# -----------------------------------------------------------------------------------------------------------
def NLStradamus(seq, TMP_PATH, SCRIPT_PATH):
    """ Function: NLStradamus()

        Purpose:  Given a protein sequence, run NLStradamus to predict NLS motifs.
              
        Input:    Protein sequence.
    
        Return:   List of NLStradamus NLS motifs that were found in the sequence
                  with score, start/end position and motif.
    """ 
    input_file = TMP_PATH + 'NLS_input.fasta'
    output_file = TMP_PATH + 'Nucleus.NLStradamus'

    write_FASTA(input_file, ['>input\n'], [seq])

    ParamList = ['perl', SCRIPT_PATH + '/nlstradamus.pl', '-tab', '-i', 
             input_file]

    with open(output_file, 'wb') as out:
        try:
            Process = subprocess.Popen(ParamList, shell=False, stdout=out)
            sts = Process.wait()
            cstdout, cstderr = Process.communicate()

            if Process.returncode:
                raise Exception("Calling NLStradamus returned %s"%Process.returncode)
            if cstdout:
                pass
            elif cstderr:
                sys.exit()
        except:
            e = sys.exc_info()[1]
            print("Error calling NLStradamus: %s" % e)
            sys.exit(1)

    NLStradamus = []

    f = open(output_file, 'r')
    content = f.readlines()    
    f.close()

    for line in content:
        if not line.startswith('#'):
            identifier = line.split('\t')[0]
            score = float(line.split('\t')[2])
            start, stop = int(line.split('\t')[3]), int(line.split('\t')[4])
            sequence = line.split('\t')[5]
            sequence = sequence.strip()

            # Motif needs to be >= 4 aas
            if stop - start + 1 >= 4.0:
                if identifier in NLStradamus:
                    NLStradamus.append((score, start, stop, sequence))
                else:
                    NLStradamus = [(score, start, stop, sequence)]

    return NLStradamus
# -----------------------------------------------------------------------------------------------------------
def map_coordinates(transit_peptide, seq, full_seq):
    """ Function: map_coordinates()

        Purpose:  Given a transit peptide with its coordinates in the mature sequence, 
                  calculate the coordinates in the full sequence.
              
        Input:    Transit peptide coordinates, mature sequence, full sequence
    
        Return:   Transit peptide coordinates in full sequence
    """ 
    coordinates = []
    # Now map the coordinates from the mature sequence to the full sequence
    for item in transit_peptide:
        start, end, prob, marker = item[0], item[1], item[2], item[3]
        SP_end = full_seq.index(seq)
        coordinates = [int(start + SP_end), int(end + SP_end), prob, marker]

    return coordinates  
# -----------------------------------------------------------------------------------------------------------
def write_FASTA(f_output, IDENTIFIERS, SEQUENCES):

    f = open(f_output, 'w')

    for ident, seq in zip(IDENTIFIERS, SEQUENCES):
        f.writelines(ident)
        f.writelines(seq + '\n')
    f.close()

    return
# -----------------------------------------------------------------------------------------------------------
def chloro_classifier_allwindows(pepstats_dic_aas, pepstats_dic_aas_short, IDENTIFIERS, SEQUENCES, TMP_PATH, WEKA_PATH, SCRIPT_PATH):

    predicted_chloro = []
    weka = TMP_PATH + 'chloroplast.arff'

    with open(weka, 'w') as f:

        # Create a list of features for each protein
        X = [[] for __ in range(len(IDENTIFIERS))]

        for protein_position, TARGET_ID in enumerate(IDENTIFIERS):
            TARGET_ID = TARGET_ID.replace('>', '')
            TARGET_ID = TARGET_ID.strip()
            sequence = SEQUENCES[protein_position]

            molecular_weight, charge, isoelectric, amino_acid_classes, amino_acid_frequencies = pepstats_dic_aas[TARGET_ID]

            #prot = ProtParam.ProteinAnalysis(sequence.replace('*',''))            
            
            molecular_weight_short, charge_short, isoelectric_short, amino_acid_classes_short, amino_acid_frequencies_short = pepstats_dic_aas_short[TARGET_ID]

            helix = (sequence.count('V') + sequence.count('I') + sequence.count('Y') + sequence.count('F') + sequence.count('W') + sequence.count('L'))/len(sequence.replace('*',''))
            turn = (sequence.count('N') + sequence.count('P') + sequence.count('G') + sequence.count('S'))/len(sequence.replace('*',''))
            sheet = (sequence.count('E') + sequence.count('M') + sequence.count('A') + sequence.count('L'))/len(sequence.replace('*',''))

            aromaticity = sequence.count('F')/len(sequence.replace('*',''))
            aromaticity += sequence.count('W')/len(sequence.replace('*',''))
            aromaticity += sequence.count('Y')/len(sequence.replace('*',''))

            X[protein_position] = [charge, isoelectric] + amino_acid_classes + amino_acid_frequencies + [GRAVY(sequence)] + [helix, turn, sheet, aromaticity] + [charge_short, isoelectric_short] + amino_acid_frequencies_short

            #X[protein_position] = [charge, isoelectric] + amino_acid_classes + amino_acid_frequencies + [GRAVY(sequence)] + [prot.secondary_structure_fraction()[0], prot.secondary_structure_fraction()[1], prot.secondary_structure_fraction()[2], prot.aromaticity()] + [charge_short, isoelectric_short] + amino_acid_frequencies_short   

        f.writelines(parameters.ARFF_CHLOROPLAST_HEADER)
        for index, vector in enumerate(X):
            for feature in vector:
                f.writelines(str(feature) + ',')
            f.writelines('?\n')

    ParamList = ['java', '-cp', WEKA_PATH, 'weka.classifiers.functions.SMO',
             '-l', SCRIPT_PATH + '/MODEL_FILES/CHLOROPLAST_NOTCHLOROPLAST.model',
             '-T', weka, '-p', 'first-last']

    with open(TMP_PATH + 'Chloroplast_Predictions.txt', 'wb') as out:
        try:
            Process = subprocess.Popen(ParamList, shell=False, stdout=out)
            sts = Process.wait()
            cstdout, cstderr = Process.communicate()

            if Process.returncode:
                raise Exception("Calling WEKA returned %s"%Process.returncode)
            if cstdout:
                pass
            elif cstderr:
                sys.exit()
        except:
            e = sys.exc_info()[1]
            print("Error calling WEKA: %s" % e)
            sys.exit(1)

    file_input = TMP_PATH + 'Chloroplast_Predictions.txt'
    file_output = TMP_PATH + 'Chloroplast_Predictions.fasta'
    predicted_chloro = parse_weka_output(file_input, IDENTIFIERS, SEQUENCES, 'Chloroplast', 'Non-Chloroplast')

    return predicted_chloro
# -----------------------------------------------------------------------------------------------------------
def mito_classifier_allwindows(pepstats_dic_aas, pepstats_dic_aas_short, IDENTIFIERS, SEQUENCES, TMP_PATH, WEKA_PATH, SCRIPT_PATH):

    predicted_mito = []
    weka = TMP_PATH + 'mito.arff'

    with open(weka, 'w') as f:

        # Create a list of features for each protein
        X = [[] for __ in range(len(IDENTIFIERS))]

        for protein_position, TARGET_ID in enumerate(IDENTIFIERS):
            TARGET_ID = TARGET_ID.replace('>', '')
            TARGET_ID = TARGET_ID.strip()
            sequence = SEQUENCES[protein_position]

            molecular_weight, charge, isoelectric, amino_acid_classes, amino_acid_frequencies = pepstats_dic_aas[TARGET_ID]

            #prot = ProtParam.ProteinAnalysis(sequence.replace('*',''))

            molecular_weight_short, charge_short, isoelectric_short, amino_acid_classes_short, amino_acid_frequencies_short = pepstats_dic_aas_short[TARGET_ID]

            helix = (sequence.count('V') + sequence.count('I') + sequence.count('Y') + sequence.count('F') + sequence.count('W') + sequence.count('L'))/len(sequence.replace('*',''))
            turn = (sequence.count('N') + sequence.count('P') + sequence.count('G') + sequence.count('S'))/len(sequence.replace('*',''))
            sheet = (sequence.count('E') + sequence.count('M') + sequence.count('A') + sequence.count('L'))/len(sequence.replace('*',''))

            aromaticity = sequence.count('F')/len(sequence.replace('*',''))
            aromaticity += sequence.count('W')/len(sequence.replace('*',''))
            aromaticity += sequence.count('Y')/len(sequence.replace('*',''))

            X[protein_position] = [charge, isoelectric] + amino_acid_classes + amino_acid_frequencies + [GRAVY(sequence)] + [helix, turn, sheet, aromaticity] + [charge_short, isoelectric_short] + amino_acid_frequencies_short

            #X[protein_position] = [charge, isoelectric] + amino_acid_classes + amino_acid_frequencies + [GRAVY(sequence)] + [prot.secondary_structure_fraction()[0], prot.secondary_structure_fraction()[1], prot.secondary_structure_fraction()[2], prot.aromaticity()] + [charge_short, isoelectric_short] + amino_acid_frequencies_short  
 
        f.writelines(parameters.ARFF_MITOCHONDRIA_HEADER)
        for index, vector in enumerate(X):
            for feature in vector:
                f.writelines(str(feature) + ',')
            f.writelines('?\n')

    ParamList = ['java', '-cp', WEKA_PATH, 'weka.classifiers.functions.SMO',
             '-l', SCRIPT_PATH + '/MODEL_FILES/MITOCHONDRIA_NOTMITOCHONDRIA.model',
             '-T', weka, '-p', 'first-last']

    with open(TMP_PATH + 'Mitochondria_Predictions.txt', 'wb') as out:
        try:
            Process = subprocess.Popen(ParamList, shell=False, stdout=out)
            sts = Process.wait()
            cstdout, cstderr = Process.communicate()

            if Process.returncode:
                raise Exception("Calling WEKA returned %s"%Process.returncode)
            if cstdout:
                pass
            elif cstderr:
                sys.exit()
        except:
            e = sys.exc_info()[1]
            print("Error calling WEKA: %s" % e)
            sys.exit(1)

    file_input = TMP_PATH + 'Mitochondria_Predictions.txt'
    file_output = TMP_PATH + 'Mitochondria_Predictions.fasta'
    predicted_mito = parse_weka_output(file_input, IDENTIFIERS, SEQUENCES, 'Mitochondria', 'Non-Mitochondria')

    return predicted_mito
# -----------------------------------------------------------------------------------------------------------
def parse_weka_output(file_input, IDENTIFIERS, SEQUENCES, CLASS1, CLASS2):
    
    predicted_proteins = []

    with open(file_input) as f:

        content = f.readlines()            
        content = content[5:]

        for line in content:
            if line.strip():
                prediction = line.split()[2]
                position = line.split()[0]
                # WEKA output counts from position 1, our identifiers are counted from zero
                identifier = IDENTIFIERS[int(position)-1]
                sequence = SEQUENCES[int(position)-1]

                if '2:' in prediction:
                    prob = float(line.split()[3])                
                    negative = identifier.strip()
                    negative = negative.replace('>', '')  

                elif '1:' in prediction:
                    prob = float(line.split()[3])                
                    positive = identifier.strip()
                    positive = positive.replace('>', '')                
                
                    predicted_proteins.append((positive, prob, sequence))

    return predicted_proteins
# -----------------------------------------------------------------------------------------------------------
def pepstats_aas(WINDOW_IDS, WINDOW_SEQS, TMP_PATH, PEPSTATS_PATH):
    """ Function: pepstats_aas()

        Purpose:  Given a set that contains the list of identifiers and 
                  the corresponding list of sequences, scan the given 
                  pepstats result file to extract protein properties.
              
        Input:    Set that contains the list of identifiers and 
                  the corresponding list of sequences and peptstats 
                  result file. 
    
        Return:   Dictionary of protein properties for each protein in the 
                  set.
    """
    pepstats_dic = {}
    
    f = open(TMP_PATH + 'temp_pepstats.fasta', 'w')

    for ident, seq in zip(WINDOW_IDS, WINDOW_SEQS):
        f.writelines('>' + ident.replace('>','') + '\n')
        f.writelines(seq + '\n')
    f.close()

    # Call pepstats
    ProcessExe = PEPSTATS_PATH + 'pepstats'
    ParamList = [ProcessExe, '-sequence', TMP_PATH + 'temp_pepstats.fasta', 
              '-outfile', TMP_PATH + 'temp_pepstats.pepstats']

    try:
        FNULL = open(os.devnull, 'w')
        Process = subprocess.Popen(ParamList, shell=False, stdout=FNULL, stderr=subprocess.STDOUT)
        sts = Process.wait()
        cstdout, cstderr = Process.communicate()

        if Process.returncode:
            raise Exception("Calling pepstats returned %s"%Process.returncode)
        if cstdout:
            pass
        elif cstderr:
            sys.exit()
    except:
        e = sys.exc_info()[1]
        print("Error calling pepstats: %s" % e)
        sys.exit(1)

    with open(TMP_PATH + 'temp_pepstats.pepstats') as f:         
        content = f.readlines()
        for start, line in enumerate(content):
            if 'PEPSTATS of ' in line:
                TARGET_ID = line.split('PEPSTATS of ')[1]
                TARGET_ID = TARGET_ID.split('from 1 to')[0]
                TARGET_ID = TARGET_ID.strip()
                sequence = None
                for (identifier, seq) in zip(WINDOW_IDS, WINDOW_SEQS):
                    if identifier.strip() == TARGET_ID:
                        sequence = seq.strip()
                         
                if sequence:
                    length = float(len(sequence))
                    # Amino acid frequencies in the sequence
                    amino_acid_frequencies = []
                    amino_acid_frequencies.append(100.0*sequence.count('A')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('C')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('D')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('E')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('F')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('G')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('H')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('I')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('K')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('L')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('M')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('N')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('P')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('Q')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('R')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('S')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('T')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('V')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('W')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('Y')/length)
                    # Extract molecular weight
                    mwline = content[start + 2:start + 3]
                    molecular_weight = float(re.findall("\d+.\d+", str(mwline))[0])
                    # Extract charge
                    charge_line = content[start + 3:start + 4]
                    charge = float(re.findall("[-+]?\d+.\d+", str(charge_line))[1])

                    isoelectric_line = content[start + 4:start + 5]
                    # In rare cases the isoelectric point calculation returns 'None'
                    try:
                        isoelectric = float(re.findall("\d+.\d+", str(isoelectric_line))[0])
                    except: 
                        isoelectric = 0.0

		            # Watch out, in the pepstats software, if isoelectric point == None, an 
                    # extra line will be introduced
                    start_aas = content[start:].index('Property\tResidues\t\tNumber\t\tMole%\n')
                    perline = content[start + start_aas + 1:start + start_aas + 10]

                    tiny = float(re.findall("\d+.\d+", str(perline[0]))[-1])
                    small = float(re.findall("\d+.\d+", str(perline[1]))[-1])
                    aliphatic = float(re.findall("\d+.\d+", str(perline[2]))[-1])
                    aromatic = float(re.findall("\d+.\d+", str(perline[3]))[-1])
                    non_polar = float(re.findall("\d+.\d+", str(perline[4]))[-1])
                    polar = float(re.findall("\d+.\d+", str(perline[5]))[-1])
                    charged = float(re.findall("\d+.\d+", str(perline[6]))[-1])
                    basic = float(re.findall("\d+.\d+", str(perline[7]))[-1])
                    acidic = float(re.findall("\d+.\d+", str(perline[8]))[-1])
                    amino_acid_classes = []
                    amino_acid_classes.append(tiny)
                    amino_acid_classes.append(small)
                    amino_acid_classes.append(aliphatic)
                    amino_acid_classes.append(aromatic)
                    amino_acid_classes.append(non_polar)
                    amino_acid_classes.append(polar)
                    amino_acid_classes.append(charged)
                    amino_acid_classes.append(basic)
                    amino_acid_classes.append(acidic)
                    # Store values in dictionary
                    pepstats_dic[TARGET_ID] = molecular_weight, charge, isoelectric, amino_acid_classes, amino_acid_frequencies

                else:
                    print('There was an error scanning the pepstats file.')
                    print('Could not find corresponding sequence for identifier', TARGET_ID)
                    sys.exit(1)

    return pepstats_dic
# -----------------------------------------------------------------------------------------------------------
def predict_localization(input_data):
    seq, full_seq, RESULTS_PATH, OPTION, WEKA_PATH, PEPSTATS_PATH, SCRIPT_PATH = input_data[0], input_data[1], input_data[2], input_data[3], input_data[4], input_data[5], input_data[6]
    # -----------------------------------------------------------------------------------------------------------
    # Temporary folder name identifier that will be used to store results as multiple runs are done at the same time
    TMP = str(uuid.uuid4())
    TMP_PATH = RESULTS_PATH + 'pooledruns/' + TMP + '/'
    if not os.path.exists(TMP_PATH):
        os.makedirs(TMP_PATH)
    # -----------------------------------------------------------------------------------------------------------
    result_chloro_mito = []
    result_chloro = []
    result_mito = []
    result_chloro_mito_coordinates, result_chloro_coordinates, result_mito_coordinates = [], [], []
    # -----------------------------------------------------------------------------------------------------------
    pepstats_dic_aas = {}
    pepstats_dic_aas_short = {}
    # -----------------------------------------------------------------------------------------------------------
    if len(seq) >= 20.0:
        # If full sequences are used, only scan from the first amino acids as start position (e.g. in plant proteins)
        if OPTION == '-p':
            WINDOW_IDS, WINDOW_SEQS = functions.sliding_windows_plantmode(seq, parameters.WINDOW_SIZE, parameters.STEP_SIZE)              
            if WINDOW_IDS:
                pepstats_dic_aas = pepstats_aas(WINDOW_IDS, WINDOW_SEQS, TMP_PATH, PEPSTATS_PATH)
                # Now get the data for the first 15 aas
                WINDOW_SEQS_SHORT = [sequence[:15] for sequence in WINDOW_SEQS]
                pepstats_dic_aas_short = pepstats_aas(WINDOW_IDS, WINDOW_SEQS_SHORT, TMP_PATH, PEPSTATS_PATH)      
                
        if OPTION == '-e': 
            WINDOW_IDS, WINDOW_SEQS = functions.sliding_windows(seq, parameters.WINDOW_SIZE, parameters.STEP_SIZE)
            
            if WINDOW_IDS:
                pepstats_dic_aas = pepstats_aas(WINDOW_IDS, WINDOW_SEQS, TMP_PATH, PEPSTATS_PATH)
                # Now get the data for the first 15 aas
                WINDOW_SEQS_SHORT = [sequence[:15] for sequence in WINDOW_SEQS]
                pepstats_dic_aas_short = pepstats_aas(WINDOW_IDS, WINDOW_SEQS_SHORT, TMP_PATH, PEPSTATS_PATH)            
        # -----------------------------------------------------------------------------------------------------------
        boxes_predicted = []
        boxes_predicted_chloro = []
        boxes_predicted_mito = []
        # -----------------------------------------------------------------------------------------------------------
        predicted_chloro = chloro_classifier_allwindows(pepstats_dic_aas, pepstats_dic_aas_short, WINDOW_IDS, WINDOW_SEQS, TMP_PATH, WEKA_PATH, SCRIPT_PATH)  

        if OPTION == '-e' and len(predicted_chloro) <= 5:
            predicted_chloro = []
        
        for item in predicted_chloro:
            prob = float(item[1])
            window_id = item[0]
            coords = window_id.split('_')
            coords_start = coords[1]
            coords_end = coords[2].strip()        
            if prob > 0.6:   
                boxes_predicted_chloro.append((float(coords_start), float(coords_end), prob, 0))    
        # -----------------------------------------------------------------------------------------------------------
        predicted_mito = mito_classifier_allwindows(pepstats_dic_aas, pepstats_dic_aas_short, WINDOW_IDS, WINDOW_SEQS, TMP_PATH, WEKA_PATH, SCRIPT_PATH)                    

        if OPTION == '-e' and len(predicted_mito) <= 5:
            predicted_mito = []
            
        for item in predicted_mito:
            prob = float(item[1])
            window_id = item[0]
            coords = window_id.split('_')
            coords_start = coords[1]
            coords_end = coords[2].strip()       
            if prob > 0.6:    
                boxes_predicted_mito.append((float(coords_start), float(coords_end), prob, 1))         
        # -----------------------------------------------------------------------------------------------------------
        boxes_predicted += boxes_predicted_chloro
        boxes_predicted += boxes_predicted_mito

        # Pick the interval for all chloro and mito predictions with highest probability
        if boxes_predicted: 
            result_chloro_mito = [max(boxes_predicted, key=itemgetter(2))]

        if boxes_predicted_chloro:
            result_chloro = [max(boxes_predicted_chloro, key=itemgetter(2))]

        if boxes_predicted_mito:
            result_mito = [max(boxes_predicted_mito, key=itemgetter(2))]

        # These are the results for potential dual-localization
        result_chloro = [interval for interval in result_chloro if interval not in result_chloro_mito]
        result_mito = [interval for interval in result_mito if interval not in result_chloro_mito]

        # Now map the coordinates from the mature sequence to the full sequence
        result_chloro_coordinates = map_coordinates(result_chloro, seq, full_seq)
        result_mito_coordinates = map_coordinates(result_mito, seq, full_seq)
        result_chloro_mito_coordinates = map_coordinates(result_chloro_mito, seq, full_seq)    
        # -----------------------------------------------------------------------------------------------------------
        # Now try and determine the cleavage site by re-scanning the predicted 
        # transit peptide with smaller windows and step size 1       
        boxes_predicted = []
        
        if result_chloro_mito_coordinates:

            transit_peptide_seq = seq[int(result_chloro_mito[0][0]): int(result_chloro_mito[0][1]) + 1]

            if OPTION == '-e':
                WINDOW_IDS, WINDOW_SEQS = functions.sliding_windows_cleavage(transit_peptide_seq, parameters.WINDOW_SIZE_CLEAVAGE, parameters.STEP_SIZE)        
            if OPTION == '-p':
                WINDOW_IDS, WINDOW_SEQS = functions.sliding_windows_plantmode_cleavage(transit_peptide_seq, parameters.WINDOW_SIZE_CLEAVAGE, parameters.STEP_SIZE)  

            if WINDOW_IDS:
                pepstats_dic_aas = pepstats_aas(WINDOW_IDS, WINDOW_SEQS, TMP_PATH, PEPSTATS_PATH)
                # Now get the data for the first 15 aas
                WINDOW_SEQS_SHORT = [sequence[:15] for sequence in WINDOW_SEQS]
                pepstats_dic_aas_short = pepstats_aas(WINDOW_IDS, WINDOW_SEQS_SHORT, TMP_PATH, PEPSTATS_PATH) 

                # If the region was a cTP
                if result_chloro_mito_coordinates[3] == 0:     
                    prediction = chloro_classifier_allwindows(pepstats_dic_aas, pepstats_dic_aas_short, WINDOW_IDS, WINDOW_SEQS, TMP_PATH, WEKA_PATH, SCRIPT_PATH)  
                # If the region was a mTP
                if result_chloro_mito_coordinates[3] == 1:     
                    prediction = mito_classifier_allwindows(pepstats_dic_aas, pepstats_dic_aas_short, WINDOW_IDS, WINDOW_SEQS, TMP_PATH, WEKA_PATH, SCRIPT_PATH)                    

                for item in prediction:
                    prob = item[1]
                    window_id = item[0]
                    coords = window_id.split('_')
                    coords_start = coords[1]
                    coords_end = coords[2].strip()           
                    if result_chloro_mito_coordinates[3] == 0:     
                        boxes_predicted.append((float(coords_start), float(coords_end), prob, 0))   
                    if result_chloro_mito_coordinates[3] == 1:     
                        boxes_predicted.append((float(coords_start), float(coords_end), prob, 1))   

                if boxes_predicted: 
                    result_chloro_mito = [max(boxes_predicted, key=itemgetter(2))]

                # If the probability of a different window is higher
                if result_chloro_mito[0][2] > result_chloro_mito_coordinates[2]:
                    result_chloro_mito_coordinates = [int(result_chloro_mito_coordinates[0] + result_chloro_mito[0][0]), int(result_chloro_mito_coordinates[0] + result_chloro_mito[0][1]), result_chloro_mito[0][2], result_chloro_mito[0][3]]
        # -----------------------------------------------------------------------------------------------------------
        # Now try and determine the cleavage site of the dual localization by re-scanning 
        # the predicted transit peptide with smaller windows and step size 1  
        boxes_predicted = []

        if result_chloro_coordinates:

            transit_peptide_seq = seq[int(result_chloro[0][0]): int(result_chloro[0][1]) + 1]

            if OPTION == '-e':
                WINDOW_IDS, WINDOW_SEQS = functions.sliding_windows_cleavage(transit_peptide_seq, parameters.WINDOW_SIZE_CLEAVAGE, parameters.STEP_SIZE)        
            if OPTION == '-p':
                WINDOW_IDS, WINDOW_SEQS = functions.sliding_windows_plantmode_cleavage(transit_peptide_seq, parameters.WINDOW_SIZE_CLEAVAGE, parameters.STEP_SIZE)  

            if WINDOW_IDS:
                pepstats_dic_aas = pepstats_aas(WINDOW_IDS, WINDOW_SEQS, TMP_PATH, PEPSTATS_PATH)
                # Now get the data for the first 15 aas
                WINDOW_SEQS_SHORT = [sequence[:15] for sequence in WINDOW_SEQS]
                pepstats_dic_aas_short = pepstats_aas(WINDOW_IDS, WINDOW_SEQS_SHORT, TMP_PATH, PEPSTATS_PATH)      

                predicted_chloro = chloro_classifier_allwindows(pepstats_dic_aas, pepstats_dic_aas_short, WINDOW_IDS, WINDOW_SEQS, TMP_PATH, WEKA_PATH, SCRIPT_PATH)  

                for item in predicted_chloro:
                    prob = item[1]
                    window_id = item[0]
                    coords = window_id.split('_')
                    coords_start = coords[1]
                    coords_end = coords[2].strip()           
                    boxes_predicted.append((float(coords_start), float(coords_end), prob, 0))  
     
                if boxes_predicted: 
                    result_chloro = [max(boxes_predicted, key=itemgetter(2))]

                # If the probability of a different window is higher
                if result_chloro[0][2] > result_chloro_coordinates[2]:
                    result_chloro_coordinates = [int(result_chloro_coordinates[0] + result_chloro[0][0]), int(result_chloro_coordinates[0] + result_chloro[0][1]), result_chloro[0][2], result_chloro[0][3]]
        # -----------------------------------------------------------------------------------------------------------
        # Now try and determine the cleavage site of the dual localization by re-scanning 
        # the predicted transit peptide with smaller windows and step size 1 
        boxes_predicted = []

        if result_mito_coordinates:

            transit_peptide_seq = seq[int(result_mito[0][0]): int(result_mito[0][1]) + 1]

            if OPTION == '-e':
                WINDOW_IDS, WINDOW_SEQS = functions.sliding_windows_cleavage(transit_peptide_seq, parameters.WINDOW_SIZE_CLEAVAGE, parameters.STEP_SIZE)        
            if OPTION == '-p':
                WINDOW_IDS, WINDOW_SEQS = functions.sliding_windows_plantmode_cleavage(transit_peptide_seq, parameters.WINDOW_SIZE_CLEAVAGE, parameters.STEP_SIZE)  

            if WINDOW_IDS:
                pepstats_dic_aas = pepstats_aas(WINDOW_IDS, WINDOW_SEQS, TMP_PATH, PEPSTATS_PATH)
                # Now get the data for the first 15 aas
                WINDOW_SEQS_SHORT = [sequence[:15] for sequence in WINDOW_SEQS]
                pepstats_dic_aas_short = pepstats_aas(WINDOW_IDS, WINDOW_SEQS_SHORT, TMP_PATH, PEPSTATS_PATH)
          
                predicted_mito = mito_classifier_allwindows(pepstats_dic_aas, pepstats_dic_aas_short, WINDOW_IDS, WINDOW_SEQS, TMP_PATH, WEKA_PATH, SCRIPT_PATH)  
                
                for item in predicted_mito:
                    prob = item[1]
                    window_id = item[0]
                    coords = window_id.split('_')
                    coords_start = coords[1]
                    coords_end = coords[2].strip()           
                    boxes_predicted.append((float(coords_start), float(coords_end), prob, 1))   

                if boxes_predicted:     
                    result_mito = [max(boxes_predicted, key=itemgetter(2))]

                # If the probability of a different window is higher
                if result_mito[0][2] > result_mito_coordinates[2]:
                    result_mito_coordinates = [int(result_mito_coordinates[0] + result_mito[0][0]), int(result_mito_coordinates[0] + result_mito[0][1]), result_mito[0][2], result_mito[0][3]]        
        # -----------------------------------------------------------------------------------------------------------
        # After re-scanning the intervals there could be the possibility that the dual 
        # localization has switched from e.g. cTP/possible mTP to mTP/possible cTP

        if result_chloro_mito_coordinates and result_mito_coordinates:
            # Check that the order is correct for cTP/possible mTP and change if necessary
            if result_chloro_mito_coordinates[2] < result_mito_coordinates[2]:
                # cTP becomes possible cTP
                result_chloro_coordinates = [result_chloro_mito_coordinates[0], result_chloro_mito_coordinates[1], result_chloro_mito_coordinates[2], result_chloro_mito_coordinates[3]]
                # possible mTP becomes mTP
                result_chloro_mito_coordinates = [result_mito_coordinates[0], result_mito_coordinates[1], result_mito_coordinates[2], result_mito_coordinates[3]]

        if result_chloro_mito_coordinates and result_chloro_coordinates:
            # Check that the order is correct for mTP/possible cTP and change if necessary
            if result_chloro_mito_coordinates[2] < result_chloro_coordinates[2]:
                # mTP becomes possible mTP
                result_mito_coordinates = [result_chloro_mito_coordinates[0], result_chloro_mito_coordinates[1], result_chloro_mito_coordinates[2], result_chloro_mito_coordinates[3]]
                # possible cTP becomes cTP
                result_chloro_mito_coordinates = [result_chloro_coordinates[0], result_chloro_coordinates[1], result_chloro_coordinates[2], result_chloro_coordinates[3]]
        # -----------------------------------------------------------------------------------------------------------
        # If the input sequence is shorter than the predicted transit peptides 
        if result_chloro_mito_coordinates:
            start, end, prob, classification = result_chloro_mito_coordinates[0], result_chloro_mito_coordinates[1], result_chloro_mito_coordinates[2], result_chloro_mito_coordinates[3]
            if end > len(seq):
                result_chloro_mito_coordinates = [start, int(len(seq) - 1.0), prob, classification]

        if result_chloro_coordinates:
            start, end, prob, classification = result_chloro_coordinates[0], result_chloro_coordinates[1], result_chloro_coordinates[2], result_chloro_coordinates[3]
            if end > len(seq):
                result_chloro_coordinates = [start, int(len(seq) - 1.0), prob, classification]
        
        if result_mito_coordinates:
            start, end, prob, classification = result_mito_coordinates[0], result_mito_coordinates[1], result_mito_coordinates[2], result_mito_coordinates[3]
            if end > len(seq):
                result_mito_coordinates = [start, int(len(seq) - 1.0), prob, classification]
    # -----------------------------------------------------------------------------------------------------------
    # Search for bipartite NLS
    bipart_motif = bipart_NLS(seq)
    # Search for NLS
    NLS_motifs = predict_NLS(seq)
    # Call NLSstradamus
    NLStradamus_result = NLStradamus(seq, TMP_PATH, SCRIPT_PATH)
    # -----------------------------------------------------------------------------------------------------------
    # Clean up and delete temporary folder that was created
    shutil.rmtree(TMP_PATH)
    # -----------------------------------------------------------------------------------------------------------
    return NLS_motifs, bipart_motif, NLStradamus_result, result_chloro_mito_coordinates, result_chloro_coordinates, result_mito_coordinates
# -----------------------------------------------------------------------------------------------------------

