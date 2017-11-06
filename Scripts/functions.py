#! /usr/bin/python
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
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import io
import getopt
import random
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def usage():
    """ Function: usage()

        Purpose:  Print helpful information for the user.        
        
        Input:    None.
    
        Return:   Print options for running LOCALIZER.       
    """
    print '''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LOCALIZER :: predicting subcellular localization of plant and eukaryotic effector proteins
# LOCALIZER 1.0.3 (November 2017); http://localizer.csiro.au/
# Copyright (C) 2016-2017 Jana Sperschneider, CSIRO.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    '''
    print "Usage for LOCALIZER: ", 
    print "python LOCALIZER.py [-options] -i <input_file>"
    print 
    print
    print "required option:"
    print "    -p      : Run in plant mode. This means LOCALIZER will search for transit peptides starting at the first amino acid."
    print "  OR"
    print "    -e      : Run in effector mode. This means LOCALIZER will search for transit peptides after the signal peptide has been removed."
    print
    print "other options in effector mode:"    
    print "    -M      : in effector mode, do not remove the signal peptide. Use this if you are providing mature effector sequences."
    print "    -S <x>  : in effector mode, remove the signal peptide by deleting the first x aas (default: 20)."
    print
    print "other options:"    
    print "    -o <f>  : save a tab-separated output table (Results.txt) and proteins with predicted localizations to FASTA files in folder <f>"        
    print "    -h      : show brief help on version and usage"     
    print
    sys.exit()    

    return
# -----------------------------------------------------------------------------------------------------------
def scan_arguments(commandline):
    """ Function: scan_arguments()

        Purpose:  Scan the input options given to the LOCALIZER program.        
        
        Input:    Input options given by the user.
    
        Return:   Parsed options.
    """
    try:
        opts, args = getopt.getopt(commandline, "peMhS:o:i:", ["help"])        
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    OPTION = None
    FASTA_FILE = None
    noSP_option = False
    SP_length = None
    output_folder = None

    option_count, i_count, o_count, noSP_count, Slength_count = 0, 0, 0, 0, 0

    for opt, arg in opts:

        if opt in ("-i"):
            FASTA_FILE = arg
            i_count += 1

        elif opt in ("-o"):
            output_folder = arg
            o_count += 1

        elif opt in ("-p"):
            OPTION = opt
            option_count += 1

        elif opt in ("-e"):
            OPTION = opt
            option_count += 1

        elif opt in ("-M"):
            noSP_option = True
            noSP_count += 1

        elif opt in ("-S"):
            SP_length = arg
            Slength_count += 1
            try:
                SP_length = int(SP_length)
            except:
                usage()            

        elif opt in ("-h", "--help"):
            usage()

        else:
            assert False, "unhandled option"   

    if i_count > 1 or o_count > 1 or noSP_count > 1 or option_count > 1 or Slength_count > 1:
        usage()

    if SP_length and noSP_option:
        print '''# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'''
        print '''You have specified that mature sequences are provided with '-M' AND that the first ''', SP_length, '''aas should be deleted as the signal peptide region with '-S'.'''
        print
        print '''Please use only one of these arguments.'''
        print '''# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'''
        usage()    

    if SP_length and OPTION == "-p":
        print '''# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'''
        print '''You have specified that the first ''', SP_length, '''aas should be deleted as the signal peptide region with '-S'.'''
        print
        print '''Please use effector mode '-e' and not plant mode 'p'.'''
        print '''# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'''
        usage()    

    if noSP_option and OPTION == "-p":
        print '''# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'''
        print '''You have specified that mature effector sequences are provided with '-M'.'''
        print
        print '''Please use effector mode '-e' and not plant mode 'p'.'''
        print '''# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'''
        usage()    

    return FASTA_FILE, OPTION, noSP_option, SP_length, output_folder
# -----------------------------------------------------------------------------------------------------------
def get_seqs_ids_fasta(FASTA_FILE):
    """ Function: get_seqs_ids_fasta()

        Purpose:  Given a FASTA format file, this function extracts
                  the list of identifiers and the list of sequences 
                  in the order in which they appear in the FASTA file.
              
        Input:    Path to FASTA format file.
    
        Return:   List of identifiers and list of sequences in the order 
                  in which they appear in the FASTA file.
    """ 
    identifiers = []
    sequences = []

    with open(FASTA_FILE) as f: 
        content = f.readlines()

        for position, line in enumerate(content):
            if '>' in line:
                identifiers.append(line)
                seq = []
                following_lines = content[position + 1:]
                for next_line in following_lines:
                    if '>' not in next_line:
                        seq.append(next_line.strip())
                    else:
                        break    
                sequence = "".join(seq)
                sequence = sequence.replace('*', '')
                sequences.append(sequence)

    return identifiers, sequences
# -----------------------------------------------------------------------------------------------------------
def filterX(IDENTIFIERS, SEQUENCES):
    # Replace ambiguous amino acids because ProtParam can't deal with them later on 
    # B = Asx = Aspartic acid or Asparagine
    # Z = Glx = Glutamic acid or Glutamine
    # X = any amino acid
    replaceB = ['D', 'N']
    replaceZ = ['E', 'Q']
    replaceX = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    SEQUENCES_REPLACED = []
    BLACKLISTED = []

    for identifier, seq in zip(IDENTIFIERS, SEQUENCES):
        if len(seq):
            if seq.count('B')/float(len(seq)) >= 0.1 or seq.count('Z')/float(len(seq)) >= 0.1 or seq.count('X')/float(len(seq)) >= 0.1:
                BLACKLISTED.append(identifier)

    for sequence in SEQUENCES:

        if 'B' in sequence or 'Z' in sequence or 'X' in sequence:

            sequence_replaced = []
            for char in sequence:
                char = char.replace('B', random.choice(replaceB))
                char = char.replace('Z', random.choice(replaceZ))
                char = char.replace('X', random.choice(replaceX))
                sequence_replaced += char
            sequence_replaced = "".join(sequence_replaced)

            SEQUENCES_REPLACED.append(sequence_replaced.replace('*',''))
        else:
            SEQUENCES_REPLACED.append(sequence.replace('*',''))    

    return SEQUENCES_REPLACED, BLACKLISTED
# -----------------------------------------------------------------------------------------------------------
def write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES):
    """ Function: write_FASTA_short_ids()

        Purpose:  Given a list of identifiers and the corresponding list 
                  of sequence, write these to a FASTA file using short
                  identifiers such as protein1, protein2, .... This is 
                  done because some programs like pepstats do not like 
                  long identifier names as input.
              
        Input:    Path to desired FASTA format output file, list of 
                  identifiers and list of corresponding sequences.
    
        Return:   List of short identifiers.
    """

    with open(f_output, 'w') as f:
        SHORT_IDENTIFIERS = []
        # Change identifiers to protein1, protein2, ...
        # and write to temporary file
        SET = zip(ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES)
    
        for index,  (identifier, sequence) in enumerate(SET):
            short_id = '>protein' + str(index)
            SHORT_IDENTIFIERS.append(short_id)
            f.writelines(short_id + '\n')
            f.writelines(sequence + '\n')

    return SHORT_IDENTIFIERS
# -----------------------------------------------------------------------------------------------------------
def sliding_windows_plantmode(seq, WINDOW_SIZES, STEP_SIZE):

    WINDOW_IDS, WINDOW_SEQS = [], []

    for WINDOW_SIZE in WINDOW_SIZES:

        window_seq = seq[0:WINDOW_SIZE]

        WINDOW_IDS.append('protein_' + str(0) + '_' + str((WINDOW_SIZE)) + '\n')
        WINDOW_SEQS.append(window_seq)  

    return WINDOW_IDS, WINDOW_SEQS
# -----------------------------------------------------------------------------------------------------------
def sliding_windows_plantmode_cleavage(seq, WINDOW_SIZES, STEP_SIZE):

    # The N-terminus will be scanned up to this start position of the sliding windows
    # i.e.
    # 1:30 
    # ...
    # 50:80
    WINDOW_IDS, WINDOW_SEQS = [], []   

    for WINDOW_SIZE in WINDOW_SIZES:

        window_seq = seq[0:WINDOW_SIZE]

        if len(window_seq) == WINDOW_SIZE:
            WINDOW_IDS.append('protein_' + str(0) + '_' + str((WINDOW_SIZE)) + '\n')
            WINDOW_SEQS.append(window_seq)  

    return WINDOW_IDS, WINDOW_SEQS
# -----------------------------------------------------------------------------------------------------------
def sliding_windows(seq, WINDOW_SIZES, STEP_SIZE):

    # The N-terminus will be scanned up to this start position of the sliding windows
    # i.e.
    # 1:30 
    # ...
    # 50:80
    WINDOW_IDS, WINDOW_SEQS = [], []   

    # Leave at least 40 aas in the C-terminus
    TransitPeptide_STOP = len(seq) - 40

    # The pro-domain can occur in the first 50 aas after the signal peptide
    PRODOMAIN_STOP = 50

    for WINDOW_SIZE in WINDOW_SIZES:

        for START in xrange(0, int(TransitPeptide_STOP) + 1, STEP_SIZE):

            window_seq = seq[START:START + WINDOW_SIZE]

            if START < PRODOMAIN_STOP and (START + WINDOW_SIZE) < TransitPeptide_STOP + 1 and len(window_seq) == WINDOW_SIZE:
                WINDOW_IDS.append('protein_' + str(START) + '_' + str((START + WINDOW_SIZE)) + '\n')
                WINDOW_SEQS.append(window_seq)

    return WINDOW_IDS, WINDOW_SEQS
# -----------------------------------------------------------------------------------------------------------
def sliding_windows_cleavage(seq, WINDOW_SIZES, STEP_SIZE):

    # The N-terminus will be scanned up to this start position of the sliding windows
    # i.e.
    # 1:30 
    # ...
    # 50:80
    WINDOW_IDS, WINDOW_SEQS = [], []   

    for WINDOW_SIZE in WINDOW_SIZES:

        for START in xrange(0, len(seq), STEP_SIZE):

            window_seq = seq[START:START + WINDOW_SIZE]

            if len(window_seq) == WINDOW_SIZE:
                WINDOW_IDS.append('protein_' + str(START) + '_' + str((START + WINDOW_SIZE)) + '\n')
                WINDOW_SEQS.append(window_seq)

    return WINDOW_IDS, WINDOW_SEQS
# -----------------------------------------------------------------------------------------------------------

