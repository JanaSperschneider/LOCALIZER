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
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import errno
import uuid
import shutil
import functions
import output
import localization
import multiprocessing
from multiprocessing import Pool
import tempfile
# -----------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Main Program starts here
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':    
    SCRIPT_PATH = sys.path[0]
    # -----------------------------------------------------------------------------------------------------------
    # Change the path to WEKA to the appropriate location on your computer
    WEKA_PATH = SCRIPT_PATH + '/weka-3-6-12/weka.jar'
    PEPSTATS_PATH = SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'
    # -----------------------------------------------------------------------------------------------------------
    # Check that the path to the WEKA software exists
    path_exists = os.access(WEKA_PATH, os.F_OK)
    if not path_exists:
        print()
        print("Path to WEKA software does not exist!")
        print("Check the installation and the given path to the WEKA software %s in LOCALIZER.py (line 62)." % WEKA_PATH)
        print()
        sys.exit(1)
    # -----------------------------------------------------------------------------------------------------------
    # Check that the path to the EMBOSS software exists for pepstats
    path_exists = os.access(PEPSTATS_PATH, os.F_OK)
    if not path_exists:
        print()
        print("Path to EMBOSS software does not exist!")
        print("Check the installation and the given path to the EMBOSS software %s in LOCALIZER.py (line 63)." % PEPSTATS_PATH)
        print()
        sys.exit(1)
    # -----------------------------------------------------------------------------------------------------------
    commandline = sys.argv[1:]
    # -----------------------------------------------------------------------------------------------------------
    if commandline:
        FASTA_FILE, OPTION, noSP_option, SP_length, output_folder = functions.scan_arguments(commandline)
        # If no FASTA file was provided with the -i option
        if not FASTA_FILE:
            print()
            print('Please specify a FASTA input file using the -i option!')
            functions.usage()
        if not OPTION:
            print()
            print('Please specify either plant (-p) or effector (-e) mode!')
            functions.usage()
    else:
        functions.usage()
    # -----------------------------------------------------------------------------------------------------------
    # Temporary folder 
    RESULTS_PATH = tempfile.mkdtemp() + '/'
    # -----------------------------------------------------------------------------------------------------------
    # Output folder if desired by user
    if output_folder:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
    # -----------------------------------------------------------------------------------------------------------
    # Check if FASTA file exists
    try:
        open(FASTA_FILE, 'r') 
    except OSError as e:
        print("Unable to open FASTA file:", FASTA_FILE)  #Does not exist OR no read permissions
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit(1)
    # -----------------------------------------------------------------------------------------------------------
    # Try to create folder where results will be stored
    try:
        os.mkdir(RESULTS_PATH)
    except OSError as exception:        
        if exception.errno != errno.EEXIST:
            raise
    # -----------------------------------------------------------------------------------------------------------
    ORIGINAL_IDENTIFIERS, SEQUENCES = functions.get_seqs_ids_fasta(FASTA_FILE)
    SEQUENCES = [seq.upper() for seq in SEQUENCES]
    PROTEINS_CLASSIFIED = len(ORIGINAL_IDENTIFIERS)
    # -----------------------------------------------------------------------------------------------------------
    # Check that all sequences are >= 40 aas
    '''for ident, seq in zip(ORIGINAL_IDENTIFIERS, SEQUENCES):   
        if len(seq.strip()) < 40.0:
            print 'Please provide sequences >= 40 aas'
            print 'This sequence caused problems:', ident.strip(), ' with length:', len(seq.strip()) 
            shutil.rmtree(RESULTS_PATH)
            sys.exit()'''
    # -----------------------------------------------------------------------------------------------------------
    # Replace ambiguous amino acids and return a list of identfiers that should not be used
    # later: those with more than 10% unknown bases in sequence 
    ORIGINAL_SEQUENCES = SEQUENCES
    SEQUENCES, BLACKLISTED = functions.filterX(ORIGINAL_IDENTIFIERS, SEQUENCES)
    # -----------------------------------------------------------------------------------------------------------
    # Write new FASTA file with short identifiers for SignalP and pepstats
    f_output = RESULTS_PATH + '_short_ids.fasta'
    SHORT_IDENTIFIERS = functions.write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, SEQUENCES)
    # -----------------------------------------------------------------------------------------------------------
    if OPTION == '-e' and SP_length == None:    
        SP_length = 20
    # Eukaryotic effector mode
    if OPTION == '-e':
        # Produce mature sequences
        MATURE_SEQUENCES = [seq[SP_length:] for seq in SEQUENCES]
    # Plant mode
    if OPTION == '-p' or noSP_option:
        # Use full sequences
        MATURE_SEQUENCES = SEQUENCES
    # -----------------------------------------------------------------------------------------------------------
    input_list = [[seq, full_seq, RESULTS_PATH, OPTION, WEKA_PATH, PEPSTATS_PATH, SCRIPT_PATH] for ident, seq, full_seq in zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES)]

    localizations = []

    printed10, printed20, printed30, printed40, printed50 = False, False, False, False, False
    printed60, printed70, printed80, printed90 = False, False, False, False
    
    for i, (seq, full_seq, RESULTS_PATH, OPTION, WEKA_PATH, PEPSTATS_PATH, SCRIPT_PATH) in enumerate(input_list):
        # The format is: NLS_motifs, bipart_motif, NLStradamus_result, result_chloro_mito        
        localizations.append(localization.predict_localization([seq, full_seq, RESULTS_PATH, OPTION, WEKA_PATH, PEPSTATS_PATH, SCRIPT_PATH]))
        if float(len(input_list)) > 10.0:
            if i/float(len(input_list)) > 0.10 and not printed10:
                print('Over 10% are done...')
                printed10 = True
            if i/float(len(input_list)) > 0.20 and not printed20:
                print('Over 20% are done...')
                printed20 = True
            if i/float(len(input_list)) > 0.30 and not printed30:
                print('Over 30% are done...')
                printed30 = True
            if i/float(len(input_list)) > 0.40 and not printed40:
                print('Over 40% are done...')
                printed40 = True
            if i/float(len(input_list)) > 0.50 and not printed50:
                print('Over 50% are done...')
                printed50 = True
            if i/float(len(input_list)) > 0.60 and not printed60:
                print('Over 60% are done...')
                printed60 = True
            if i/float(len(input_list)) > 0.70 and not printed70:
                print('Over 70% are done...')
                printed70 = True
            if i/float(len(input_list)) > 0.80 and not printed80:
                print('Over 80% are done...')
                printed80 = True
            if i/float(len(input_list)) > 0.90 and not printed90:
                print('Over 90% are done...')
                printed90 = True  

    #try:
    #    p = Pool(multiprocessing.cpu_count())
    #except:
    #    p = Pool(2)    
    #localizations = p.map(localization.predict_localization, input_list)

    # Do not use sequences with more than 10% unknown bases
    for index, localization in enumerate(localizations):
        if ORIGINAL_IDENTIFIERS[index] in BLACKLISTED:
            localizations[index] = [], [], [], [], [], []             
    # -----------------------------------------------------------------------------------------------------------
    # Write output FASTA files if user wants this    
    count_chloro_only = output.chloroplast_only(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)
    count_chloro_mito = output.chloroplast_mitochondria(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)
    count_chloro_mito_nucleus = output.chloroplast_mitochondria_nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)
    count_chloro_nucleus = output.chloroplast_nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)

    count_mito_only = output.mitochondria_only(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)
    count_mito_chloro = output.mitochondria_chloroplast(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)
    count_mito_chloro_nucleus = output.mitochondria_chloroplast_nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)
    count_mito_nucleus = output.mitochondria_nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)

    count_nucleus_only = output.nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, ORIGINAL_SEQUENCES, localizations, output_folder)

    if output_folder:
        for each_file in os.listdir(output_folder):
            # check size and delete if 0
            if os.path.getsize(output_folder + '/' + each_file) == 0:
                os.remove(output_folder + '/' + each_file)
    # -----------------------------------------------------------------------------------------------------------
    output_screen, output_file = output.result_table(localizations, ORIGINAL_IDENTIFIERS, OPTION)
    # -----------------------------------------------------------------------------------------------------------
    summary = '# Proteins analyzed: ' + str(PROTEINS_CLASSIFIED) + ' from file: ' + FASTA_FILE + '\n\n'
    summary += '# Number of proteins with cTP: ' + str(int(count_chloro_only))
    summary += ' (' + str(round(100.0*count_chloro_only/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary += '# Number of proteins with cTP & possible mTP: ' + str(int(count_chloro_mito))
    summary += ' (' + str(round(100.0*count_chloro_mito/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary += '# Number of proteins with cTP & NLS: ' + str(int(count_chloro_nucleus))
    summary += ' (' + str(round(100.0*count_chloro_nucleus/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary += '# Number of proteins with cTP & possible mTP & NLS: ' + str(int(count_chloro_mito_nucleus))
    summary += ' (' + str(round(100.0*count_chloro_mito_nucleus/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'

    summary += '# Number of proteins with mTP: ' + str(int(count_mito_only))
    summary += ' (' + str(round(100.0*count_mito_only/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary += '# Number of proteins with mTP & possible cTP: ' + str(int(count_mito_chloro))
    summary += ' (' + str(round(100.0*count_mito_chloro/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary += '# Number of proteins with mTP & NLS: ' + str(int(count_mito_nucleus))
    summary += ' (' + str(round(100.0*count_mito_nucleus/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary += '# Number of proteins with mTP & possible cTP & NLS: ' + str(int(count_mito_chloro_nucleus))
    summary += ' (' + str(round(100.0*count_mito_chloro_nucleus/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary += '# Number of proteins with NLS and no transit peptides: '+ str(int(count_nucleus_only))
    summary += ' (' + str(round(100.0*count_nucleus_only/PROTEINS_CLASSIFIED,1)) + '%)'
    # -----------------------------------------------------------------------------------------------------------
    summary_collapsed = '# Summary statistics\n\n'
    summary_collapsed += '# Number of proteins with chloroplast localization (cTP, cTP & possible mTP, cTP & NLS, cTP & possible mTP & NLS): ' + str(int(count_chloro_only) + int(count_chloro_mito) + int(count_chloro_nucleus) + int(count_chloro_mito_nucleus))
    summary_collapsed += ' (' + str(round(100.0*(int(count_chloro_only) + int(count_chloro_mito) + int(count_chloro_nucleus) + int(count_chloro_mito_nucleus))/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary_collapsed += '# Number of proteins with mitochondrial localization (mTP, mTP & possible cTP, mTP & NLS, mTP & possible cTP & NLS): ' + str(int(count_mito_only) + int(count_mito_chloro) + int(count_mito_nucleus) + int(count_mito_chloro_nucleus))
    summary_collapsed += ' (' + str(round(100.0*(int(count_mito_only) + int(count_mito_chloro) + int(count_mito_nucleus) + int(count_mito_chloro_nucleus))/PROTEINS_CLASSIFIED,1)) + '%)' + '\n'
    summary_collapsed += '# Number of proteins with nuclear localization and no transit peptides: ' + str(int(count_nucleus_only))
    summary_collapsed += ' (' + str(round(100.0*(int(count_nucleus_only))/PROTEINS_CLASSIFIED,1)) + '%)'  + '\n'
    summary_collapsed += '# Number of proteins with nuclear localization and with transit peptides: ' + str(int(count_chloro_nucleus) + int(count_chloro_mito_nucleus) + int(count_mito_nucleus) + int(count_mito_chloro_nucleus))
    summary_collapsed += ' (' + str(round(100.0*(int(count_chloro_nucleus) + int(count_chloro_mito_nucleus) + int(count_mito_nucleus) + int(count_mito_chloro_nucleus))/PROTEINS_CLASSIFIED,1)) + '%)'
    # -----------------------------------------------------------------------------------------------------------   
    print(output_screen)
    print('--------------------------------------')
    print('--------------------------------------')
    print(summary)
    print('--------------------------------------')
    print('--------------------------------------')
    print(summary_collapsed)
    print('--------------------------------------')
    print('--------------------------------------')
    # -----------------------------------------------------------------------------------------------------------   
    if output_folder:
        f_output = output_folder + '/Results.txt'
        f = open(f_output, 'w')
        f.writelines(output_file + '\n')
        f.writelines('#--------------------------------------\n')
        f.writelines('#--------------------------------------\n')
        f.writelines(summary + '\n')
        f.writelines('#--------------------------------------\n')
        f.writelines('#--------------------------------------\n')
        f.writelines(summary_collapsed + '\n')
        f.writelines('#--------------------------------------\n')
        f.writelines('#--------------------------------------\n')
        f.close()
    # -----------------------------------------------------------------------------------------------------------
    # Clean up and delete temporary folder that was created
    shutil.rmtree(RESULTS_PATH)
    # -----------------------------------------------------------------------------------------------------------
