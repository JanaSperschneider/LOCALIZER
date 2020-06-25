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
import os
import sys
import textwrap
# -----------------------------------------------------------------------------------------------------------
def result_table(localizations, ORIGINAL_IDENTIFIERS, OPTION):
    """ Function: result_table()

        Purpose:  Given the predicted localizations and identifiers,  
                  write the output table as a string.
              
        Input:    Predicted localizations and identifiers.                  
    
        Return:   String that contains list of predicted localizations with probabilities
                  and the longest NLSs that were predicted.
    """
    # Identify the longest ID for padding the columns
    # Split the identifier to include only a short version
    #padding = max([len(ident.replace('>','').split()[0].strip()) for ident in ORIGINAL_IDENTIFIERS])
    padding = max([len(ident.replace('>','').strip()) for ident in ORIGINAL_IDENTIFIERS])
    if padding < len('Identifier'):
        padding = len('Identifier')

    padding_screen = padding
    if padding_screen > 70:
        padding_screen = 70

    # Print a result table
    output_screen =  '# -----------------' + '\n'
    output_screen += '# LOCALIZER 1.0.5 Predictions (' + str(OPTION) + ' mode) \n'
    output_screen += '# -----------------' + '\n'
    output_screen += 'Identifier'.ljust(padding_screen) + '\t' + 'Chloroplast'.ljust(18) + '\t' + 'Mitochondria'.ljust(18) + '\t' + 'Nucleus' + '\n'

    output_file =  '# -----------------' + '\n'
    output_file += '# LOCALIZER 1.0.5 Predictions (' + str(OPTION) + ' mode) \n'
    output_file += '# -----------------' + '\n'
    output_file += 'Identifier'.ljust(padding) + '\t' + 'Chloroplast'.ljust(18) + '\t' + 'Mitochondria'.ljust(18) + '\t' + 'Nucleus' + '\n'


    for index, prediction in enumerate(localizations):
        NLS_motifs, bipart_motif, NLStradamus_result = prediction[0], prediction[1], prediction[2]
        result_coordinates = prediction[3]

        additional_chloro, additional_mito = prediction[4], prediction[5]

        nucleus_output = '-'.ljust(7)
        chloro_output = '-'.ljust(18)
        mito_output = '-'.ljust(18)   

        if NLS_motifs or bipart_motif or NLStradamus_result:
            nucleus_motifs = []
            
            if NLS_motifs:
                for motif in NLS_motifs:
                    nucleus_motifs.append(motif)
            if bipart_motif:
                for motif in bipart_motif:
                    nucleus_motifs.append(motif)
            if NLStradamus_result:
                for motif in NLStradamus_result:
                    nucleus_motifs.append(motif[3])

            nucleus_motifs.sort(key=len)

            if len(nucleus_motifs) > 1.0:
                nucleus_motifs_reduced = []
                for pos, NLS in enumerate(nucleus_motifs):
                    contained = False
                    for other_NLSs in nucleus_motifs[pos+1:]:
                        if NLS in other_NLSs:
                            contained = True
                    if not contained:
                        nucleus_motifs_reduced.append(NLS)
            else:
                nucleus_motifs_reduced = nucleus_motifs

            nucleus_output = 'Y' + ' (' 
            for motif in nucleus_motifs_reduced:
                nucleus_output += motif + ','

            nucleus_output = nucleus_output[:-1]
            nucleus_output += ')'                

        # Print any chloroplast localization
        if result_coordinates:        
            # This is chloroplast as the main prediction
            if result_coordinates[3] == 0:
                chloro_output = 'Y' + ' (' + str(result_coordinates[2]) + ' | ' + str(result_coordinates[0] + 1) + '-' + str(result_coordinates[1] + 1) + ')'

        if additional_mito:
            if additional_mito[3] == 1:
                mito_output = 'Y' + ' (' + str(additional_mito[2]) + ' | ' + str(additional_mito[0] + 1) + '-' + str(additional_mito[1] + 1) + ')'

        # Print any mitochondrial localization
        if result_coordinates:        
            # This is mitochondria as the main prediction
            if result_coordinates[3] == 1:
                mito_output = 'Y' + ' (' + str(result_coordinates[2]) + ' | ' + str(result_coordinates[0] + 1) + '-' + str(result_coordinates[1] + 1) + ')'

        if additional_chloro:
            if additional_chloro[3] == 0:
                chloro_output = 'Y' + ' (' + str(additional_chloro[2]) + ' | ' + str(additional_chloro[0] + 1) + '-' + str(additional_chloro[1] + 1) + ')'

        identifier_output = ORIGINAL_IDENTIFIERS[index]
        identifier_output = identifier_output.replace('>','')
        # Split the identifier to include only a short version
        #identifier_output = identifier_output.split()[0]
        identifier_output = identifier_output.strip()
        identifier_output = identifier_output.ljust(padding)    

        # This is the output line for a particular protein
        output_file += identifier_output + '\t' + chloro_output.ljust(18) + '\t' + mito_output.ljust(18) + '\t' + nucleus_output + '\n'
           
        new_identifier = textwrap.wrap(identifier_output)[:-1] + [textwrap.wrap(identifier_output)[-1].ljust(padding_screen)]
        new_identifier = "\n".join(new_identifier)
        output_screen += new_identifier + '\t' + chloro_output.ljust(18) + '\t' + mito_output.ljust(18) + '\t' + nucleus_output + '\n'

    output_screen = output_screen.rstrip()

    return output_screen, output_file
# -----------------------------------------------------------------------------------------------------------
def scan_localization(localizations, index):

    transit_peptide = localizations[index][3]
    additional_result_chloro, additional_result_mito = localizations[index][4], localizations[index][5]

    if transit_peptide:        
        start, end, prob, loc = transit_peptide[0], transit_peptide[1], transit_peptide[2], transit_peptide[3]
    else:
        start, end, prob, loc = None, None, None, None

    if localizations[index][0] or localizations[index][1] or localizations[index][2]:
      	nucleus_found = True        
    else:
        nucleus_found = False

    return start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found  
# -----------------------------------------------------------------------------------------------------------
def chloroplast_only(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/chloroplast_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins exclusively targeted to the chloroplast
        if loc == 0.0 and not nucleus_found and not additional_result_mito:
            if output_folder:            
                f.writelines(ident.strip() + ' | cTP coordinates: ' + str(start) + ':' + str(end) + ' | Probability: ' + str(prob) + '\n')
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------
def chloroplast_mitochondria(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/chloroplast_mitochondria_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins exclusively targeted to the chloroplast
        if loc == 0.0 and not nucleus_found and additional_result_mito:
            ident = ident.strip() + ' | cTP coordinates: ' + str(start) + ':' + str(end) + ' | Probability: ' + str(prob)
            ident += ' | Possible mTP coordinates: ' + str(additional_result_mito[0]) + ':' + str(additional_result_mito[1]) + ' | Probability: ' + str(additional_result_mito[2])
            if output_folder:
                f.writelines(ident + '\n')
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------
def chloroplast_mitochondria_nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/chloroplast_mitochondria_nucleus_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins exclusively targeted to the chloroplast
        if loc == 0.0 and nucleus_found and additional_result_mito:
            ident = ident.strip() + ' | cTP coordinates: ' + str(start) + ':' + str(end) + ' | Probability: ' + str(prob)
            ident += ' | Possible mTP coordinates: ' + str(additional_result_mito[0]) + ':' + str(additional_result_mito[1]) + ' | Probability: ' + str(additional_result_mito[2])
            if output_folder:
                f.writelines(ident + '\n')
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------
def chloroplast_nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/chloroplast_nucleus_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins exclusively targeted to the chloroplast
        if loc == 0.0 and nucleus_found and not additional_result_mito:
            ident = ident.strip() + ' | cTP coordinates: ' + str(start) + ':' + str(end) + ' | Probability: ' + str(prob)
            if output_folder:           
                f.writelines(ident + '\n')
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------
def mitochondria_only(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/mitochondria_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins exclusively targeted to mitochondria
        if loc == 1.0 and not nucleus_found and not additional_result_chloro:
            if output_folder:
                f.writelines(ident.strip() + ' | mTP coordinates: ' + str(start) + ':' + str(end) + ' | Probability: ' + str(prob) + '\n')
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------
def mitochondria_chloroplast(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/mitochondria_chloroplast_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins targeted to mitochondria and possibly chloroplast
        if loc == 1.0 and not nucleus_found and additional_result_chloro:
            ident = ident.strip() + ' | mTP coordinates: ' + str(start) + ':' + str(end) + ' | Probability: ' + str(prob)
            ident += ' | Possible cTP coordinates: ' + str(additional_result_chloro[0]) + ':' + str(additional_result_chloro[1]) + ' | Probability: ' + str(additional_result_chloro[2])
            if output_folder:            
                f.writelines(ident + '\n')
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------
def mitochondria_chloroplast_nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/mitochondria_chloroplast_nucleus_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins targeted to mitochondria and possibly chloroplast
        if loc == 1.0 and nucleus_found and additional_result_chloro:
            ident = ident.strip() + ' | mTP coordinates: ' + str(start) + ':' + str(end) + ' | Probability: ' + str(prob)
            ident += ' | Possible cTP coordinates: ' + str(additional_result_chloro[0]) + ':' + str(additional_result_chloro[1]) + ' | Probability: ' + str(additional_result_chloro[2])
            if output_folder:
                f.writelines(ident + '\n')
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------
def mitochondria_nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/mitochondria_nucleus_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins targeted to mitochondria and possibly chloroplast
        if loc == 1.0 and nucleus_found and not additional_result_chloro:
            ident = ident.strip() + ' | mTP coordinates: ' + str(start) + ':' + str(end) + ' | Probability: ' + str(prob)
            if output_folder:
                f.writelines(ident + '\n')
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------
def nucleus(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES, SEQUENCES, localizations, output_folder):

    count = 0.0

    if output_folder:
        f = open(output_folder + '/nucleus_predicted.fasta', 'w')

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, MATURE_SEQUENCES)):

        start, end, prob, loc, additional_result_chloro, additional_result_mito, nucleus_found = scan_localization(localizations, index)
        #--------------------------------------
        # Proteins targeted to nucleus
        if loc == None and nucleus_found and not loc and not additional_result_chloro and not additional_result_mito:
            if output_folder:
                f.writelines(ident)
                f.writelines(SEQUENCES[index] + '\n')
            count += 1.0
        #--------------------------------------
    if output_folder:
        f.close()

    return count
# -----------------------------------------------------------------------------------------------------------


