"""
Author: Iain P. S. Smith
Date: 14:46 Tuesday 09/07/19

Overview:
---------
    This script takes in a structure data file with degenerate residue numbers (i.e more than one residue with the same index)
    and updates the residue numbers in the file with an incremental counter.

    This counter starts at i = 1 and updates i += 1 every time it registers a change in the residue numbers in the original file.
    Once this check has been done it then sets the residue number on the current line to be equal to the value of the counter, i.

    This assumes that the only problem in the original file was the degenerate residue numbers, and that all atoms within a 
    contiguous block of equal residue numbers are indeed grouped correctly (i.e by the monomer unit of the system).

    The code then continues this process over all solvent and ion groups to ensure that their residue numbers are updated accordingly.


    THIS CODE ALSO WORKS FOR .PDB FILES - However the output file is not named correctly:
        
        currently the filename.strip function removes first letter from .pdb source file, not sure how to fix this
"""


import sys
from tqdm import tqdm
import numpy as np

from copy import deepcopy

arglist = sys.argv

filename = arglist[1]

if '.pqr' in filename:
    filename_split = filename.strip('.pqr')
    filetype = ".pqr"
    

elif '.pdb' in filename:
    filename_split = filename.strip(".pdb")
    filetype = ".pdb"

def get_lines(filename):
    """
    Performs a read lines method on a file of title "filename.txt" and returns
    the result.
    
    Parameters
    ----------
    filename : string
        name of file to be read
        
    Output : list of strings
        list containing strings of each line from file
    """
    try:
        type(filename) == str
    except TypeError:
        print('Expected string input, got {} filename'.format(type(filename)))
    try:
        with open('{}'.format(filename),'r') as f1:
            lines = f1.readlines()
    except FileNotFoundError:
        print('Could not find file "{}", please choose a valid file name.'.format(filename))
    return lines
    
def set_incremental_resnum(filename):
    
    """
    
    Parameters
    ----------
    filename : string
        name of original .itp file to use as source
        
    Output
    ------
    None
        saves a new file to current directory: 
    """

    lines = get_lines(filename)

    num_lines = len(lines)

    old_lines = np.empty(num_lines, dtype="object")

    for i in range(num_lines):
        old_lines[i] = lines[i].split()

    new_lines = deepcopy(old_lines)

    residue_counter = 1

    for i in range(num_lines):

        if old_lines[i][0] == "ATOM":
            
            # Skip first ATOM entry
            if old_lines[i][1] == "1":
                continue

            # Every time a new group is found, increment the res_id value
            if int(old_lines[i][4]) != int(old_lines[i-1][4]):
                residue_counter += 1

                # THIS SECTION IS NOT GENERAL: USED TO ACCOUNT FOR INCORRECTLY GROUPED RESIDUES IN DRY_PEPTIDOGLYCAN.PDB
                if old_lines[i-1][2] == 'C4G':
                    new_lines[i-1][4] = residue_counter

            new_lines[i][4] = residue_counter      
            
    fileout = filename_split + "_uniqueres" + filetype

    savePDB(fileout, lines, new_lines)

    # with open(fileout, 'w+') as f:
    #     for line in new_lines:
    #         for element in line:
    #             f.write('{}\t'.format(element))
    #         f.write('\n')
            
    #     f.close()
    
    return None

def savePDB(fileout, old_lines, new_lines_split):
    """
    Function to format data into an output .pdb or .pqr file

    Requires that the new_lines_split argument contains elements for every column of every line for the new file

    Currently only the ATOM section is updated and so the old lines are passed to the function.

    to be used for all other sections of the file. This is done to streamline the formatting process, generalisation to be done later.

    Parameters:
    -----------

        fileout: string
            Name of output file

        old_lines: list of strings
            list of strings corresponding to the lines of the original .pdb file

        new_lines_split: list of lists of strings
            list containing the lists of updated column entries for each line in the new file
    """

    with open(fileout, 'w+') as f:
        for index, line in enumerate(new_lines_split):
            
            if line[0] != "ATOM":
                f.write(old_lines[index])

            elif line[0] == "ATOM":
                line[0] = line[0].ljust(6)#atom#6s
                # print("line[0]:", line[0], len(line[0]))
                line[1] = line[1].rjust(5)#aomnum#5d
                # print("line[1]:", line[1], len(line[1]))
                # print("line[2]:", line[2], len(line[2]))
                if len(line[2]) < 4:
                    line[2] = ' ' + line[2]
                line[2] = line[2].center(4)#atomname$#4s
                # print("line[2]:", line[2], len(line[2]))
                line[3] = line[3].ljust(6)#resname#1s
                # print("line[3]:", line[3], len(line[3]))
                # line[4] = line[4].ljust(3) #Astring : THIS IS THE STANDARD FORMAT BUT OUR GET_LINES LUMPS G4MP & X INTO G4MPX BECAUSE THERE IS NO WHITESPACE
                line[4] = str('%3d' % float(line[4])).rjust(3) #resnum
                # print("line[4]:", line[4], len(line[4]))
                line[5] = str('%8.3f' % float(line[5])).rjust(12) #x
                # print("line[5]:", line[5], len(line[5]))
                line[6] = str('%8.3f' % float(line[6])).rjust(8)#y
                line[7] = str('%8.3f' % float(line[7])).rjust(8) #z\
                line[8] =str('%6.2f'% float(line[8])).rjust(6)#occ
                line[9]=str('%6.2f'% float(line[9])).ljust(6)#temp
                # line[10]=line[10].rjust(12)#elname : THIS IS ALSO NOT INCLUDED IN OUR INPUT PDB   
                f.write("{0}{1} {2} {3}{4}{5}{6}{7}{8}{9}\n".format(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9]))
            

set_incremental_resnum(filename)

