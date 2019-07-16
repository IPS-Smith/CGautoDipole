#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 11:31:57 2019

@author: iain
"""

import sys

arglist = sys.argv

filename = arglist[1]

filename_split = filename.strip(".pqr")

directory = arglist[2]

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


def formatPQR(filename):
    
    lines = get_lines(filename)
    
    lines_split = list()

    for i in range(len(lines)):
        lines_split.append(lines[i].split())
        
    fileout = filename_split + "_formatted.pqr"
    
    fileout = directory + '/' + fileout
     
    with open(fileout, 'w+') as f:
        
        for index, line in enumerate(lines):
            if(lines_split[index][0] == 'ATOM'):
                f.write("%4s %6s %4s %4s %3s %11s %7s %7s %5s %5s\n" % (lines_split[index][0], lines_split[index][1], 
                                                                      lines_split[index][2], lines_split[index][3],
                                                                      lines_split[index][4], lines_split[index][5], 
                                                                      lines_split[index][6], lines_split[index][7], 
                                                                      lines_split[index][8], lines_split[index][9]))
            else:
                f.write(line)
            
        f.close()

    print("Formatted file saved to: ", fileout)
    
    return 0

formatPQR(filename)
    
    
