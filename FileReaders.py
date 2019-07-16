"""
Author: Iain P. S. Smith
Date: 11:11 Monday 08/07/2019

Overview:
---------

    This script contains class objects that provide ease of access to the data contained within .pqr files 
    and a custom file format designed to contain the coarse-grained mapping schemes (atom groupings/partitioning 
    schemes, there are many nomenclatures).

    It is also designed to calculate the dipoles associated with each bead in a coarse-grained model as a result 
    of the differing electronegativities of the constituent atoms. The orientation is determined through 
    consideration of the electric field at the position of the parent bead:

        θ_i = arccos[ μ_i . E(r_i) / |μ_i||E(r_i)| ]                                                        (1)

    where:

        μ_i = q_ai(r_ai - r_i) + q_bi(r_bi - r_i) + q_ci(r_ci - r_i) + q_di(r_di - r_i)                     (2)

    Where μ_i is the dipole moment associated with the i-th bead and a, b, c, d denote the four daughter atoms that 
    comprise the i-th bead; with r_ai, q_ai being the coordinate point and charge associated with the relevant atom.

    We ignore the contribution of the static charge of the u-th bead to the local electric field to avoid division 
    by zero errors: also since the static charge contribution is radial it will have no effect on the orientation of 
    the dipole that is centred upon it.

    Once the magnitude and orientation are known, the dipoles are then placed at a distance +- 0.14 nm 
    from the parent bead (concordant with the MARTINI polarizable water model). The charge assigned to the dipole
    particles is then chose to fit the desired dipole magnitude according to the following equation:

        μ_i = q_i * r_i

    Furthermore, the contributions from the dipoles themselves are ignored when calculating the electric field
    to avoid the recursive problem of trying to define dipole orientations from the orientation of said dipoles.
    This is justified since the magnitude of the static charges is considerably larger than the magnitude of the
    dipoles and thus the relative effect of the dipoles is small.

    This script requires an atomistic structural data file (.pdb), a structural data file
    for the corresponding coarse-grained system and a file containing the grouped indices of the atoms
    that were mapped to each bead during the CG mapping process.

To be implemented:
------------------

    Need to adapt StructureData class to have two subclasses - SourceData & TargetData
    This would allow us to specify Source atomistic structures with methods for property generation
    and then a Target CG structure with methods for property manipulation.
    However in general one does not already have a CG .pqr file whilst building the system from atomistic data 
    and thus the use cases for this structure may be limited

Notes:
------
    Currently Hydrogen atoms are ignored in all of our methods, it may be desirable to consider these
    contributions, however other papers (see Drude Oscillator FF - https://pubs.acs.org/doi/full/10.1021/ct400781b)
    have deemed their results satisfactory without their inclusion.
"""

import numpy as np
from copy import deepcopy

class StructureData:
    """
    Defines a class object generated from .pqr file containing the indices, positions and charges of particles.

    Methods:

        get_lines(self): returns a list of all lines in .pqr file

        molecule_indices(self): returns the start and end bounds of the molecule data in terms of line indices

        get_positions(self): returns the position data for all particles within main molecule

        get_charges(self): returns the charge data for all particles within main molecule

        get_radii(self): returns the atomic radii data for all particles within main molecule
    
    Notes:
        The included methods assume that there is a single molecule whose data is located at the head of
        the .pqr file: i.e; solvent and ions are added below the molecule.
    
    """

    def __init__(self, filename):

        self.filename = filename

        self.lines = self.get_lines() # Pre-load data from file to reduce function calls within methods

        self.molecule_indexbounds = self.get_molecule_indices() # Pre-load molecule boundaries for the same reason

    def get_lines(self):
        """
        Performs a read lines method on a file of title "filename" and returns
        the result.
        
        Output
        ------
        lines : list of strings
            list containing strings of each line from file
        """
        try:
            type(self.filename) == str
        except TypeError:
            print('Expected string input, got {} filename'.format(type(self.filename)))
        try:
            with open('{}'.format(self.filename),'r') as f1:
                lines = f1.readlines()
        except FileNotFoundError:
            print('Could not find file "{}", please choose a valid file name.'.format(self.filename))
        return lines

    def get_molecule_indices(self):
        """
        Method to obtain the line number of the start and end of the molecule in the .pqr file. 
        Labelling starts at 0.

        Output
        ------
        bounds : list of ints - [start, end]
            line index for the start and end of data for the main molecule in pqr file
        """

        molecule_startline = 0 # First line of molecule structure data
        molecule_endline = 0   # Last line of molecule structure data

        found_start = False # Bool to track whether we know the molecule startline

        for line_index, line_string in enumerate(self.lines): # Iterate over all lines

            if line_string[0:4] == "ATOM": # Only consider the ATOM entries within the file

                # Set the startline of molecule to be at the first occurrence of ATOM
                if found_start == False:
                    molecule_startline = line_index
                    found_start = True
        
                line_split = line_string.split() # Obtain list of entry elements in current line

                # Set the endline of molecule to be the first occurrence of solvent data
                if line_split[3] == "W" or line_split[3] == "SOL":
                    molecule_endline = line_index
                    break
        
        return [molecule_startline, molecule_endline]

    def get_positions(self):
        """
        Method to obtain a list of positions for all constituent particles within the molecule
        """

        num_atoms = self.molecule_indexbounds[1] - self.molecule_indexbounds[0] # Calculate number of atoms

        positions = np.zeros((num_atoms, 3)) # Initialise numpy container for position data

        j = 0 # Dummy counter for assigning position data into relevant container index

        # Iterate over relevant data and assign the x, y, z coords into numpy container
        for i in range(self.molecule_indexbounds[0], self.molecule_indexbounds[1]):
            line_split = self.lines[i].split()
            positions[j,0] = line_split[5]
            positions[j,1] = line_split[6]
            positions[j,2] = line_split[7]
            j += 1
        
        return positions

    def get_charges(self):
        """
        Method to obtain a list of charges for all constituent particles within the molecule
        """

        num_atoms = self.molecule_indexbounds[1] - self.molecule_indexbounds[0] # Calculate number of atoms

        charges = np.zeros(num_atoms) # Initialise numpy container for charge data

        j = 0 # Dummy counter for assigning charge data into relevant container index

        # Iterate over relevant data and assign the scalar charge value into numpy container
        for i in range(self.molecule_indexbounds[0], self.molecule_indexbounds[1]):
            line_split = self.lines[i].split()
            charges[j] = line_split[8]
            j += 1
        
        return charges

    def get_radii(self):
        """
        Method to obtain a list of atomic radii for all constituent particles within the molecule
        """

        num_atoms = self.molecule_indexbounds[1] - self.molecule_indexbounds[0] # Calculate number of atoms

        radii = np.zeros(num_atoms) # Initialise numpy container for radii data

        j = 0 # Dummy counter for assigning atomic radii into relevant container index

        # Iterate over relevant data and assign the atomic radii into numpy container
        for i in range(self.molecule_indexbounds[0], self.molecule_indexbounds[1]):
            line_split = self.lines[i].split()
            radii[j] = line_split[9]
            j += 1
        
        return radii

    def get_atom_names(self):
        """
        Method to obtain an ordered list of atom names for use in labelling the output arrays from other methods.

        This method assumes that the maximum length of an atom name is 4 characters long, the last characters will
        be truncated if this assumption is not upheld.
        """

        num_atoms = self.molecule_indexbounds[1] - self.molecule_indexbounds[0] # Calculate number of atoms

        atom_names = np.empty(num_atoms, dtype="<U4") # Initialise numpy container for charge data

        j = 0 # Dummy counter for assigning charge data into relevant container index

        # Iterate over relevant data and assign the scalar charge value into numpy container
        for i in range(self.molecule_indexbounds[0], self.molecule_indexbounds[1]):
            line_split = self.lines[i].split()
            atom_names[j] = str(line_split[2])
            j += 1
        
        return atom_names

    def get_residue_numbers(self):
        """
        Method to obtain an ordered list of residue numbers for use in labelling the output arrays from other methods.
        """

        num_atoms = self.molecule_indexbounds[1] - self.molecule_indexbounds[0] # Calculate number of atoms

        residue_numbers = np.zeros(num_atoms, dtype="int") # Initialise numpy container for residue numbers

        j = 0 # Dummy counter for assigning charge data into relevant container index

        # Iterate over relevant data and assign the scalar charge value into numpy container
        for i in range(self.molecule_indexbounds[0], self.molecule_indexbounds[1]):
            line_split = self.lines[i].split()
            residue_numbers[j] = str(line_split[4])
            j += 1
        
        return residue_numbers

    def get_atom_indices(self):
        """
        Method to obtain an ordered list of atom indices for use in labelling the output arrays from other methods.

        This method assumes that the maximum length of an atom name is 4 characters long, the last characters will
        be truncated if this assumption is not upheld.
        """

        num_atoms = self.molecule_indexbounds[1] - self.molecule_indexbounds[0] # Calculate number of atoms

        atom_indices = np.empty(num_atoms, dtype="int") # Initialise numpy container for charge data

        j = 0 # Dummy counter for assigning charge data into relevant container index

        # Iterate over relevant data and assign the scalar charge value into numpy container
        for i in range(self.molecule_indexbounds[0], self.molecule_indexbounds[1]):
            line_split = self.lines[i].split()
            atom_indices[j] = int(line_split[1])
            j += 1
        
        return atom_indices

    def group_atoms_by_residue(self):
        """
        Method to generate an array of lists, one for each residue, containing the indices of all atoms within that residue
        """

        atom_indices = self.get_atom_indices()
        atom_residues = self.get_residue_numbers()

        # Find the total number of main molecule residues in the system
        num_residues = max(self.get_residue_numbers()) 
        # This method using max() works since the molecule by our assumptions starts at resnum=1

        num_atoms = self.molecule_indexbounds[1] - self.molecule_indexbounds[0] # Calculate number of atoms

        residue_atom_indices = np.empty(num_residues, dtype="object") # numpy container for atom indices grouped by residue number

        for i in range(num_residues):
            residue_atom_indices[i] = list() # Initialise the numpy container with a set of empty lists

        # The following code appends the atom index of each atom into the list contained within the residue_atom_indices
        # container whose index is the same as the atom's residue number (minus 1 since index counts from zero)
        for i in range(num_atoms):
            residue_atom_indices[atom_residues[i]-1].append(atom_indices[i])

        return residue_atom_indices

    def get_beads_by_index(self, PartitionScheme):

        """
        THIS METHOD SHOULD ONLY BE APPLIED TO .PQR FILES WITH RESIDUE NUMBERS THAT CORRESPOND TO MONOMER UNITS.
        Currently our coarse-grain mapping process generates structure files whose residue numbers follow the atomic indices.

        However, given that one should only ever be generating dipoles from an atomistic structure, this shouldn't be a problem.
        """

        named_atom_groups = PartitionScheme.get_atom_groups() # Atom groups for CG beads in terms of their atom names

        atom_names = self.get_atom_names()
        atom_indices = self.get_atom_indices()
        atom_residues = self.get_residue_numbers()

        central_atom_names = PartitionScheme.get_central_atoms()

        residue_atom_indices = self.group_atoms_by_residue()
        
        indexed_atom_groups = list()

        # Iterate through the data file entries to find the central atoms - preserving their order
        for atom_index in atom_indices:
            if atom_names[atom_index-1] in central_atom_names:
                indexed_atom_groups.append([atom_index])

        num_beads = len(indexed_atom_groups)

        # iterate over each central atom then
        # append the atoms that exist in the same bead as the central atom into the list containing the central atom index
        for bead_number in range(num_beads):

            central_atom_index = indexed_atom_groups[bead_number][0]
            central_atom_resnum = atom_residues[central_atom_index-1]
            central_atom_name = atom_names[central_atom_index-1]

            for group_index in range(len(named_atom_groups)): # iterate over all groups of atom names given in the PartitionScheme.txt file
                if central_atom_name in named_atom_groups[group_index]: # find which group the central atom belongs to
                    bead_atom_names = named_atom_groups[group_index] # save the list of atom names belonging to the central atom group

            for atom_index in residue_atom_indices[central_atom_resnum-1]: # iterate over all atom indices in the central atom residue
                atom_name = atom_names[atom_index-1] # save the current atom name as the atomline_index-1th element of the atom_names array
                if atom_name in bead_atom_names and atom_name != central_atom_name: # if current atom is in the central atom group and is not itself the central atom then do
                    indexed_atom_groups[bead_number].append(atom_index) # attach current atom index into central atom bead

        return indexed_atom_groups

    def get_beads_by_name(self, PartitionScheme):
        """
        Method to obtain the list of central atom names in the same order as the output from generate_dipoles.
        This will be used to check that the dipoles are being assigned to the correct atoms in the CG topology.
        """

        atom_names = self.get_atom_names() # Indexed list of atom names as they are ordered in the StructureData file
        bead_atom_indices = self.get_beads_by_index(PartitionScheme) # numpy array of length num_beads, each element is a list of atom indices contained within each bead

        num_beads = len(bead_atom_indices)

        bead_atom_names = np.empty(num_beads, dtype="object") # Create numpy container for bead atom names

        for i in range(num_beads):
            bead_atom_names[i] = list() # Initialise each container element as a list to be populated with the atom names in each bead

        for bead_id in range(num_beads):
            for atom_id in bead_atom_indices[bead_id]:
                    
                bead_atom_names[bead_id].append(atom_names[atom_id-1]) # StructureData file counts from 1, hence the -1 factor

        return bead_atom_names
    
    def get_beads_by_centralatoms(self, PartitionScheme):
        """
        Method to obtain the ordered list of beads labelled by their central atom names
        """

        bead_atom_names = self.get_beads_by_name(PartitionScheme)

        central_atoms = PartitionScheme.get_central_atoms()

        num_beads = len(bead_atom_names)

        bead_centralatom_names = np.empty(num_beads, dtype = "<U4")

        for bead_id in range(num_beads):

            for atom_name in bead_atom_names[bead_id]:
                if atom_name in central_atoms:
                    bead_centralatom_names[bead_id] = atom_name

        return bead_centralatom_names


    def generate_dipoles(self, PartitionScheme):
        """
        For use ONLY on the source data (e.g atomistic system)

        Method to generate the dipoles associated with each bead as a function of the constituent atom charges and positions
        The dipole is assumed to be centred on the central heavy atom.
        """

        central_atoms = PartitionScheme.get_central_atoms()

        bead_atom_indices = self.get_beads_by_index(PartitionScheme)
        atom_charges = self.get_charges()
        atom_positions = self.get_positions()
        atom_names = self.get_atom_names()

        num_beads = len(bead_atom_indices)

        dipoles = np.zeros((num_beads, 3), dtype="float") # Create vector of dipole moments for future population

        missing_central_atoms = list()
        
        for bead_number in range(num_beads):
            
            found_central_bead = False
            # Find index of central atom of the current bead
            for atom_index in bead_atom_indices[bead_number]:
                if atom_names[atom_index-1] in central_atoms: 
                    # Remember that the array indices count from 0, atom indices from data count from 1: hence the -1 above
                    central_atom_index = atom_index
                    found_central_bead = True
                
            if found_central_bead == False:
                # Skip calculation of dipole for any bead whose central atom is unknown
                continue

            cumulative_contribution = 0

            # Calculates the dipole moment resulting from the distribution of charge within the bead according to eq(2) in the docstring
            for atom_index in bead_atom_indices[bead_number]:

                # dx = atom_positions[atom_index-1, 0] - atom_positions[central_atom_index-1, 0]
                # dy = atom_positions[atom_index-1, 1] - atom_positions[central_atom_index-1, 1]
                # dz = atom_positions[atom_index-1, 2] - atom_positions[central_atom_index-1, 2]

                displacement = np.array(atom_positions[atom_index-1, :] - atom_positions[central_atom_index-1, :])

                # Ensure the particles are looking at their "nearest neighbour"
                # BEWARE: The necessity of this step may elude to an as of yet unforseen bug in the treatment of atom positions
                for i in range(len(displacement)):
                    if displacement[i] > 50: displacement[i] -= 100
                    elif displacement[i] < -50: displacement[i] +=100

                cumulative_contribution += atom_charges[atom_index-1] * displacement
                # Currently this does not account for any beads with missing constituent atoms, need to implement a predicted moment method

            dipoles[bead_number,:] = cumulative_contribution[:]

        return dipoles

    def generate_bead_charges(self, PartitionScheme):
        """
        Method to obtain the summed charges of constituent atoms in each bead.
        Values are returned in an ordered array similar to that returned by get_beads_by_index
        """

        bead_atom_indices = self.get_beads_by_index(PartitionScheme)
        atom_charges = self.get_charges()

        num_beads = len(bead_atom_indices)

        bead_charges = np.zeros(num_beads, dtype="float")

        for bead_id in range(num_beads):
            for atom_id in bead_atom_indices[bead_id]:
                bead_charges[bead_id] += atom_charges[atom_id-1]
        
        return bead_charges

    def show_data_for_group(self, PartitionScheme, group_id):
        """
        Method to print all charge and position data for the atoms in a given bead, indexed by group_id
        """

        bead_atom_indices = self.get_beads_by_index(PartitionScheme)
        atom_charges = self.get_charges()
        atom_positions = self.get_positions()
        atom_names = self.get_atom_names()

        if len(bead_atom_indices[group_id]) == 0:
            print("No atoms in group ", group_id)
            return 0

        else:
            print("Showing Data for bead group ", group_id)
            print("---------------------------------")

            for index in bead_atom_indices[group_id]:
                print("Atom name ", atom_names[index-1], ": Index ", index, ":\t Charge ", atom_charges[index-1], \
                      "\t Pos-xyz \t", atom_positions[index-1, 0], "\t", atom_positions[index-1, 1], "\t", atom_positions[index-1, 2])
            
            print("---------------------------------")
        
        return 0
        

class PartitionScheme:
    """
    Reads in a custom file format in the form:

        # Header - documentation
        # Header - documentation

        #GROUPDATA CentralAtom   Natoms  AtomNames
        GROUP data data data data data data
        GROUP data data data data data (data)
        GROUP data data data data data data
        GROUP data data data data data data (data)

    Note that the number of AtomNames per group may vary, this is why the number of atoms is also
    included in the data file so that one may easily iterate through the list of known length.

    """

    def __init__(self,filename):

        self.filename = filename

        self.lines = self.get_lines()

    def get_lines(self):
        """
        Performs a read lines method on a file of title "filename" and returns
        the result.
        
        Output
        ------
        lines : list of strings
            list containing strings of each line from file
        """
        try:
            type(self.filename) == str
        except TypeError:
            print('Expected string input, got {} filename'.format(type(self.filename)))
        try:
            with open('{}'.format(self.filename),'r') as f1:
                lines = f1.readlines()
        except FileNotFoundError:
            print('Could not find file "{}", please choose a valid file name.'.format(self.filename))
        return lines

    def data_linebounds(self):
        """
        Method to obtain the indices of the start and end lines of the partitioning group data
        """

        partitions_startline = 0 # First line of atom group partition data
        partitions_endline = 0   # Last line of atom group partition data

        found_start = False # Bool to track whether we know the startline of the grouping data

        for line_index, line_string in enumerate(self.lines): # Iterate over all lines

            if line_string[0:5] == "GROUP": # Only consider the GROUP entries within the file

                # Set the startline of molecule to be at the first occurrence of GROUP
                if found_start == False:
                    partitions_startline = line_index
                    found_start = True

                # Set the endline of molecule to be the last occurrence of GROUP data
                partitions_endline = line_index + 1 
                # The +1 allows for other methods to iterate up to and including the last entry
        
        return [partitions_startline, partitions_endline]

    def get_atom_groups(self):
        """
        Method to obtain the atom groups of the given partitioning scheme
        """

        partition_databounds = self.data_linebounds()

        num_groups = partition_databounds[1] - partition_databounds[0]

        # Initialise a numpy container for num_groups lists of varying size
        atom_groups = np.empty(num_groups, dtype="object")

        j = 0 # Dummy counter for assigning atom groups into relevant container index

        # Iterate over relevant data and assign the names of atoms in each group into numpy container
        for i in range(partition_databounds[0], partition_databounds[1]):

            current_group = list() # Empty list to contain an individual atom group
            line_split = self.lines[i].split()

            num_atoms = int(line_split[2]) # Read the number of atoms in current group from data file

            for name in range(3, 3 + num_atoms): # AtomName entries start at index 3 and run to 3 + num_atoms in data file
                current_group.append(line_split[name])

            atom_groups[j] = current_group
            j += 1

        return atom_groups

    def get_central_atoms(self):
        """
        Method to obtain the central atoms of each group in the given partitioning scheme
        """

        partition_databounds = self.data_linebounds()

        num_groups = partition_databounds[1] - partition_databounds[0]
        
        # Initialise a numpy container for num_groups strings of varying length
        central_atoms = np.empty(num_groups, dtype="<U4")

        j = 0 # Dummy counter for assigning atom groups into relevant container index

        # Iterate over relevant data and assign the names of atoms in each group into numpy container
        for i in range(partition_databounds[0], partition_databounds[1]):

            line_split = self.lines[i].split()

            central_atoms[j] = line_split[1] # CentralAtom name is index 1 on the relevant lines in data file
            j += 1

        return central_atoms


class TargetTopology:
    """
    Class object to contain methods for manipulation of topology properties

    This class is designed for use with the StructureData and PartitionScheme classes 
    to allow for the assignment of charges and dipoles into a Coarse Grained system.
    """

    def __init__(self, filename):

        self.filename = filename

        self.lines = self.get_lines() # Pre-load data from file to reduce function calls within methods

        self.atom_linebounds = self.get_atom_linebounds() # Pre-load molecule boundaries for the same reason

    def get_lines(self):
        """
        Performs a read lines method on a file of title "filename" and returns
        the result.
        
        Output
        ------
        lines : list of strings
            list containing strings of each line from file
        """
        try:
            type(self.filename) == str
        except TypeError:
            print('Expected string input, got {} filename'.format(type(self.filename)))
        try:
            with open('{}'.format(self.filename),'r') as f1:
                lines = f1.readlines()
        except FileNotFoundError:
            print('Could not find file "{}", please choose a valid file name.'.format(self.filename))
        return lines

    def get_atom_linebounds(self):

        """
        Method to obtain the line index of the start and end of the atom directive entries for the
        molecule (ignoring solvent and ions) in the topology file. 
        Labelling starts at 0.

        Output
        ------
        bounds : list of ints - [start, end]
            line index for the start and end of data for the main molecule atom data in topology
        """

        atom_startline = 0 # First line of atom data in topology
        atom_endline = 0   # Last line of atom data in topology

        found_start = False # Bool to track whether we know the atom startline - also used to "turn off" algorithm at the end of atom directive

        for line_index, line_string in enumerate(self.lines): # Iterate over all lines
            
            if "atoms" in line_string:
                atom_startline = line_index+1
                found_start = True

            elif found_start == True and line_string.strip() !='': # Iterate while in the atoms directive, until an empty line is found
                atom_endline = line_index # Keep updating atom_endline to the current line for as long as the above conditions hold
            
            else:
                found_start = False # Turn off algorithm once atom directive has been parsed and necessary information has been extracted.
        
        return [atom_startline, atom_endline]

    def set_molecule_atom_charges(self, StructureData, PartitionScheme):
        """
        Method to set the bead charges in the TargetTopology to be equal to the summed charges
        of the constituent atoms defined through the CG-mapping PartitionScheme and the relevant
        atomistic StructureData.
        """ 

        new_charges = StructureData.generate_bead_charges(PartitionScheme)

        lines = self.lines

        num_lines = len(lines)

        old_lines_split = np.empty(num_lines, dtype="object")

        for i in range(num_lines):
            old_lines_split[i] = lines[i].split()

        new_lines_split = deepcopy(old_lines_split)

        bead_index = 0

        for line_id in range(self.atom_linebounds[0], self.atom_linebounds[1]):
            
            new_lines_split[line_id][6] = new_charges[bead_index]

            bead_index += 1

        # THE FOLLOWING CODE FORMATS AND SAVES THE CHARGES TO A NEW TOPOLOGY FILE 
        #                  - THIS SHOULD BE UPDATE TO A GENERAL FUNCTION LATER ON
        # -----------------------------------------------------------------------

        filehandle = self.filename.strip('.top')
        fileout = filehandle + '_newcharges.top'

        with open(fileout, 'w+') as f:
            for line_id, line in enumerate(new_lines_split):
        
                if line_id in range(self.atom_linebounds[0], self.atom_linebounds[1]+1):
                    line[0] = line[0].ljust(8)#beadnum#8d
                    line[1] = line[1].ljust(8)#beadname#8s
                    line[2] = line[2].ljust(8)#resnum#8d
                    line[3] = line[3].ljust(8)#resname#8s
                    line[4] = line[4].ljust(8)#atomname#8s
                    line[5] = line[5].ljust(8)#chargegroup#8d
                    line[6] = str('%8.3f' % float(line[6])).ljust(8)#charge#8d

                    f.write("        {0}{1}{2}{3}{4}{5}{6}\n".format(line[0],line[1],line[2],line[3],line[4],line[5],line[6]))

                else:
                    f.write(lines[line_id])

        return None