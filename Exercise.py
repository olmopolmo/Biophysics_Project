import subprocess
import argparse
import os
import Bio
import math
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import parse_pdb_header
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from math import e
###################################################################################################################
#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset

parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)

parser.add_argument(
    '--naccess',
    action='store',
    dest='naccess_bin',
    default=os.path.dirname(os.path.abspath(__file__)) + '/soft/NACCESS/naccess',
    help='Vdw parameters'
)
parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default=os.path.dirname(os.path.abspath(__file__)) + '/data/vdwprm',
    help='Vdw parameters'
)

parser.add_argument('pdb_file',help='Input PDB', type=open)
parser.add_argument('pdbqt_file',help='Input PDBQT', type=open)

args = parser.parse_args()

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Getting Charges and Atom type from PDBQT

params=[{}]

#Fix aton numbers when they do not start in 1
i = 1
for at in st.get_atoms():
    at.serial_number = i
    i += 1

for line in args.pdbqt_file:
    line = line.rstrip()
    #Skip TER records from PDBQT
    if line.find('TER') != -1:
        continue
    params.append({'charge': line[69:76], 'type': line[77:].replace(' ','')})

total_charge = 0.
for at in st.get_atoms():
    at.xtra['atom_type'] = params[at.serial_number]['type']
    at.xtra['charge'] = float(params[at.serial_number]['charge'])
    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
    total_charge += at.xtra['charge']


# Calculating surfaces
# Srf goes to .xtra['EXP_NACCESS'] field
srf = NACCESS_atomic(st[0], naccess_binary=args.naccess_bin)
###################################################################################################################

parser = PDBParser(PERMISSIVE=1)
def bash_command(c):
    subprocess.Popen(c,shell=True, executable="/bin/bash")
bash_command("cd /home/marc/Desktop/Biophysics_project")

#Parameters
MAXDIST = 3.4

#PDB list:

select = []

for at in st.get_atoms():
    select.append(at)


#calculation of E:

def elec(at1,at2,state):
    r=at1-at2
    if state=='vaccum':
        epsilon=1
        
    elif state=='water':
        epsilon=80
    else:
        epsilon=((86.9525)/(1-7.7839*e**(-0.3153*r)))-8.5525

    elec = (332.16*at1.xtra['charge']*at2.xtra['charge'])/(r*epsilon)
    return(elec)
def vdw(at1,at2,r):
    sig_i = vars(at1.xtra['vdw'])['sig']
    sig_j = vars(at2.xtra['vdw'])['sig']
    Ei = vars(at1.xtra['vdw'])['eps']
    Ej = vars(at2.xtra['vdw'])['eps']
    Ai = 2*((Ei)**0.5)*sig_i**6
    Aj = 2*((Ej)**0.5)*sig_j**6
    Ci = 2*((Ei)**0.5)*sig_i**3
    Cj = 2*((Ej)**0.5)*sig_j**3
    Evdw = ((Ai*Aj)/(r**12))-((Ci*Cj)/(r**6))  
    return(Evdw)  

def energy(atoms):
    counter=1
    Mutation=0
    Mvdw=0
    Melec=0
    Msolvation=0
    total_vdw=0
    c=0
    total_elec=0
    total_solvation=0
    i=1
    list_of_mutated_energies=[]
    for at1 in select:
        if at1.get_parent().get_id()[1]>counter:
            Mutation=Mvdw+Melec+Msolvation
            list_of_mutated_energies.append(Mutation)
            counter+=1
            Mutation=0
            Mvdw=0
            Melec=0
            Msolvation=0

        if at1.xtra['atom_type'] not in ('HD','H'):
            total_solvation= total_solvation + float(at1.xtra['EXP_NACCESS'])*vars(at1.xtra['vdw'])['fsrf']
        for at2 in select[i:]:
            r=at1-at2
            if vars(at2)['full_id'][3] != vars(at1)['full_id'][3]:
                #vdw parameters
                Evdw = vdw(at1,at2,r)
                #elec parameters  
                electrostatic = elec(at1,at2,'a')   
                if 2<r: 
                    total_vdw=total_vdw + Evdw
                if (str(at1)=='<Atom C>' and str(at2)== '<Atom N>') and abs(at1.get_parent().get_id()[1] - at2.get_parent().get_id()[1])==1:
                    c+=1
                else:
                    total_elec=total_elec + electrostatic
    
                #calculating the energy of the mutation
                if counter==at1.get_parent().get_id()[1]  and at1.get_parent().get_resname() not in ('GLY', 'ALA') and (str(at1) not in ('<Atom N>','<Atom C>','<Atom O>','<Atom CA>','<Atom CB>')) and 2<r:
                    Mvdw = Mvdw + vdw(at1,at2,r)
                    Melec = Melec + elec(at1,at2,'a')
                    if at.xtra['atom_type'] not in ('H','HD'):
                        Msolvation = Msolvation + float(at1.xtra['EXP_NACCESS'])*vars(at1.xtra['vdw'])['fsrf']
        i+=1
    print('Vdw: ',total_vdw)
    print('Elec: ',total_elec)
    print('Solvation: ', total_solvation)
    print('total E: ',total_elec+total_vdw+total_solvation)
    print()
    print('energy mutation: ', list_of_mutated_energies)

energy(select)

###########################################
#ASA
G_solvent=0
for at in st.get_atoms():
    if at.xtra['atom_type'] not in ('H','HD'):
        G_solvent= G_solvent + float(at.xtra['EXP_NACCESS'])*vars(at.xtra['vdw'])['fsrf']
print('∆G_solvent: ', G_solvent)

#writes a file with the G solvation of each aminoacid
f=open('data/AA_PDB/∆G_aminoacids.txt','a+')
f.write(str(G_solvent))


residuos={}
ids=[]
contador=0
for a in select:
    ids.append(vars(a)['full_id'][3])
    if (a.get_parent().get_resname()) not in residuos and ids[contador]!=ids[contador-1]:
        residuos[(a.get_parent().get_resname())] = 1
    elif ids[contador]!=ids[contador-1]:
        residuos[(a.get_parent().get_resname())] += 1
    contador+=1
print(residuos)

Amino = {'ALA': 2.174477, 'CYS': -5.705364000000002, 'ASP': 0.8736269999999999, 'GLU': -0.64498399999999961, 'PHE': 9.094826, 'GLY': -0.9312199999999997, 'HIS': 2.820512, 'ILE':-2.158062999999999, 'LYS': 0.6120529999999997, 'LEU': 4.758536, 'MET': 1.9830559999999995, 'ASN': 1.352348, 'PRO': -1.1933710000000002, 'GLN': -1.404772, 'ARG': -12.657179, 'SER': -2.0582459999999996, 'THR': -4.945886999999998, 'VAL': -3.2336639999999991, 'TRP': 5.5763230000000022, 'TYR': 0.276618}

unfolded=0

unfolded_mutations=[]

for i in residuos:
    if i not in ('HIE', 'NME','ACE'):
        unfolded= unfolded + residuos[i]*Amino[i]

for i in residuos:
    if i not in ('HIE','GLY', 'NME','ACE'):
        unfolded_mutations.append(unfolded - Amino[i] + Amino['ALA'])
    
print('unfolded E:', unfolded)
print('undfolded_mutations:', unfolded_mutations)

