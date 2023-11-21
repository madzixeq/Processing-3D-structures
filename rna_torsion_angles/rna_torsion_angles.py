import numpy as np
import pandas as pd
import argparse
from Bio.PDB import *
from calculate_angles import *

def getPDBFile(filepath):
    isFileCorrect = True
    try:
        parser = PDBParser()
        PDBFile = parser.get_structure("Struct", filepath)
    except FileNotFoundError:
        print("Error! File not found. Please make sure that filename and path to file are correct.")
        isFileCorrect = False  
        PDBFile = ""
    return isFileCorrect, PDBFile

def getTorsionAngleMatrix(PDBFile, includeSugarAngles):
    residues= [residue for residue in PDBFile.get_residues()]
    residuesWithoutHet = []
    for residue in residues:
        if(residue.get_resname() in ['A', 'C', 'U', 'G']):
            residuesWithoutHet.append(residue)

    angles = np.array([])
    anglesFromResidue = []
    for id, residue in enumerate(residuesWithoutHet):
        if id != 0:
            anglesFromResidue.append(alphaAngle(residuesWithoutHet[id-1],residue))
        else:
            anglesFromResidue.append(None)
        anglesFromResidue.append(betaAngle(residue))
        anglesFromResidue.append(gammaAngle(residue))
        anglesFromResidue.append(deltaAngle(residue))
        if id < (len(residuesWithoutHet) - 1):
            anglesFromResidue.append(epsilonAngle(residue, residuesWithoutHet[id+1]))
            anglesFromResidue.append(zetaAngle(residue, residuesWithoutHet[id+1]))
        else:
            anglesFromResidue.append(None)
            anglesFromResidue.append(None)
        anglesFromResidue.append(chiAngle(residue))
        if(includeSugarAngles):
            anglesFromResidue.append(v0Angle(residue))
            anglesFromResidue.append(v1Angle(residue))
            anglesFromResidue.append(v2Angle(residue))
            anglesFromResidue.append(v3Angle(residue))
            anglesFromResidue.append(v4Angle(residue))
        if id == 0:
            angles = np.array(anglesFromResidue)
        else:
            angles = np.vstack([angles,anglesFromResidue])
        anglesFromResidue=[]
    return angles

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input PDB file with RNA structure.", required=True)
    parser.add_argument("-o", "--output", help="CSV output file.", required=True)
    parser.add_argument("-s", "--includeSugar", help="A flag for including ribose torsion angles. If not provided, ribose torsion angles will not be included.", default=False, type=bool)
    args = parser.parse_args()
    isFileCorrect, PDBFile = getPDBFile(args.input)
    if isFileCorrect:
        angles = getTorsionAngleMatrix(PDBFile, args.includeSugar)
        DF = pd.DataFrame(angles)
        if not args.output.endswith(".csv"):
            args.output+=".csv"
        DF.to_csv(args.output, header=False, index=False)