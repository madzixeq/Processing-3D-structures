import math
import numpy as np
import matplotlib.pyplot as plt
import argparse
from Bio.PDB import *

def contacts_map(PDBFile, dist, chain_id, model_id, output):
    parser = PDBParser()
    try:
        structure = parser.get_structure("structure", PDBFile)
        if model_id == None:
            for i in structure.get_models():
                model = i
                break
        else:
            model = structure[model_id]
        if chain_id == None:
            for i in model:
                chain = i
                break
        else:
            chain = model[chain_id]
        list_of_amin=[]
        for residue in chain:
            if "CA" in residue:
                list_of_amin.append(residue["CA"])
        matrix = np.zeros((len(list_of_amin),len(list_of_amin)))
        for i, x in enumerate(list_of_amin):
            for j, y in enumerate(list_of_amin):
                if (x-y)<=dist:
                    matrix[i][j] = 1
        plot_map(matrix, output)
    except FileNotFoundError:
        print("Error! File not found. Please make sure that filename and path to file are correct.")
    except KeyError:
        print("Error! Model number or chain code are incorrect.")

def plot_map(matrix, output):
    xpoints=[]
    ypoints=[]
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j]==1:
                xpoints.append(i)
                ypoints.append(j)
    plt.scatter(xpoints,ypoints,s=1)
    axes=plt.gca()
    axes.get_xaxis().set_visible(False)
    axes.get_yaxis().set_visible(False)
    if output:
        if not output.endswith(".jpg"):
            output += ".jpg"
        plt.savefig(output)
    plt.show()


parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="Input file in PDB format", required=True)
parser.add_argument("-d","--distance",help="Distance threshold in Angstrom. Default value: 8", type = int, default=8)
parser.add_argument("-c","--chain",help="Chain code from PDB file. If not provided, program takes first chain from model.")
parser.add_argument("-m","--model",help="Model number from PDB file. If not provided, program takes first model from file.")
parser.add_argument("-o","--output",help="Output .jpg file in which the resulting plot will be saved. If not provided, plot will not be saved to your computer!")
args = parser.parse_args()

contacts_map(args.input, args.distance, args.chain, args.model, args.output)


