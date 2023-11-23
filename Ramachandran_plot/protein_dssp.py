from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


path = "/home/magdazyd/torsion_angles/1prr.pdb"
p = PDBParser()
structure=p.get_structure("1MBO",path)
print(structure)

for item in structure.get_models():
    print(item)
model=structure[0]
dssp = DSSP(model,path,dssp='/usr/bin/dssp')

# (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
# NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
# NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)


kolor_mapping = {'H': 'red','B':'orange','E':'yellow', 'G':'blue','I':'violet','T': 'pink','S':'green', '-': 'black'}
for amino in list(dssp.keys()):
   
    plt.scatter(dssp[amino][4], dssp[amino][5], c=kolor_mapping[dssp[amino][2]],alpha=0.7)
    #print(dssp[amino])
    # structures=(dssp[amino][2])
    # phi=(dssp[amino][4])
    # psi=(dssp[amino][5])

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Alpha helix'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=10, label='Isolated beta-bridge'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='yellow', markersize=10, label='Strand'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='3-10 helix'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='violet', markersize=10, label='Pi helix'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='pink', markersize=10, label='Turn'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='Bend'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='None')
]

plt.legend(handles=legend_elements)
plt.title('Ramachandran plot')
plt.xlabel('Phi value')
plt.ylabel('Psi value')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.xticks([-180,-120,-60,0,60,120, 180]) 
plt.yticks([-180,-120,-60,0,60,120, 180]) 
plt.savefig('my_plot.png')