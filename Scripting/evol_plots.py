from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Define a base path
base_path=Path.cwd()
thesis_chapter_path_vector = Path.cwd().parent.parent / 'Thesis' / 'phd-thesis-template-2.4' / 'Chapter3' / 'Figs' / 'Vector'
thesis_chapter_path_raster = Path.cwd().parent.parent / 'Thesis' / 'phd-thesis-template-2.4' / 'Chapter3' / 'Figs' / 'Raster'
evol_5mm_path=base_path / 'evol_5mm.csv'
evol_3mm_path=base_path / 'evol_3mm.csv'
#evol_2mm_path=base_path / 'evol_2mm.csv'
evol_10mm_path=base_path / 'evol_10mm.csv'

# Read the data
evol_5mm=pd.read_csv(evol_5mm_path)
evol_3mm=pd.read_csv(evol_3mm_path)
#evol_2mm=pd.read_csv(evol_2mm_path)
evol_10mm=pd.read_csv(evol_10mm_path)

# fake data for testing

#evol_3mm=evol_5mm.copy()
#evol_3mm['EVOL Sum']=0.5*evol_5mm['EVOL Sum']

evol_2mm=evol_5mm.copy()
evol_2mm['EVOL Sum']=1.05*evol_5mm['EVOL Sum']

#evol_10mm=evol_5mm.copy()
#evol_10mm['EVOL Sum']=1.1*evol_5mm['EVOL Sum']

#print(evol_5mm.head())

#ys.exit()

# Find peak EVOL 
peak_5mm_evol=evol_5mm['EVOL Sum'].max()
peak_3mm_evol=evol_3mm['EVOL Sum'].max()
peak_2mm_evol=evol_2mm['EVOL Sum'].max()
peak_10mm_evol=evol_10mm['EVOL Sum'].max()
"""
#fake data for testing
peak_3mm_evol=0.5*peak_5mm_evol
peak_2mm_evol=0.25*peak_5mm_evol
peak_10mm_evol=1.1*peak_5mm_evol
"""
# Mesh sizes
mesh_sizes=[2, 3, 5, 10]

# plotting peak EVOL vs mesh size
plt.figure(figsize=(10,6))
plt.scatter(mesh_sizes, [peak_2mm_evol, peak_3mm_evol, peak_5mm_evol, peak_10mm_evol], color='b', label='Peak EVOL')
plt.xlabel('Mesh size (mm)')
plt.ylabel('Peak Cerebrum Volume (mm³)')
plt.xlim(0)
#plt.ylim(0)
plt.savefig(thesis_chapter_path_vector / 'peak_cerebrum_vol_vs_mesh_size.png')
plt.savefig(thesis_chapter_path_raster / 'peak_cerebrum_vol_vs_mesh_size.png')
plt.show()


# plotting EVOL vs time for each mesh
plt.figure(figsize=(10,6))
plt.plot(evol_2mm['Frame'], evol_2mm['EVOL Sum'], color='r', label='2mm')
plt.plot(evol_3mm['Frame'], evol_3mm['EVOL Sum'], color='g', label='3mm')
plt.plot(evol_5mm['Frame'], evol_5mm['EVOL Sum'], color='b', label='5mm')
plt.plot(evol_10mm['Frame'], evol_10mm['EVOL Sum'], color='y', label='10mm')
plt.xlabel('Frame')
plt.ylabel('Cerebrum Volume (mm³)')
plt.xlim(1)
#plt.ylim(0)
plt.legend()
plt.savefig(thesis_chapter_path_vector / 'cerebrum_vol_vs_time.png')
plt.savefig(thesis_chapter_path_raster / 'cerebrum_vol_vs_time.png')
plt.show()




