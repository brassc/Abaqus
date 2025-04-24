from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import sys

def set_publication_style():
	"""Set matplotlib parameters for publication-quality figures."""
	plt.rcParams.update({
		'font.family': 'serif',
		'font.serif': ['Times New Roman'],
		'mathtext.fontset': 'stix',
		'font.size': 12,
		'axes.labelsize': 14,
		'axes.titlesize': 16,
		'axes.titleweight': 'bold', # This makes titles bold
		'xtick.labelsize': 12,
		'ytick.labelsize': 12,
		'legend.fontsize': 10,
		'figure.dpi': 150,
		'savefig.dpi': 300,
		'savefig.format': 'png',
		'savefig.bbox': 'tight',
		'axes.grid': True,
		'grid.alpha': 0.3,
		'grid.linestyle': '-',
		'axes.spines.top': False,
		'axes.spines.right': False,
	})

# Define a base path
base_path=Path.cwd()
print(base_path)

thesis_chapter_path_vector = base_path.parent / 'Thesis' / 'phd-thesis-template-2.4' / 'Chapter3' / 'Figs' / 'Vector'
thesis_chapter_path_raster = base_path.parent / 'Thesis' / 'phd-thesis-template-2.4' / 'Chapter3' / 'Figs' / 'Raster'
evol_10mm_path=base_path / 'Scripting/evol_10mm.csv'
evol_7mm_path=base_path / 'Scripting/evol_7mm.csv'
evol_5mm_path=base_path / 'Scripting/evol_5mm.csv'
evol_4mm_path=base_path / 'Scripting/evol_4mm.csv'
evol_3mm_path=base_path / 'Scripting/evol_3mm.csv'




# Read the data
evol_10mm=pd.read_csv(evol_10mm_path)
evol_7mm=pd.read_csv(evol_7mm_path)
evol_5mm=pd.read_csv(evol_5mm_path)
evol_4mm=pd.read_csv(evol_4mm_path)
evol_3mm=pd.read_csv(evol_3mm_path)
#evol_2mm=pd.read_csv(evol_2mm_path)


# fake data for testing

#evol_3mm=evol_5mm.copy()
#evol_3mm['EVOL Sum']=0.5*evol_5mm['EVOL Sum']

#evol_2mm=evol_5mm.copy()
#evol_2mm['EVOL Sum']=1.05*evol_5mm['EVOL Sum']

#evol_10mm=evol_5mm.copy()
#evol_10mm['EVOL Sum']=1.1*evol_5mm['EVOL Sum']

#print(evol_5mm.head())

#ys.exit()

# Find evol at frame 25 (index 26)
index=21 # for frame 20
evol_10mm_frame25=evol_10mm.loc[index, 'EVOL Sum']
evol_7mm_frame25=evol_7mm.loc[index, 'EVOL Sum']
evol_5mm_frame25=evol_5mm.loc[index, 'EVOL Sum']
evol_4mm_frame25=evol_4mm.loc[index, 'EVOL Sum']
evol_3mm_frame25=evol_3mm.loc[index, 'EVOL Sum']


#peak EVOL
peak_10mm_evol=evol_10mm['EVOL Sum'].max()
peak_7mm_evol=evol_7mm['EVOL Sum'].max()
peak_5mm_evol=evol_5mm['EVOL Sum'].max()
peak_4mm_evol=evol_4mm['EVOL Sum'].max()
peak_3mm_evol=evol_3mm['EVOL Sum'].max()
#peak_2mm_evol=evol_2mm['EVOL Sum'].max()


# Mesh sizes
mesh_sizes=[3, 4, 5, 7, 10]

# plotting peak EVOL vs mesh size
plt.figure(figsize=(10,6))
#plt.scatter(mesh_sizes, [peak_3mm_evol, peak_5mm_evol, peak_10mm_evol], color='b', label='Peak EVOL')
plt.scatter(mesh_sizes, [evol_3mm_frame25, evol_4mm_frame25, evol_5mm_frame25, evol_7mm_frame25, evol_10mm_frame25], color='b', label='EVOL at Frame 25', s=15)

plt.xlabel('Mesh Size (mm)')
plt.ylabel('Cerebrum Volume (mm³)')
plt.xlim(0, 12)
plt.ylim(0, 1300000)
plt.title(f'Cerebrum Volume (mm³) vs. Mesh Size in Step 2, Frame {index}')
plt.savefig(thesis_chapter_path_vector / 'peak_cerebrum_vol_vs_mesh_size_more.png')
plt.savefig(thesis_chapter_path_raster / 'peak_cerebrum_vol_vs_mesh_size_more.png')
plt.show()


# plotting EVOL vs time for each mesh
plt.figure(figsize=(10,6))
#plt.plot(evol_2mm['Frame'], evol_2mm['EVOL Sum'], color='r', label='2mm')
plt.plot(evol_3mm['Frame'], evol_3mm['EVOL Sum'], color='g', label='3mm')
plt.plot(evol_4mm['Frame'], evol_4mm['EVOL Sum'], color='m', label='4mm')
plt.plot(evol_5mm['Frame'], evol_5mm['EVOL Sum'], color='b', label='5mm')
plt.plot(evol_7mm['Frame'], evol_7mm['EVOL Sum'], color='c', label='7mm')
plt.plot(evol_10mm['Frame'], evol_10mm['EVOL Sum'], color='y', label='10mm')
plt.xlabel('Frame')
plt.ylabel('Cerebrum Volume (mm³)')
plt.xlim(1)
#plt.ylim(0)
plt.title('Cerebrum Volume (mm³) over Time')
plt.legend()
plt.savefig(thesis_chapter_path_vector / 'cerebrum_vol_vs_time_more.png')
plt.savefig(thesis_chapter_path_raster / 'cerebrum_vol_vs_time_more.png')
plt.show()




