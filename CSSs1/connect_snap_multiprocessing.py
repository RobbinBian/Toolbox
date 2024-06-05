import pandas as pd
from tqdm import tqdm
import h5py
import csv
import sys
sys.path.append("/home/bianyuan/workspace/Grab_TNGgalaxy/")
import numpy as np
import pandas as pd
import h5py
import csv
import os
from multiprocessing import Pool
import illustris_python as il
from tqdm import tqdm
# %matplotlib widget
basePath = '/media/bianyuan/data-TNG-1/TNG50-1/output'
# from astroML.plotting import setup_text_plots
#Lets text in plots use latex
# setup_text_plots(usetex=True)

def process_snap(snap):
    full_path = desktop_path + 'MRI_output_'+str(snap)+'.csv'

    # 确保目录存在，如果不存在则创建
    os.makedirs(os.path.dirname(full_path), exist_ok=True)
    
    with open(full_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ID', 'Snap', 'Galaxy_ID', 'Number_of_Stars'])

    with h5py.File(f'/media/bianyuan/data-TNG-3/StarID/starID_1e7_gals_snap_{snap}.hdf5', 'r') as f:
        for ids in IDs:
            star_ids = il.snapshot.loadSubhalo(basePath, 99, ids, 'star', fields='ParticleIDs')
            for galaxy_id, group in f.items():
                galaxy_star_ids = group['StarID'][:]
                common_stars = set(star_ids).intersection(galaxy_star_ids)

                if common_stars:
                    with open(full_path, 'a', newline='') as file:
                        writer = csv.writer(file)
                        row = [ids, snap, galaxy_id, len(common_stars)]
                        writer.writerow(row)
                        print(f"Found {len(common_stars)} matching stars in galaxy ID {galaxy_id} for CSSs3 ID {ids} in snap {snap}.")

# specified_snaps = [68, 67, 65, 63, 62, 61, 60, 59, 59, 58, 57, 56]
# specified_snaps = [98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40]
specified_snaps = [68, 67]

CSSs3 = pd.read_csv('/home/bianyuan/workspace/data/CSSs3.csv')
# IDs = np.array(CSSs3[(CSSs3['Traced or Not'] == 0)]['ID'])
IDs = np.array(CSSs3['ID'])[np.where( ((13.8-CSSs3['MWted-Age']) >= 2.9) & ((13.8-CSSs3['MWted-Age']) <= 6.4) & (np.array(CSSs3['[Z/H]_r']) > 0.1) & (CSSs3['Number'] > 1) & (CSSs3['Distance_grpvir'] > 0) ) ]
desktop_path = "/home/bianyuan/workspace/data/connect_snap/" 

# if __name__ == "__main__":
#     with Pool() as pool:
#         pool.map(process_snap, specified_snaps)
        
if __name__ == "__main__":
    # 限制进程池的大小以减少内存使用
    with Pool(processes=8) as pool:  # 减少进程数量以降低内存占用
        pool.map(process_snap, specified_snaps)