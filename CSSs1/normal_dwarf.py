import sys
sys.path.append("/home/bianyuan/workspace/Grab_TNGgalaxy/")
from Project_galaxy import deIllustrisTNG_galaxy
# from mpl_toolkits import mplot3d
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il
# import pandas as pd
import csv
import os
# import seaborn as sns

run = 'TNG50'
basePath = '/media/bianyuan/data-TNG-1/' + run + '-1/output'

G = 4.3009*1e-6
a = 1
H0 = 67.74
h = H0/100.
WM = 0.3089
WV = 0.6911
Tyr = 977.8
WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
WK = 1-WM-WR-WV
n = 1000

GroupFirstSub = il.groupcat.loadHalos(basePath, 99, fields=['GroupFirstSub'])
GroupNsubs = il.groupcat.loadHalos(basePath, 99, fields=['GroupNsubs'])
Group_M_Crit200 = il.groupcat.loadHalos(basePath, 99, fields=['Group_M_Crit200'])
Re = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloHalfmassRadType'])[
    :, 4] * a * 1000 / h
mass = il.groupcat.loadSubhalos(
    basePath, 99, fields=['SubhaloMassInHalfRadType']) * 1e10 / h
star_mass = mass[:, 4]

# metal = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetallicity'])
C = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,2]
N = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,3]
O = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,4]
Ne = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,5]
Mg = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,6]
Si = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,7]
Fe = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,8]
total = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,9]
H = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloStarMetalFractions'])[:,0]

Pos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos'])


desktop_path = "/home/bianyuan/workspace/data/"   # 新创建的txt文件的存放路径 
full_path = desktop_path + 'Normal_Dwarf_2.csv'          #也可以创建一个.doc的word文档

with open('/home/bianyuan/workspace/data/Normal_Dwarf_2.csv','w') as file:
    writer = csv.writer(file)
    writer.writerow(["ID","HostID","[Z/H]","MWted-Age","Star mass","Stellarmass","Re","Distance","Traced or Not"])
    file = open(full_path, 'w')                  # w 的含义为可进行读写                           
    file.close() 

# cggal_id, fgal_id = [], []
age_c, age_g, age_f = [], [], []
Z_c, Z_g, Z_f = [], [], []

dis = np.zeros(len(Pos))
dis_norm = np.zeros(len(Pos))
# for j in range(2):
for j in range(len(GroupFirstSub)):
    # print(j)
    try:
        if GroupFirstSub[j] == -1:
            continue

        all_fields_of_cent = il.groupcat.loadSingle(basePath,99,subhaloID=GroupFirstSub[j])

        # M = all_fields_of_cent['SubhaloMass'] * 1e10 / h
        # rho_crit = (2.7754*1e11 * 1e-9) * (h**2)
        # Rvir = ((3*M/(4*np.pi * 200*rho_crit)) ** (1/3))

        M = Group_M_Crit200[j] * 1e10 / h
        # rho_crit = (2.7754*1e11 * 1e-9) * (h**2)
        rho_crit0 = 3 * (H0*1e-3)**2 / (8*np.pi*G)
        Rvir = ((3*M/(4*np.pi * 200*rho_crit0)) ** (1/3))
        if Rvir == 0:
            continue
        x1 = all_fields_of_cent['SubhaloPos'][0]
        y1 = all_fields_of_cent['SubhaloPos'][1]
        z1 = all_fields_of_cent['SubhaloPos'][2]

        for i in range(GroupFirstSub[j],GroupFirstSub[j]+GroupNsubs[j]):
            if ( Re[i] > 600) & (1e8 < star_mass[i] < 1e10):
                galaxy = deIllustrisTNG_galaxy(i, deproject=True, align_with='star', radius_align_max=10.,
                                            radius_align_min=0., centmode='pot', snapshot=99, basePath=basePath, run=run)
                az = galaxy.s['tform']  # scale factor
                m = galaxy.s['mass']
                age = 0.
                for jj in range(n):
                    a = az*(j+0.5)/n
                    adot = (WK+(WM/a)+(WR/(a*a))+(WV*a*a))**0.5
                    age = age + adot**(-1)

                zage = az*age/n
                zage_Gyr = (Tyr/H0)*zage  # this is calculate age
                # zage_LBT = 13.8 - zage_Gyr
                MWage = np.dot(m,zage_Gyr)/np.sum(m)

                all_fields_of_sat = il.groupcat.loadSingle(basePath,99,subhaloID=i)
                mass = all_fields_of_sat['SubhaloMassInHalfRadType'][4] * 1e10 / h
                x2= all_fields_of_sat['SubhaloPos'][0]
                y2= all_fields_of_sat['SubhaloPos'][1]
                z2= all_fields_of_sat['SubhaloPos'][2]
                dis = (np.sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))/h) / Rvir
                print(dis)
                RowNum, _, _ = il.sublink.treeOffsets(basePath, 99, i, "SubLink")
                print(np.int(RowNum>0))

                with open('/home/bianyuan/workspace/data/Normal_Dwarf_2.csv','a') as file:
                    row = [i, GroupFirstSub[j], np.log10((np.sum(C[i]+N[i]+O[i]+Ne[i]+Mg[i]+Si[i]+Fe[i]+total[i])/H[i])/(0.0127/0.7381)), MWage, star_mass[i], np.sum(m), Re[i], dis, np.int(RowNum>0)]
                    out = open("/home/bianyuan/workspace/data/Normal_Dwarf_2.csv", "a")
                    csv_writer = csv.writer(out, dialect = "excel")
                    csv_writer.writerow(row)
                    file.close() 


                # if GroupNsubs[j] > 50:
                #     # cgal_id.append(i)
                #     age_c.append(MWage)
                #     Z_c.append(np.log10((np.sum(C[i]+N[i]+O[i]+Ne[i]+Mg[i]+Si[i]+Fe[i]+total[i])/H[i])/(0.0196/0.7)))
                #     with open('/media/bianyuan/data-TNG-1/data/CSSs.csv','a') as file:
                #         row = [i, np.log10((np.sum(C[i]+N[i]+O[i]+Ne[i]+Mg[i]+Si[i]+Fe[i]+total[i])/H[i])/(0.0196/0.7)), MWage, star_mass[i], Re[i], dis, np.int(RowNum>0), 'cluster galaxy']
                #         out = open("/media/bianyuan/data-TNG-1/data/CSSs.csv", "a")
                #         csv_writer = csv.writer(out, dialect = "excel")
                #         csv_writer.writerow(row)
                #         file.close() 
                #     print(i)
                # elif 0 < GroupNsubs[j] < 50:
                #     age_g.append(MWage)
                #     Z_g.append(np.log10((np.sum(C[i]+N[i]+O[i]+Ne[i]+Mg[i]+Si[i]+Fe[i]+total[i])/H[i])/(0.0196/0.7)))
                #     with open('/media/bianyuan/data-TNG-1/data/CSSs.csv','a') as file:
                #         row = [i, np.log10((np.sum(C[i]+N[i]+O[i]+Ne[i]+Mg[i]+Si[i]+Fe[i]+total[i])/H[i])/(0.0196/0.7)), MWage, star_mass[i], Re[i], dis, np.int(RowNum>0), 'group galaxy']
                #         out = open("/media/bianyuan/data-TNG-1/data/CSSs.csv", "a")
                #         csv_writer = csv.writer(out, dialect = "excel")
                #         csv_writer.writerow(row)
                #         file.close() 
                # else:
                #     # fgal_id.append(i)
                #     age_f.append(MWage)
                #     Z_f.append(np.log10((np.sum(C[i]+N[i]+O[i]+Ne[i]+Mg[i]+Si[i]+Fe[i]+total[i])/H[i])/(0.0196/0.7)))
                #     with open('/media/bianyuan/data-TNG-1/data/CSSs.csv','a') as file:
                #         row = [i, np.log10((np.sum(C[i]+N[i]+O[i]+Ne[i]+Mg[i]+Si[i]+Fe[i]+total[i])/H[i])/(0.0196/0.7)), MWage, star_mass[i], Re[i], dis, np.int(RowNum>0), 'field galaxy']
                #         out = open("/media/bianyuan/data-TNG-1/data/CSSs.csv", "a")
                #         csv_writer = csv.writer(out, dialect = "excel")
                #         csv_writer.writerow(row)
                #         file.close()
                #     print(i)
            else:
                continue
    except:
        continue

# plt.figure(figsize=(12, 8))
# # plt.scatter(age_c, Z_c, label='cluster CSSs')
# # plt.scatter(age_g, Z_g, label='group CSSs')
# # plt.scatter(age_f, Z_f, label='field CSSs')
# # sns.jointplot(age_cg, Z_cg, label='cluster & group CSSs')
# # sns.jointplot(age_f, Z_f, label='field CSSs')

# plt.xlabel('Age[Gyr]')
# plt.ylabel('[Z/H]')
# plt.savefig('/media/bianyuan/data-TNG-1/pic/z_age'+'.jpg')