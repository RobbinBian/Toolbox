import sys
sys.path.append("/media/bianyuan/data-TNG-1/Grab_TNGgalaxy/")
from Project_galaxy import deIllustrisTNG_galaxy
import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il
import csv
import os

run = 'TNG50'
basePath = '/media/bianyuan/data-TNG-1/' + run + '-1/output'

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
GroupMass = il.groupcat.loadHalos(basePath,99,fields=['GroupMass'])
Re = il.groupcat.loadSubhalos(basePath, 99, fields=['SubhaloHalfmassRadType'])[
    :, 4] * a * 1000 / h
mass = il.groupcat.loadSubhalos(
    basePath, 99, fields=['SubhaloMassInHalfRadType']) * 1e10 / h
star_mass = mass[:, 4]

desktop_path = "/media/bianyuan/data-TNG-1/data/"   # 新创建的txt文件的存放路径 
full_path = desktop_path + 'Groupmass.csv'          #也可以创建一个.doc的word文档

with open('/media/bianyuan/data-TNG-1/data/Groupmass.csv','a') as file:
    writer = csv.writer(file)
    writer.writerow(["ID","HostID","Groupmass"])
    file = open(full_path, 'w')                  # w 的含义为可进行读写                           
    file.close() 

for j in range(len(GroupFirstSub)):
    try:
        if GroupFirstSub[j] == -1:
            continue

        all_fields_of_cent = il.groupcat.loadSingle(basePath,99,subhaloID=GroupFirstSub[j])

        M = all_fields_of_cent['SubhaloMass'] * 1e10 / h
        rho_crit = (2.7754*1e11 * 1e-9) * (h**2)
        Rvir = ((3*M/(4*np.pi * 200*rho_crit)) ** (1/3))
        if Rvir == 0:
            continue
        x1 = all_fields_of_cent['SubhaloPos'][0]
        y1 = all_fields_of_cent['SubhaloPos'][1]
        z1 = all_fields_of_cent['SubhaloPos'][2]

        for i in range(GroupFirstSub[j],GroupFirstSub[j]+GroupNsubs[j]):
            if (10 < Re[i] < 600) & (1e8 < star_mass[i] < 1e10):
                with open('/media/bianyuan/data-TNG-1/data/Groupmass.csv','a') as file:
                    row = [i, GroupFirstSub[j], GroupMass[j] * 1e10 / h]
                    out = open("/media/bianyuan/data-TNG-1/data/Groupmass.csv", "a")
                    csv_writer = csv.writer(out, dialect = "excel")
                    csv_writer.writerow(row)
                    file.close() 
            else:
                continue
    except:
        continue