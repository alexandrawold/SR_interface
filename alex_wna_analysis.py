
import numpy as np
import datetime as dt
import os
from constants_settings import *
from convert_coords import convert2
from satellites import sat
from bfield import getBdir
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile, read_damp_simple
from ray_plots import plotray2D, plotrefractivesurface, plot_plasmasphere_2D
import PyGeopack as gp
from spacepy import coordinates as coord
import matplotlib.pyplot as plt



C = 3e8 #m/s
Re = 6371e3

rayfile_directory = '/home/alwo2026/Workspace/SR_output'

wnas_all = []
freq_all = []
coor_all = []

for lats in [-60, -40, -20, 0, 20, 40, 60]:
    for lons in [-180, -120, -60, 0, 60, 120, 180]:
        pos_array = [[500.0e3+R_E, lats, lons]]
        dt_array = [dt.datetime(2014,1,1,12,0,0, tzinfo=dt.timezone.utc)]
        crs1 = 'GEO'
        crs2 = 'SM'
        carsph1 = 'sph'
        carsph2 = 'car'
        units1 = ['m','deg','deg']
        units2 = ['m','m','m']
        ray_start = convert2(pos_array, dt_array, crs1, carsph1, units1, crs2, carsph2, units2)

        # setting frequeinces
        freqs = [1.0e3, 2.0e3, 3.0e3, 4.0e3, 5.0e3, 6.0e3, 7.0e3, 8.0e3] # Hz

        # Which plasmasphere model should we run?
        #   1 - Ngo model - do NOT use
        #   6 - Simplified GCPM from Austin Sousa's thesis
        #   7 - New! Diffusive Equilibrium AT64ThCh (see docs)
        md = 7

        ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(dt_array[0], '%Y-%m-%d_%H_%M_%S')

        file_titles = os.listdir(ray_out_dir)

        raylist = []
        mode_name = 'mode' + str(md)
        fn_str = 'lat%s_lon%s'%(pos_array[0][1], pos_array[0][2])

        # use mode name to avoid workers of the same label
        x = sorted(file_titles)

        # print(len(x))

        for filename in x:
            if '.ray' in filename and mode_name in filename and fn_str in filename:
                raylist += read_rayfile(os.path.join(ray_out_dir, filename))
                # print(filename)

        # print(raylist[0]['w'])

         # Calculate WNAs
        
        for r in range(0, len(raylist)):
            ray = raylist[r]
            w = ray['w']
            tmp_kcoords = list(zip((w/C) * ray['n'].x, (w/C) * ray['n'].y, (w/C) * ray['n'].z))
            kcoords = [(tmp_kcoords[s][0], tmp_kcoords[s][1], tmp_kcoords[s][2]) / np.sqrt(tmp_kcoords[s][0]**2 + tmp_kcoords[s][1]**2 + tmp_kcoords[s][2]**2) for s in range(len(tmp_kcoords))]
            coords = list(zip(ray['pos'].x/Re, ray['pos'].y/Re, ray['pos'].z/Re))
            car_cvals = coord.Coords(coords, 'SM', 'car')
            sph_cvals = car_cvals.convert('SM', 'sph')
            for ti in range(0, len(ray['time'])):
                if np.abs(sph_cvals.lati[ti])<=10:
                    B   =  ray['B0'].iloc[ti]
                    Bmag = np.linalg.norm(B)
                    bunit = B/Bmag
                    kunit = np.array([kcoords[ti][0], kcoords[ti][1], kcoords[ti][2]])
                    alpha = np.arccos(np.dot(kunit, bunit))
                    wnas_all.append(alpha)
                    freq_all.append(w)
                    coor_all.append(sph_cvals[ti])

# print(wnas_all)

wnas = np.asarray(wnas_all)
print(np.shape(wnas))
freq = np.asarray(freq_all)
print(np.shape(freq))
# coor = np.asarray(coor_all)
radi_all = []
lati_all = []
long_all = []
for i in range(0, len(wnas_all)):
    radi_all.append(coor_all[i].radi)
    lati_all.append(coor_all[i].lati)
    long_all.append(coor_all[i].long)
radi = np.asarray(radi_all)
lati = np.asarray(lati_all)
long = np.asarray(long_all)


# sns.set(style=“ticks”)
fig, (ax0) = plt.subplots(1, 1, figsize=(8,5), dpi=250, constrained_layout=True)
ax0.scatter(long, radi)#, s=2, c=np.rad2deg(wnas), cmap=‘gnuplot2’)
ax0.set_ylabel('Earth Radii')
ax0.set_xlabel('Longitude')
# cb = fig.colorbar(c, ax=ax0)
# cb.set_label(r’$\theta$ [deg]’)
plt.show()

# sns.set(style=“ticks”)
fig, (ax0) = plt.subplots(1, 1, figsize=(8,5), dpi=250, constrained_layout=True)
ax0.scatter(lati, radi)#, s=2, c=np.rad2deg(wnas), cmap=‘gnuplot2’)
ax0.set_ylabel('Earth Radii')
ax0.set_xlabel('Latitude')
# cb = fig.colorbar(c, ax=ax0)
# cb.set_label(r’$\theta$ [deg]’)
plt.show()

# plt.style.use(‘seaborn’)
# import seaborn as sns
plt.figure(figsize=(8,4), dpi = 250)
plt.hist(np.rad2deg(wnas), bins=45)
# plt.yscale('log')
plt.title('Raytracer WNA Distribution')
plt.show()