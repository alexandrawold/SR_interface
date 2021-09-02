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


# example call to ray tracer!
# --------------------------------------- set up ------------------------------------------
rayfile_directory = '/home/alwo2026/Workspace/SR_output' # store output here

# FIRST, navigate to constants_settings and make sure the settings are correct for the run

# picking start position in GEO sph and converting to SM car

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

    # how many rays? 
    nrays = len(freqs) # how many rays

    # next, define the direction of the ray
    makedate = dt_array[0].strftime('%Y%m%d')
    Date = int(makedate)
    ut = dt_array[0].hour + dt_array[0].minute/60 + dt_array[0].second/3600

    Bx,By,Bz = gp.ModelField(ray_start[0][0]/R_E,ray_start[0][1]/R_E,ray_start[0][2]/R_E,Date,ut,Model='T96',CoordIn='SM',CoordOut='SM')
    Bdir = np.array([Bx, By, Bz])
    Bunit = Bdir/np.linalg.norm(Bdir)
    if pos_array[0][1] > 0:
        Bsouth = [-1*float(Bunit[0]), -1*float(Bunit[1]), -1*float(Bunit[2])]
    else:
        Bsouth = [float(Bunit[0]), float(Bunit[1]), float(Bunit[2])]


    # parallelizing babeeeyyy
    nworkers = 8 #max cores = 16

    freqs_list = [freqs[int(i * (nrays/nworkers)):int((i+1)*nrays/nworkers)] for i in range(nworkers)]
    print('RUNNING ', nrays, nrays/nworkers)

    # same freq and starting position for all
    directions_list = [[Bsouth for p in range(len(d))] for d in freqs_list]
    positions = [[ray_start[0] for p in range(len(d))] for d in freqs_list]

    tvec = [dt_array[0] for n in range(nworkers)]
    directory_list = [rayfile_directory for i in range(len(tvec))]
    mds = [md for i in range(len(tvec))]

    fn_str = ['lat%s_lon%s'%(pos_array[0][1], pos_array[0][2]) for i in range(len(tvec))]
    print(fn_str)

    parallel_run_rays(tvec, positions, directions_list, freqs_list, directory_list, mds, fn_str)

# #  --------------------------- get ray file output ---------------------------
# ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(dt_array[0], '%Y-%m-%d_%H_%M_%S')

# file_titles = os.listdir(ray_out_dir)

# raylist = []
# mode_name = 'mode' + str(md)
# # use mode name to avoid workers of the same label
# x = sorted(file_titles)

# for filename in x:
#     if '.ray' in filename and mode_name in filename:
#         raylist += read_rayfile(os.path.join(ray_out_dir, filename))
#         print(filename)

# print(raylist[0]['w'])

# # plotting

# for r in raylist:
#     plotray2D(dt_array[0], [r], ray_out_dir, crs='GEO', carsph='car', units=['Re','Re','Re'], md=md, checklat=pos_array[0][1], t_save=None, plot_wna=False, show_plot=True, damping_vals=None)

