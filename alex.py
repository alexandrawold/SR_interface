import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from run_rays import single_run_rays, parallel_run_rays
from ray_plots import plotray2D
from bfield import getBdir
from raytracer_utils import read_rayfile, read_input_jobs, read_damp_matlab, read_bigrayfile, read_bigrayfile_in
from constants_settings import *
from convert_coords import convert2
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import FormatStrFormatter
import PyGeopack as gp

md = 6
ray_datenum = dt.datetime(2014,1,1,12,0, tzinfo=dt.timezone.utc)
frequencies = [ 0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. , 5.5,  6. ,  6.5,  7. ,  7.5,  8. ,  8.5,  9. ,  9.5, 10. , 10.5, 11. , 11.5, 12. ]
rayfile_directory = '/home/rileyannereid/workspace' # store output here

makedate = ray_datenum.strftime('%Y%m%d')
Date = int(makedate)
ut = ray_datenum.hour + ray_datenum.minute/60 + ray_datenum.second/3600

for freq in frequencies:
    print(freq)
    latitudes = [-70., -63., -56., -49., -42., -35., -28., -21., -14.,  -7.,   0., 7.,  14.,  21.,  28.,  35.,  42.,  49.,  56.,  63.,  70.]
    longitudes = [-180., -160., -140., -120., -100.,  -80.,  -60.,  -40.,  -20.,  0.,   20.,   40.,   60.,   80.,  100.,  120.,  140.,  160., 180.]
    freqs = [freq*1e3 for n in latitudes] # create list of frequencies for each
    positions = []
    directions = []
    for la, lo in zip(latitudes,longitudes):
        start_pt = convert2([[R_E+500e3, la,lo]], [ray_datenum],'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])
        positions.append(start_pt[0])
    
        # get Bfield direction
        Bx,By,Bz = gp.ModelField(start_pt[0][0]/R_E,start_pt[0][1]/R_E,start_pt[0][2]/R_E,Date,ut,Model='T96',CoordIn='SM',CoordOut='SM')
        Bdir = np.array([Bx, By, Bz])
        Bunit = Bdir/np.linalg.norm(Bdir)
        if la > 0:
            Bsouth = [-1*float(Bunit[0]), -1*float(Bunit[1]), -1*float(Bunit[2])]
        else:
            Bsouth = [float(Bunit[0]), float(Bunit[1]), float(Bunit[2])]
            
        directions.append(Bsouth)
    
    #single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md, runmodeldump=True)

# that's it! let's look at output
for freq in frequencies:
    ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')+ '/' + str(int(freq*1e3/10))
    file_titles = os.listdir(ray_out_dir)
    mode_name = 'mode'+str(md)

    nraylist = []
    x = sorted(file_titles)
    for filename in x:
        if '.ray' in filename and mode_name in filename:
            nraylist += read_rayfile(os.path.join(ray_out_dir, filename))
            print(filename)
    for ri,rr in enumerate(nraylist):
        plotray2D(ray_datenum, [rr], ray_out_dir, 'GEO', 'car', ['Re','Re','Re'], md, show_plot=False,plot_density=False,damping_vals=None,t_save=ri)
