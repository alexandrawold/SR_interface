# import needed packages
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import datetime as dt
from constants_settings import *
from convert_coords import convert2
from raytracer_utils import readdump
from bfield import getBline
import os 
from spacepy import coordinates as coord
from spacepy.time import Ticktock
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap
from satellites import sat

#------------------------------- rotation for plots --------------------------------------------
def rotateplane(rc, tvec_datetime, crs, carsph, units):
    # rotate to be at prime meridian for plotting purposes
    rot_crs = convert2(rc, tvec_datetime, crs, carsph, units, crs, 'sph', ['Re','deg','deg'])

    rot_crs_new = [[rcs[0], rcs[1], 0] for rcs in rot_crs]
    final_crs = convert2(rot_crs_new, tvec_datetime, crs,'sph', ['Re','deg','deg'], crs, carsph, units)

    return final_crs
#---------------------------------------------------------------------------


# ---------------------------------------- STIX PARAM --------------------------------------------
def stix_parameters(ray, t, w):
    # from austins code

    # get b field
    B   =  ray['B0'].iloc[t]
    Bmag = np.linalg.norm(B)
    # get plasma species
    Q    = np.abs(np.array(ray['qs'].iloc[t,:]))
    M    = np.array(ray['ms'].iloc[t,:])
    Ns   = np.array(ray['Ns'].iloc[t,:])
    # cyclotron and plasma freq
    Wcs   = Q*Bmag/M
    Wps2  = Ns*pow(Q,2)/EPS0/M
    # stix param calculations
    R = 1.0 - np.sum(Wps2/(w*(w + Wcs)))
    L = 1.0 - np.sum(Wps2/(w*(w - Wcs)))
    P = 1.0 - np.sum(Wps2/(w*w))
    S = (R+L)/2.0
    D = (R-L)/2.0

    return R, L, P, S, D
# ---------------------------------------------------------------------------------------------


#-----------------------------  L SHELL TRACING  ----------------------------------------------
def get_Lshells(lshells, tvec_datetime, crs, carsph, units):
    Lshell_fline = []
    for lsh in lshells:
        # needs to be cartesian, output will be in Re
        # plot l shell as if at equator and prime merid
        newt = convert2([[lsh*R_E,0,0]], tvec_datetime, crs, 'sph', ['m','deg','deg'], 'SM', 'car', ['m','m','m'])
        T = getBline(newt[0], tvec_datetime[0],100)
        
        # repack
        T_repackx = T.x
        T_repacky = T.y
        T_repackz = T.z
        T_repack = [[tx, ty, tz] for tx,ty,tz in zip(T_repackx, T_repacky, T_repackz)]
        LshellT = [tvec_datetime[0] for i in range(len(T_repack))]
        
        T_convert = convert2(T_repack, LshellT, 'SM','car',['Re','Re','Re'], crs, carsph, units)
        
        T_repackx = [tt[0] for tt in T_convert]
        T_repackz = [tt[2] for tt in T_convert]
        
        Lshell_fline.append([T_repackx, T_repackz])
    return Lshell_fline
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
def plotray2D(ray_datenum, raylist, ray_out_dir, crs, carsph, units, md, t_save=None, show_plot=True, damping_vals=None):

    # !!!!!!! ONLY suports cartesian plotting!
    # make figures directory if doesnt exist
    if not os.path.exists(ray_out_dir+'/figures'):
        os.makedirs(ray_out_dir+'/figures')

    # first, find the max index of refraction (useful for ray spot plots)
    findmax_deg = []
    # convert to desired coordinate system into vector list rays
    ray_coords = []

    # loop through the rays -- NEED ALL THE RAYS
    for ri, r in enumerate(raylist):
        B   =  r['B0'].iloc[0]
        Bmag = np.linalg.norm(B)
        bunit = B/Bmag

        # get the initial wna (already confirmed this is equiv to the initial)
        w = r['w']
        kveci = [(w/C)*r['n'].x.iloc[0],(w/C)*r['n'].y.iloc[0],(w/C)*r['n'].z.iloc[0]]
        kunit = np.array(kveci)/np.linalg.norm(kveci)

        # alpha is wna
        alpha = np.arccos(np.dot(kunit, bunit))
        alphaedg = float(alpha)*R2D
        findmax_deg.append(alphaedg)
        
        # only want to plot 100 rays to visualize
        if ri < 100:
            # comes in as SM car in m 
            tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
            nmag = np.sqrt(r['n'].x.iloc[0]**2 +r['n'].y.iloc[0]**2+r['n'].z.iloc[0]**2)
            tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
            new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)

            # save it
            ray_coords.append(new_coords) 

    # find max 
    max_wna = max(findmax_deg)  
    max_wna_ind = findmax_deg.index(max(findmax_deg))
    rm = raylist[max_wna_ind]
    nmax = np.sqrt(rm['n'].x.iloc[0]**2 +rm['n'].y.iloc[0]**2+rm['n'].z.iloc[0]**2)
    
    # -------------------------------- PLOTTING --------------------------------
    ray1 = ray_coords[0][0]
    # this is GEO for the Earth plot
    ray1_sph = convert2([ray1], tvec_datetime, crs, carsph, units, 'GEO', 'sph', ['m','deg','deg'])

    # find ray starting long just for Earth
    long_r = ray1_sph[0][2]
    # convert for cartopy issues
    if long_r > 180:
        plane_long = 180-(long_r-180)
        catch = -1
    else:
        plane_long = long_r
        catch = 1
    try:
        import cartopy
        import cartopy.crs as ccrs
        # make this into a an image (sneaky)

        # chunk of code for the Earth
        ax = plt.axes(projection=ccrs.Orthographic(central_longitude=(catch*plane_long)-90, central_latitude=0.0))
        ax.add_feature(cartopy.feature.OCEAN, zorder=0)
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
        ax.coastlines()
        plt.savefig(ray_out_dir+'/figures/ccrs_proj.png')
        plt.cla()

        import matplotlib.patches as patches
        import matplotlib.cbook as cbook

        img = plt.imread(ray_out_dir+'/figures/ccrs_proj.png', format='png')
        fig, ax = plt.subplots(figsize=(12,12))
        im = ax.imshow(img, extent=[-1.62,1.62,-1.3,1.3],zorder=2) # not the most scientific
        patch = patches.Circle((0, 0), radius=1.02, transform=ax.transData)
        im.set_clip_path(patch)

    except: # or just a patch for the Earth
        fig, ax = plt.subplots(figsize=(12,12))
        earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=100)
        iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)
        ax.add_artist(earth)
        ax.add_artist(iono)

    # plot the rays
    rotated_rcoords = []
    # rotate to be in plane (XZ plane)
    for rayc in ray_coords:
        for rc in rayc:
            rotated = rotateplane([rc], tvec_datetime, crs, carsph, units)
            rotated_rcoords.append(rotated[0])

    # repack in xz plane
    rotated_rcoords_x = [rcs[0] for rcs in rotated_rcoords]
    rotated_rcoords_z = [rcs[2] for rcs in rotated_rcoords]

    # plot!
    ax.scatter(rotated_rcoords_x[0], rotated_rcoords_z[0], c = 'Blue', s=10)
    ax.scatter(rotated_rcoords_x, rotated_rcoords_z, s = 1, zorder = 103)
    #if t_save:
    #    ax.scatter(rotated_rcoords_x[t_save], rotated_rcoords_z[t_save], c = 'r', s=20, zorder=104)

    # final clean up
    # plot field lines (from IGRF13 model)
    L_shells = [2, 3, 4]  # Field lines to draw
    Lshell_flines = get_Lshells(L_shells, tvec_datetime, crs, carsph, units)
    
    for lfline in Lshell_flines:
        ax.plot(lfline[0], lfline[1], color='Black', linewidth=1, linestyle='dashed',zorder=3)

    # find DSX
    dsx = sat()             # define a satellite object
    dsx.catnmbr = 44344     # provide NORAD ID
    dsx.time = ray_datenum  # set time
    dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

    dsx.propagatefromTLE(sec=0, orbit_dir='future', crs=crs, carsph='car', units=['Re','Re','Re'])
    dsx_spot = dsx.pos
    rotated_dsx = rotateplane(dsx_spot, tvec_datetime, crs, carsph, units)
    plt.scatter(-1*dsx_spot[0][0],dsx_spot[0][2],marker='*',s=20,zorder=107)

    plt.xlabel('L (R$_E$)')
    plt.ylabel('L (R$_E$)')
    #plt.xlim([0, max(L_shells)])
    plt.xlim([0, 3])
    plt.ylim([-2, 2])
    
    if t_save:
        plt.savefig(ray_out_dir + '/figures/' + dt.datetime.strftime(ray_datenum, '%Y_%m_%d_%H%M%S') + '_' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz' +'_2Dview' + str(md) + '_'+str(t_save)+ '.png',bbox_inches='tight')
    else:    
        plt.savefig(ray_out_dir + '/figures/' + dt.datetime.strftime(ray_datenum, '%Y_%m_%d_%H%M%S') + '_' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz' +'_2Dview' + str(md) + '.png',bbox_inches='tight')

    if show_plot:
        plt.show()
    plt.close()

    return nmax
# --------------------------------------------------------------------------------------


# ---------------------------------- PLOT REFRACTIVE SURFACE ----------------------------------------------------
def plotrefractivesurface(ray_datenum, ray, ray_out_dir, t_save):
    # set up phi vec
    phi_vec = np.linspace(0,360,int(1e5))*D2R
    w = ray['w']
    
    # comes in as SM car in m 
    tmp_kcoords = list(zip((w/C) * ray['n'].x, (w/C) * ray['n'].y, (w/C) * ray['n'].z))

    # convert to a unit vector first
    unitk = [(tmp_kcoords[s][0], tmp_kcoords[s][1], tmp_kcoords[s][2]) / np.sqrt(tmp_kcoords[s][0]**2 + tmp_kcoords[s][1]**2 + tmp_kcoords[s][2]**2) for s in range(len(tmp_kcoords))]
    kcoords = unitk

    #int_plt = len(kcoords)//10

    for ti, t in enumerate(ray['time']):
        #if ti % int_plt == 0:
        if ti == t_save: 
            print('plotting a refractive surface')
            # get stix param
            R, L, P, S, D = stix_parameters(ray, ti, w)
            root = -1 # whistler solution

            eta_vec = np.zeros_like(phi_vec)

            # solution from antenna white paper!
            #resangle = np.arctan(np.sqrt(-P/S)) 
            #print(resangle)
            # find phi
            B   =  ray['B0'].iloc[ti]
            Bmag = np.linalg.norm(B)

            wce = Q_EL*Bmag/M_EL

            gendrin_angle = 180 - (R2D*np.arccos(2*w/wce))

            bunit = B/Bmag

            kunit = np.array([kcoords[ti][0], kcoords[ti][1], kcoords[ti][2]])

            alpha = np.arccos(np.dot(kunit, bunit))
            alphaedg = float(alpha)*R2D

            alpha_diff_last = 1

            for phi_ind, phi  in enumerate(phi_vec):

                # Solve the cold plasma dispersion relation
                cos2phi = pow(np.cos(phi),2)
                sin2phi = pow(np.sin(phi),2)

                A = S*sin2phi + P*cos2phi
                B = R*L*sin2phi + P*S*(1.0+cos2phi)

                discriminant = B*B - 4.0*A*R*L*P
                n1sq = (B + np.sqrt(discriminant))/(2.0*A)
                n2sq = (B - np.sqrt(discriminant))/(2.0*A)

                # negative refers to the fact that ^^^ B - sqrt
                n1 = np.sqrt(n1sq)
                n2 = np.sqrt(n2sq)

                eta_vec[phi_ind] = n2

                # save the angle near alpha (easier for plotting)
                alpha_diff = np.abs(alpha - phi)
                if alpha_diff < alpha_diff_last:
                    phi_ind_save = phi_ind
                alpha_diff_last = alpha_diff

            # ---------- plot it --------------------------
            fig, ax = plt.subplots(1,1, figsize=(8,5))

            # plot the surface
            ax.plot(eta_vec*np.sin(phi_vec), eta_vec*np.cos(phi_vec), 'gray', LineWidth = 1, label = 'e + ions')
            
            quiver_mag = eta_vec[phi_ind_save]
            ax.plot([0, quiver_mag*np.sin(float(alpha))], [0, quiver_mag*np.cos(float(alpha))], 'r--')

            xlim1 = -100
            xlim2 = -xlim1
            ylim1 = -100
            ylim2 = -ylim1

            ax.set_xlim([xlim1, xlim2])
            ax.set_ylim([ylim1, ylim2])
            ax.set_xlabel('Transverse Refractive Component')
            #plt.legend(loc='upper right')
            ax.text(60, 65, r'$\alpha$' +  ' = ' + str(round(alphaedg,2)))
            ax.text(60, 55, r'$\beta$' +  ' = ' + str(round(gendrin_angle,2)))
            
            #plt.title(dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S') + ' ' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz \n refractive surface ')
            plt.savefig(ray_out_dir+'/refractive_surface'+str(ti)+'.png')
            plt.clf()
            plt.close()


    return

# ------------------------------------------- END --------------------------------------------


# --------------------------- plot geometric factor along ray path! ----------------------------------------------------
def plotgeomfactor(ray_datenum, raylist, ray_out_dir, crs, carsph, units, show_plot=True):

    # convert to desired coordinate system into vector list rays
    ray_coords = []
    for r in raylist:
        w = r['w']

        # comes in as SM car in m 
        tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        
        # convert to a unit vector first
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        seconds_passed = [s for s in r['time']]

        new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)

        # save it
        ray_coords.append(new_coords)
        
    # -------------------------------- PLOTTING --------------------------------
    # flatten
    ray_coords = ray_coords[0]
    r_init = ray_coords[0][0]
    lat_init = ray_coords[0][1]
    
    fig, ax = plt.subplots(figsize=(6,6))
    
    for ti, rc in enumerate(ray_coords):
        # find geometric factor
        gf = r_init * np.cos(lat_init*D2R) / (rc[0] * np.cos(rc[1]*D2R))
        ax.scatter(seconds_passed[ti],gf,color='b')

    plt.xlabel('seconds')
    #plt.gcf().autofmt_xdate()
    plt.ylabel('Geom. Factor')
    plt.title(dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S') + ' ' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz')
    plt.savefig(ray_out_dir + '/' + dt.datetime.strftime(ray_datenum, '%Y_%m_%d_%H%M%S') + '_' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz' +'_geomfac_RE' + '.png')

    if show_plot:
        plt.show()

    return
# ------------------------------------------- END --------------------------------------------


# --------------------------- plot Ne along ray path! ----------------------------------------------------
def plotNe(raylist):
    
    for ri,r in enumerate(raylist):
        Ne = []
        for ti, t in enumerate(r['time']):
            edens = np.array(r['Ns'].iloc[ti])
            Ne.append(np.log10(edens[0]))
        plt.plot(r['time'],Ne,label=['e - model '+str(ri)])

    plt.xlabel('time [sec]')
    plt.xlim([0,0.5])
    plt.ylabel('log(density)')
    plt.show()
    plt.close()

# --------------------------------------------------------------------------------------


# ------------------------------------------- plot plasmasphere density (2D plot) --------------------------------------------
def plot_plasmasphere_2D(md,kp):
    fig, ax = plt.subplots(1,1)

    model_path = 'modeldumps/'
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_XY_'+str(kp)+'.dat')
    d_xy = readdump(plasma_model_dump)
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_XZ_'+str(kp)+'.dat')
    d_xz = readdump(plasma_model_dump)
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_YZ_'+str(kp)+'.dat')
    d_yz = readdump(plasma_model_dump)
    Ne_xy = d_xy['Ns'][0,:,:,:].squeeze().T
    Ne_xy[np.isnan(Ne_xy)] = 0
    Ne_xz = d_xz['Ns'][0,:,:,:].squeeze().T
    Ne_xz[np.isnan(Ne_xz)] = 0
    Ne_yz = d_yz['Ns'][0,:,:,:].squeeze().T
    Ne_yz[np.isnan(Ne_yz)] = 0

    # Axis spacing depends on how the modeldump was ran
    psize = 10
    px = np.linspace(-10, 10, 200)
    py = np.linspace(-10, 10, 200)

    # Colorbar limits (log space)
    clims = [-2, 5]

    # Plot background plasma 
    g = plt.pcolormesh(px, py, np.log10(Ne_xz*10**(-6)), vmin=0,vmax=4, cmap = 'jet')

    fig.colorbar(g, orientation="horizontal", pad = 0.1, label= 'Plasmasphere density #/cm^3')
    earth = plt.Circle((0,0),1,color='k',alpha=1, zorder=100)
    ax.add_patch(earth)   
    plt.ylim([-4,4])
    plt.xlim([0,8])
    ax.set_aspect('equal')

    plt.show()
    plt.close()

# ------------------------------------------- END --------------------------------------------


# ------------------------------------------- plot plasmasphere 1D --------------------------------------------
def plot_plasmasphere_1D(md,kp):
    fig, ax = plt.subplots(1,1,figsize=(4,6))
    # Plot background plasma 
    import seaborn as sns

    # Apply the default theme
    sns.set_theme()

    model_path = 'modeldumps/'
    #for kp in kps:
    #for md in mds:
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_XY_'+str(kp)+'.dat')
    d_xy = readdump(plasma_model_dump)
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_XZ_'+str(kp)+'.dat')
    d_xz = readdump(plasma_model_dump)
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_YZ_'+str(kp)+'.dat')
    d_yz = readdump(plasma_model_dump)

    labels=  ['e','o+','h']
    cs = ['b','r','g','r']

    for i in range(3):
        print(i)
        Ne_xy = d_xy['Ns'][i,:,:,:].squeeze().T
        Ne_xy[np.isnan(Ne_xy)] = 0
        Ne_xz = d_xz['Ns'][i,:,:,:].squeeze().T
        Ne_xz[np.isnan(Ne_xz)] = 0
        Ne_yz = d_yz['Ns'][i,:,:,:].squeeze().T
        Ne_yz[np.isnan(Ne_yz)] = 0
        # Axis spacing depends on how the modeldump was ran
        psize = 10
        px = np.linspace(-10, 10, 200)
        py = np.linspace(-10, 10, 200)

        plt.plot(-px,np.log10(Ne_yz[99,:]),cs[i])

    plt.xlabel('L shell')
    plt.ylabel('Log Equatorial Density (#/cm^3)')

    plt.xlim([0.5,6])
    plt.ylim([-1,12])
    #plt.legend(['ngo','diff eq','gcpm'])
    #plt.legend([0,4])
    plt.legend(labels)
    plt.grid()
    plt.show()
    plt.close()

# ------------------------------------------- END --------------------------------------------
#plot_plasmasphere_1D(7,1)


def plot_density_alongpath(ray_datenum, ray, ray_out_dir,t_save):
    fig = plt.figure(figsize=(3.5,5))
    xs = np.linspace(0,1,len(ray['time']))
    #print(ray['time'],len(ray['time']))
    for ti,t in enumerate(ray['time']):
        raydens = ray['Ns'].iloc[ti]
        Ne = raydens[28]
        Nps = raydens[29]
        Nhe = raydens[30]
        Nos = raydens[31]

        plt.scatter(xs[ti],np.log10(Ne),c='k',s=10)
        plt.scatter(xs[ti],np.log10(Nps),c='y',s=10)
        plt.scatter(xs[ti],np.log10(Nhe),c='g',s=10)
        plt.scatter(xs[ti],np.log10(Nos),c='b',s=10)

        if ti==t_save:
            plt.scatter(xs[ti],np.log10(Ne),c='r',zorder=100)
            plt.scatter(xs[ti],np.log10(Nps),c='r',zorder=100)
            plt.scatter(xs[ti],np.log10(Nhe),c='r',zorder=100)
            plt.scatter(xs[ti],np.log10(Nos),c='r',zorder=100)

    plt.ylabel('log density #/m^3')
    #plt.legend(['electrons (protons)', 'helium ions', 'oxygen ions'])
    #plt.ylim([6,9.75])
    #plt.show()
    plt.savefig(ray_out_dir+'/density_path'+str(t_save)+'.png')
    plt.close()