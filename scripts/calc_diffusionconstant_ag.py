"""Script to convert HISTORY file
AG, February 25, 2013
AG, April 25, 2013

script to determine which atoms are adsorbates in the zeolite simulations
AG, June 24, 2013

correction: do unwrapping
AG, Sept 10, 2013

DATA
time data from HISTORY
unitcell from h5
trajectory from h5
atomtypes from HISTORY
molecule sizes from FIELD
"""

import numpy as np
import saman
from saman.dlpolyfunc import *

# for matplotlib figure generation...
import matplotlib
matplotlib.use('Agg')

######################################################################

##############
# INPUT
##############

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(prog='calc_diffusionconstant_ag',
              description='Extracts the diffusion constant from DL_POLY files.')
    parser.add_argument('--field',dest='filename_field',
                        help='FIELD file from DL_POLY')
    parser.add_argument('--h5',dest='filename_h5',nargs='*',
                        help='h5 file, converted from DL_POLY')
    parser.add_argument('--history',dest='filename_history',
                        help='HISTORY file from DL_POLY')
                     # TODO can be done with CONFIG?
    parser.add_argument('-c','--combine',dest='combine',
                        default='perfile',
                        help='how to combine trajectories when taking the average: over molecules, over files, or just cumulate the files.') 
    parser.add_argument('-t','--times',dest='times',nargs=3,
                        help='t_start,t_end,t_delta in ps')
    parser.add_argument('--atomtypes',dest='atomtypes',nargs='*',
                        help='selected atomtypes')
    parser.add_argument('--unfolded',dest='unfolded',action='store_true',
                        default=False,
                        help='whether the trajectory is already unfolded [default=False (not unfolded yet)]')
    parser.add_argument('--histogram',dest='filename_histogram',
                        default='./hist',
                        help='filename (base) for histograms')
    #>>> parser = argparse.ArgumentParser(prog='PROG')
    #>>> parser.add_argument('--foo')
    #>>> parser.add_argument('command')
    #>>> parser.add_argument('args', nargs=argparse.REMAINDER)
    #>>> print parser.parse_args('--foo B cmd --arg1 XX ZZ'.split())
    #Namespace(args=['--arg1', 'XX', 'ZZ'], command='cmd', foo='B')
    #import os
    #cwd = os.getcwd()
    #index = cwd[::-1].index("/")
    #indir = cwd[:-index]+"/../results/"+cwd[-index:]
    
    args = parser.parse_args()
    print "parsed"
    print args.__dict__
    
    # define filenames
    filename_field = args.filename_field
    list_filename_h5 = args.filename_h5
    if len(list_filename_h5) == 1:
       list_filename_h5 = list_filename_h5[0].split()
    print "list_filename_h5",list_filename_h5
    filename_history = args.filename_history
    combine = args.combine # how to average over molecules/trajectories

    # whose MSD is taken into account, e.g. ['C2','H1'], or ['C2'], or ...
    selected_atomtypes = args.atomtypes
    unfolded = args.unfolded
    times = [float(t) for t in args.times]

    ######################################################################
    do_calctensorD = True
    do_gethistograms = False
    do_getksi = False
    ######################################################################
    
    ##############
    # CALCULATE
    ##############

    # read data
    list_allcoor,unitcell,seladsorbates = get_xyz(filename_field,list_filename_h5,filename_history,combine,selected_atomtypes,atoms="adsorbates")
    t_start,t_end,t_delta,t0,t1,dt,dtc,dn1,dn2,ddn = get_times_MSD(times,filename_history)


    if do_getksi:
       pass 


    if do_gethistograms:
        from mcdiff.tools.functionshistogram import plot_histogram
        filename_histogram = args.filename_histogram
    
        def create_projected(vec,list_x,list_y,list_z):
            # project
            lenv = np.linalg.norm(vec)
            v = vec/lenv
            unitcell_v[:,i,i] = lenv
            proj = [v[0]*x+v[1]*y+v[2]*z for (x,y,z) in zip(list_x,list_y,list_z)]
            # reduce to one box
            #proj = [p-lenv*np.floor(p/lenv+0.5) for p in proj]
            #unfolded = False
            return proj
    
        unitcellvecs = unitcell[0,:,:]
        projs = []
        unitcell_v = np.zeros(unitcell.shape,float)
        for i in range(3):
            vec = unitcellvecs[:,i]
            proj = create_projected(vec,list_x,list_y,list_z)
            unitcell_v[:,i,i] = np.linalg.norm(vec)
            projs.append(proj)

        # create variations
        for (i,j,k) in [[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1],[1,1,1],[1,1,-1],[1,-1,1],[-1,1,1]]:
            vec = unitcellvecs[:,0]*i + unitcellvecs[:,1]*j + unitcellvecs[:,2]*k
            proj = create_projected(vec,list_x,list_y,list_z)
            projs.append(proj)

        kwargs = {"ymin":0,"ymax":6}
        for j,proj in enumerate(projs):
            # these are not reduced to one box
            c = []
            for i in range(len(proj)): c.extend(proj[i].ravel().tolist())
           #print max(np.array(c).ravel()),min(np.array(c).ravel())
            plot_histogram(np.array(c),"%s.proj%i"%(filename_histogram,j),**kwargs)

        if True:
            # compute MSD
            outdir = "./"
            from mcdiff.tools.functionsdistance import analyze_matrixdist

            if unfolded:
                # faster routine than when unfolding is still necessary
                analyze_matrixdist(projs[0],projs[1],projs[2],dn1,outdir,dtc,dn2=dn2,ddn=ddn,)
            else:
                # this will first unfold the trajectory based on the unitcells
                analyze_matrixdist(projs[0],projs[1],projs[2],dn1,outdir,dtc,dn2=dn2,ddn=ddn,unitcell=unitcell)



#    unitcellvecs = unitcell[0,:,:].transpose()
#    projs = []
#    for i,vec in enumerate(unitcellvecs):
#        # project
#        v = vec/np.linalg.norm(vec)
#        proj = [v[0]*x+v[1]*y+v[2]*z for (x,y,z) in zip(list_x,list_y,list_z)]
#        projs.append(proj)
#
#    kwargs = {"ymin":0,"ymax":6}
#    for j,proj in enumerate(projs):
#        c = []
#        for i in range(len(proj)): c.extend(proj[i].ravel().tolist())
#        plot_histogram(np.array(c),"%s.projT%i"%(filename_histogram,j),**kwargs)


    if False:
        # create histogram
        from mcdiff.tools.functionshistogram import plot_histogram
        filename_histogram = args.filename_histogram
        kwargs = {"ymin":0,"ymax":6}
        c = []
        for i in range(len(list_x)): c.extend(list_x[i].ravel().tolist())
        plot_histogram(np.array(c),"%s.x"%(filename_histogram),**kwargs)
        c = []
        for i in range(len(list_y)): c.extend(list_y[i].ravel().tolist())
        plot_histogram(np.array(c),"%s.y"%(filename_histogram),**kwargs)
        c = []
        for i in range(len(list_z)): c.extend(list_z[i].ravel().tolist())
        plot_histogram(np.array(c),"%s.z"%(filename_histogram),**kwargs)
        #plot_histogram(np.array(list_x).ravel(),"%s.x"%(filename_histogram))
        #plot_histogram(np.array(list_y).ravel(),"%s.y"%(filename_histogram))
        #plot_histogram(np.array(list_z).ravel(),"%s.z"%(filename_histogram))


    if do_calctensorD:
        # compute MSD
        outdir = "./"
        from mcdiff.tools.functionsdistance import analyze_matrixdist
    
        if unfolded:
            # faster routine than when unfolding is still necessary
            #analyze_matrixdist(list_x,list_y,list_z,dn1,outdir,dtc,dn2=dn2,ddn=ddn,)
            analyze_matrixdist([arr[:,:,0] for arr in list_allcoor],
                 [arr[:,:,1] for arr in list_allcoor],[arr[:,:,2] for arr in list_allcoor],
                 dn1,outdir,dtc,dn2=dn2,ddn=ddn,)
        else:
            # this will first unfold the trajectory based on the unitcells
            analyze_matrixdist([arr[:,:,0] for arr in list_allcoor],
                 [arr[:,:,1] for arr in list_allcoor],[arr[:,:,2] for arr in list_allcoor],
                 dn1,outdir,dtc,dn2=dn2,ddn=ddn,unitcell=unitcell) 
    
