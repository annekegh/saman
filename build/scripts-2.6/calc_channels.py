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
from saman.readdlpolyfunc import *

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
    # RINGS
    parser.add_argument('--rings',dest='ringdir',
                        default='./',
                        help='location of files that contain ring definitions')
    
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

    ##############
    # CALCULATE
    ##############

    print "*"*20
    print "Create RingGroup"
    print "*"*20
    # STRUCTURE
    from saman.ringfunc import read_coords,read_rings,read_unitcell
    from saman.ringfunc import RingGroup
    from saman.ringfunc import check_location_bis, get_centers_planes, write_ksi_select

    # create ring using config.xyz
    #fn_xyz = "/user/home/gent/vsc400/vsc40051/mydata/germansastre/channels/ringdefinitions/sapo34/config.xyz"
    #fn_rings = "/user/home/gent/vsc400/vsc40051/mydata/germansastre/channels/ringdefinitions/sapo34/ringatoms"
    fn_xyz = "%s/config.xyz"%args.ringdir
    fn_rings = "%s/ringatoms"%args.ringdir
    fn_unitcell = "%s/cellvectors"%args.ringdir
    natom,atomtypes,pos = read_coords(fn_xyz)
    list_rings_indices = read_rings(fn_rings)
    unitcell0 = read_unitcell(fn_unitcell)
    print "unitcell0",unitcell0

    ####### !!!! while testing
    ####list_rings_indices = list_rings_indices[:2]

    # create ring group (rings that are split by pbc get new positions)
    rg = RingGroup(list_rings_indices,pos,atomtypes,unitcell=unitcell0)
    #rg.set_ellips()

    nfile = len(list_filename_h5)
    print "nfile",nfile
    for n in range(nfile):
        # do file one per one
        # extract run number:  format name h5 is xxx.traj9.h5 or xxx.unfold-traj9.h5
        print "doing file",n,list_filename_h5[n]
        i = list_filename_h5[n].find("traj")    # start of word "traj"
        j = list_filename_h5[n].find(".",i)    # will be start of ".h5" extension
        run = int(list_filename_h5[n][i+4:j])
        print "extracted run:",run

        # read data framework
        selected_atomtypes_framework = ['Si','P','O1','O2','Al','H1']
        list_allcoorframework,unitcell,selframework = get_xyz(
                 filename_field,[list_filename_h5[n]],filename_history,combine,selected_atomtypes_framework,
                 atoms='framework')

        # read data adsorbates
        list_allcoor,unitcell,seladsorbates = get_xyz(
                 filename_field,[list_filename_h5[n]],filename_history,combine,selected_atomtypes,)
        t_start,t_end,t_delta,t0,t1,dt,dtc,dn1,dn2,ddn = get_times_MSD(times,filename_history)


        assert len(list_allcoor)==1
        assert len(list_allcoorframework)==1
        print "list_allcoorframework[0]",list_allcoorframework[0].shape
        print "list_allcoor[0]",list_allcoor[0].shape


        select = np.arange(0,2*len(seladsorbates),2)  # just take one atom of both
        nsel = len(select)
        print "nsel",nsel

        print "*"*20
        print "Start ring analysis"
        print "*"*20

        ntime = len(list_allcoor[0])
        print "ntime",ntime
        centers,planes = get_centers_planes(list_allcoorframework[0],rg,atomtypes,unitcell=unitcell0)
        #print "centers,planes",centers.shape,planes.shape  #ntime x 3
        for i,ring in enumerate(rg.list_rings):
            print "---doing ring%i"%i
            pos = np.take(list_allcoor[0],select,1)  # ntime x natom x 3
            dist = np.zeros((ntime,nsel),float)
            ksi = np.zeros((ntime,nsel),float)
            for at in range(nsel):
                d,k = check_location_bis(centers[:,i,:],planes[:,i,:],pos[:,at,:],unitcell[0,:,:])   #unitcell of first time step
                dist[:,at] = d[:]
                ksi[:,at] = k[:]

            fn_ksi = "ksi.%i.run%i.dat"%(i,run)
            import numpy.ma as ma
            dist_ma = ma.masked_greater(dist,5.)
            # detect ring passage
            sign = (ksi>0)    # True if ksi>0, False if ksi<0
            signchange = np.array(sign[1:,:],int)-np.array(sign[:-1,:],int)

            write_ksi_select(fn_ksi,dist_ma,ksi,signchange)


#            def write_xyz_diffvector(filename,pos,center,unitcell):
#                from saman.ringfunc import shortest_vector
#                ntime = pos.shape[0]
#                natom = pos.shape[1]
#
#                delta = np.zeros(pos.shape)
#                delta2 = np.zeros(pos.shape)
#                for i in range(ntime):
#                    for j in range(natom):
#                        delta[i,j,:] = pos[i,j,:]-center[i,:]
#                delta2 = shortest_vector(delta,unitcell)
#
#                f = file(filename,"w+")
#                for i in range(ntime):
#                    print >> f, natom*2
#                    print >> f,"title"
#                    for j in range(natom):
#                        print >> f, "C"," ".join([str(a) for a in center[i,:]+delta[i,j,:]])
#                    for j in range(natom):
#                        print >> f, "C8"," ".join([str(a) for a in center[i,:]+delta2[i,j,:]])
#                f.close()
#                print "file written...",filename

            #filename = fn_ksi+".traj.xyz"
            #write_xyz_diffvector(filename,pos,centers[:,i,:],unitcell0)


