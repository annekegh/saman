"""Script to convert HISTORY file
AG, February 25, 2013
AG, April 25, 2013

script to determine which atoms are adsorbates in the zeolite simulations
AG, June 24, 2013

correction: do unwrapping
AG, Sept 10, 2013

split up in file with functions and file with options/args
AG, March 19, 2014

DATA
time data from HISTORY
unitcell from h5
trajectory from h5
atomtypes from HISTORY
molecule sizes from FIELD
"""

import numpy as np
import h5py as h5
from molmod import electronvolt, picosecond,angstrom
from molmod.io import DLPolyHistoryReader  #, DLPolyOutputReader

# for matplotlib figure generation...
import matplotlib
matplotlib.use('Agg')

def read_timedata_from_history(filename_history):
    """Read data from DL_POLY
    time data in p.s."""
    # read HISTORY file
    g = DLPolyHistoryReader(filename_history)
    print "Opening HISTORY...",filename_history
    
    # extract some data
    for r in g:
        dt = r['timestep']   # in atomic units
        t0 = r['time']
        #print r.keys()
        break
    for r in g:
        t1 = r['time']
        break
    dtc = (t1-t0)
    #g.close()
    # convert to ps
    t0 /= picosecond
    t1 /= picosecond
    dt /= picosecond
    dtc /= picosecond

    if True:   # print some data
        print "t0:",t0,"ps"
        print "t1:",t1,"ps"
        print "dt",dt,"ps   integration step"
        print "dtc:",dtc,"ps   time between sampled coords"

    # print np.array(r['symbols'])[ethylene[0]]
    print "Closing HISTORY...",filename_history
    return t0,t1,dt,dtc  #in ps

def read_unitcell_from_h5(filename_h5):
    """read data from h5py file
    unit cell data in angstrom
    format:
      has dimension ntime x 3 x 3
      in h5, unit cell vectors are in rows
      here, they are transposed such that each column is a unit cell vector
      such that it matches the format used by molmod package"""
    # read h5py file
    f = h5.File(filename_h5, mode='r')
    #print f['trajectory'].keys()  # cell, pos, ...
    uc = f['trajectory']['cell']
    #print f['trajectory/cell']
    #print f['trajectory/cell'].shape
    #print "unitcell shape",uc.shape
    #print uc[:10,:,:]
    unitcell = np.zeros(uc.shape,float)
    unitcell[:] = uc[:]/angstrom    # in h5 format
    f.close()

    ##### TODO  WARNING use transpose
    print "WARNING: transposing the unit cell"
    for i in range(len(unitcell)):
        unitcell[i,:,:] = unitcell[i,:,:].transpose()
    # such that unitcell[i,:,j] = unit cell vector j in timestep i
    # this is molmod format
    #####

    return unitcell  # in angstrom

def read_trajectory_from_list_h5(list_filename_h5,selected,combine):
    """
    arguments::
      selected -- list of molecules, each molecule is a list of atom indices,
                  OR None (all are selected) TODO
    """
    assert combine in ["permolecule","perfile","cumul"]
    #list_x = []
    #list_y = []
    #list_z = []
    list_allcoor = []
    for n,filename_h5 in enumerate(list_filename_h5):
        allcoor = read_trajectory_from_h5(filename_h5,selected)
        if combine == "perfile":
            #list_x.append(x)
            #list_y.append(y)
            #list_z.append(z)
            list_allcoor.append(allcoor)
        elif combine == "permolecule":
            for at in range(x.shape[1]):  #natom
                #list_x.append(np.reshape(x[:,at],(-1,1)))
                #list_y.append(np.reshape(y[:,at],(-1,1)))
                #list_z.append(np.reshape(z[:,at],(-1,1)))
                #ist_allcoor.append(np.reshape(allcoor[:,at,:],(-1,3)))
                list_allcoor.append(allcoor[:,at,:])
        elif combine == "cumul":
            raise ValueError("not implemented!")
    print "list_allcoor:",len(list_allcoor)
    print "Closing h5 file...",filename_h5
    #return list_x,list_y,list_z
    return list_allcoor

def read_trajectory_from_h5(filename_h5,adsorbates):
    """Read data h5py file

    trajectory data in angstrom

    arguments::
      adsorbates -- list of molecules, each molecule is a list of atom indices,
                  OR None (all are selected) TODO
    """
    # read h5py file
    f = h5.File(filename_h5, mode='r')
    #print f['trajectory'].keys()  # cell, pos, ...
    print "Opening h5 file...",filename_h5
    print "h5 - trajectory nstep x natom x 3",f['trajectory']['pos'].shape  # nstep x natom x 3

    nadsorbates = len(adsorbates)
    size = len(adsorbates[0])
    natom = nadsorbates*size
    ntime = f['trajectory']['pos'].shape[0]

    # gather the coordinates
    allcoor = np.zeros((ntime,natom,3),float)
    print "size,natom,ntime,nadsorbates",size,natom,ntime,nadsorbates
    print "adsorbates",adsorbates
    for m in range(nadsorbates):
        allcoor[:,m*size:(m+1)*size,:] = f['trajectory']['pos'][:,adsorbates[m],:]/angstrom

    #x = allcoor[:,:,0]
    #y = allcoor[:,:,1]
    #z = allcoor[:,:,2]
    f.close()
    print "shape allcoor:",allcoor.shape, "ntime x natom x 3"
    return allcoor  # in angstrom

def read_info(f,):
    nummols = 0
    numatoms = 0
    for line in f:
        if line.startswith("NUMMOLS"):
            nummols = int(line.split()[1])
            break
    for line in f:
        if line.startswith("ATOMS"):
            numatoms = int(line.split()[1])
            break
    return nummols, numatoms

def read_atomtypes_from_history(filename_history):
    f = file(filename_history)
    atomtypes = []
    for line in f:
        if line.strip().startswith("timestep"): break
    for line in f:
        if line.strip().startswith("timestep"): break
        words = line.split()
        if len(words) == 4:
            atomtypes.append(words[0])
    f.close()
    return atomtypes

def read_numbers_fieldfile(fieldfile):
    """read information from FIELD file

    WARNING
      assuming the following order in the FIELD file:
      1. framework (one type),
      2. adsorbate molecules (one type),
      3. product molecules (one type).
      The product molecules do not need to be present.

    nframework -- number of atoms in framework
    size -- number of atoms in one adsorbate molecule
    nadsorbates -- number of adsorbate molecules
    """
    f = file(fieldfile,"r")
    
    # framework: zeo or sapo
    nummols,numatoms = read_info(f)
    assert nummols == 1
    nframework = numatoms
    # adsorbates: ethylene, propene
    nummols,numatoms = read_info(f)
    nadsorbates = nummols
    size = numatoms
    assert (size==6 or size==9)
    # products: hmb, pyr
    nproducts = 0
    try:
       nummols,numatoms = read_info(f)
       nproducts = nummols
    except EOFError:
       pass
    
    f.close()

    if True:
        print "nframework",nframework
        print "nadsorbates",nadsorbates
        print "size",size
        print "nproducts",nproducts

    return nframework,nadsorbates,size,nproducts

def convert_t2n(t,dtc):
    # t = n*dtc  # time in ps
    return t/float(dtc)

######################################################################

def get_times_MSD(times,filename_history):
    # MSD settings for fit: timings [ps]
    t_start = times[0]
    t_end   = times[1]
    t_delta = times[2]

    # time information
    t0,t1,dt,dtc = read_timedata_from_history(filename_history)
    
    # standard MSD
    dn1 = int(np.ceil(convert_t2n(t_start,dtc)))
    dn2 = int(np.floor(convert_t2n(t_end,dtc)))
    ddn = int(np.round(convert_t2n(t_delta,dtc)))
    print "in fit: (dn1,dn2,ddn)",dn1,dn2,ddn, "sample"
    print "in fit: (t_start,t_end,t_delta)",dn1*dtc,dn2*dtc,ddn*dtc, "ps"
    print "asked: (t_start,t_end,t_delta)",t_start,t_end,t_delta,"ps"
    return t_start,t_end,t_delta,t0,t1,dt,dtc,dn1,dn2,ddn


def get_seladsorbates(filename_field,filename_history,selected_atomtypes,):
    """Create list with adsorbate indices

    Arguments::
      filename_field  --  FIELD file of DLPOLY
      filename_history  --  HISTORY file of DLPOLY
      selected_atomtypes  --  list with atomtypes to be selected,
                            e.g. ['C2','H2']

    Returns: list consisting of atom index lists, one list for
    each adsorbate, e.g. [[0,1],[10,15],[12,9]] for 3 adsorbates
    """

    # Read: All atomtypes
    atomtypes = read_atomtypes_from_history(filename_history)

    # Read: select the ethylene or propene molecules - numbering starts with 0
    nframework,nadsorbates,size,nproducts = read_numbers_fieldfile(filename_field)
    start = nframework                      # index of first adsorbate atom
    end = nframework + size*nadsorbates -1  # index of last adsorbate atom
    select = range(start,end+1)
    natom = len(select)
    assert nadsorbates == natom/size
    print "natom,nadsorbates,size", natom,nadsorbates,size
    # All adsorbates
    adsorbates = [range(start+size*i,start+size*(i+1)) for i in range(nadsorbates)]

    # Select some adsorbate atoms
    print "selected_atomtypes",selected_atomtypes
    seladsorbates = [ [k for k in range(start+size*i,start+size*(i+1)) if atomtypes[k] in selected_atomtypes] for i in range(nadsorbates)]
    # only the first atom:   TODO   # like German does
    #seladsorbates = [ [a[0]] for a in seladsorbates ]
    # framework atoms:
    #seladsorbates = [ [k for k in range(nframework) if atomtypes[k] in selected_atomtypes] ]

    if "sapo56_32ethylene_300K" in filename_history:
        print "WARNING - manually delete one ethene from specific run"
        seladsorbates = [ a for a in seladsorbates if 1240 not in a and 1243 not in a]
        print "Deleted atoms from sapo56_32ethylene_300K: 1240, 1243 (numbering starts with 0)"

    print "seladsorbates",seladsorbates   # list of lists
    return seladsorbates

def get_selframework(filename_field,filename_history,selected_atomtypes,):
    """Create list with framework indices

    Arguments::
      filename_field  --  FIELD file of DLPOLY
      filename_history  --  HISTORY file of DLPOLY
      selected_atomtypes  --  list with atomtypes to be selected,
                            e.g. ['Si','Al']

    Returns: list consisting of atom index lists, one list for
    each adsorbate, e.g. [[0,1],[10,15],[12,9]] for 3 adsorbates
    """

    # Read: All atomtypes
    atomtypes = read_atomtypes_from_history(filename_history)

    # Read: select framework atoms - numbering starts with 0
    nframework,nadsorbates,size,nproducts = read_numbers_fieldfile(filename_field)
    print "nframework",nframework
    framework = range(nframework)

    # Select some framework atoms
    print "selected_atomtypes",selected_atomtypes
    selframework = [ k for k in range(nframework) if atomtypes[k] in selected_atomtypes]
    print "selframework",selframework
    return selframework  #list of atom indices

def get_xyz(filename_field,list_filename_h5,filename_history,combine,selected_atomtypes,atoms="adsorbates"):
    """
    arguments::
      selected_atomtypes -- list of atomtypes, e.g. ['C2','H2']
      atoms -- string to decide whether to select framework atoms or adsorbate atoms

    returns::
      list_allcoor -- list of ntime x natom x 3 arrays, one array per run or per molecule
      unitcell -- ntime x 3 x 3
      selected -- list of molecules, each molecule is a list of atom indices,
                  OR None (all are selected) TODO
    """

    # which atoms to take
    if atoms == "adsorbates":
        selected = get_seladsorbates(filename_field,filename_history,selected_atomtypes,)
    elif atoms == "framework":
        selected = [get_selframework(filename_field,filename_history,selected_atomtypes,)]
    else:
        selected = None
    # unit cell ntime x 3 x 3 (in molmod format: cell vectors in columns)
    unitcell = read_unitcell_from_h5(list_filename_h5[0])
    print "unitcell",unitcell[0,:,:]
    # trajectory
    list_allcoor = read_trajectory_from_list_h5(list_filename_h5,selected,combine)

    return list_allcoor,unitcell,selected 


######################################################################

