#!/usr/bin/env python


import numpy as np
import transformations as trans


######### Functions ELLIPS #########

def align_xy(pos):
    #align a set of ring coordinates to the xy-plane, with the center on the origin
    #1)shift the coordinates to to have the center on the origin
    #2)rotate the coordinates so that the plane is in the xy-plane
  
    natom=len(pos) #number of atoms in the ring
    #shift coordinates so that the center is on the origin:
    center=np.average(pos,axis=0)
    pos=pos-center
    #calculate the mean plane:
    R1=sum(pos[i,:]*np.sin(2*np.pi*i/natom) for i in xrange(natom))
    R2=sum(pos[i,:]*np.cos(2*np.pi*i/natom) for i in xrange(natom))
    plane=np.cross(R1,R2) #plane is the vector normal to the plane
    plane=plane/np.linalg.norm(plane) #normalize
  #   print "normalized plane",plane
    zaxis=np.array([0,0,1])
    xaxis=np.array([1,0,0])
    theta=0.0
    #rotate coordinates so that the mean plane is in the xy-plane:
    if abs(np.dot(plane,zaxis)) != 1.0: #ring is not in the xy-plane
      #calculate the rotation axis:
      #rotate the ring coordinates according to the rotation matrix:
      theta=trans.angle_between_vectors(plane, zaxis)
      angle=theta
      if theta > np.pi/2:
        angle=theta-np.pi
      if theta < -np.pi/2:
        angle=np.pi-theta
  #     print "angle,theta",angle,theta
  #     M = trans.rotation_matrix(-theta, trans.vector_product(plane, zaxis))
      M = trans.rotation_matrix(angle, trans.vector_product(plane, zaxis))
  #     print "M",M[:3,:3]
      pos=np.dot(pos, M[:3,:3].T)
  #     for i in xrange(natom):
  #       pos[i]=dot(M[:3,:3],pos[i])
  #       pos[i]=dot(M,pos[i])
    return pos,theta,plane,center


def fitEllipse(x,y):
    try:
        x = x[:,np.newaxis]
        y = y[:,np.newaxis]
        D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
        S = np.dot(D.T,D)
        C = np.zeros([6,6])
        C[0,2] = C[2,0] = 2; C[1,1] = -1
#        print "S",S
#        print "C",C
        E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
#        print "np.linalg.inv(S)", np.linalg.inv(S)
#        print "np.dot(np.linalg.inv(S), C)",  np.dot(np.linalg.inv(S), C)
#        print "E",E
#        print "V",V
        n = np.argmax(np.abs(E))
        a = V[:,n]
        return a
    except:
        #import sys
        #sys.stderr.write("ellipse calculation error")
        return np.zeros(6,float)

def fitEllipse3D(pos,align=True):
    if align==True:
        pos,theta,plane,center=align_xy(pos)
    x = pos[:,0]
    y = pos[:,1]
    try:
        a = fitEllipse(x,y)
    except:
        a = np.arange(1.,7.)/np.linalg.norm(np.arange(1.,7.))
    return a,pos,theta,plane,center
        #import sys
        #sys.stderr.write("ellipse calculation error")

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(abs(up/down1)) #np.sqrt(up/down1)
    res2=np.sqrt(abs(up/down2)) #np.sqrt(up/down2)
    return np.array([res1, res2])


####################################
######### RINGS/CHANNELS   #########
####################################

class Ring(object):
    def __init__(self,indices,pos,atomtypes,unitcell=None):
        self.indices = indices
        self.pos = np.take(pos,indices,0)
        self.atomtypes = np.take(atomtypes,indices)
        self.unitcell = unitcell
        # some extra ring properties
        if unitcell is None:
            posalign,theta,plane,center = align_xy(self.pos)
            self.newpos = np.take(pos,indices,0)  # just a copy
        else:
            newpos = reposition_pbc(self.pos,unitcell)
            posalign,theta,plane,center = align_xy(newpos)
            self.newpos = newpos
        self.plane = plane
        self.center = center
        self.theta = theta
        self.posalign = posalign  # aligned positions, 2D
    def set_ellips(self):
        a = fitEllipse(self.posalign[0,:],self.posalign[1,:])
        ax1,ax2 = ellipse_axis_length(a)  # in 2D plane
        phi = ellipse_angle_of_rotation(a)  # in 2D plane
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        self.a = a
        self.b = b
        self.ax1 = ax1
        self.ax2 = ax2
        #print "theta,plane,center",theta,plane,center
    def set_passing(self,):
        self.passing = []
    def analyze_composition(self,ingredients=None):
        if ingredients is None:
            ingredients = set(self.atomtypes)
        ingredients = sorted(ingredients)
        ingred_freq = []
        for ingred in ingredients:
            freq = np.sum(self.atomtypes == ingred)
            ingred_freq.append(freq)
        self.ingred = ingredients
        self.ingred_freq = ingred_freq
        #print "ingredients",ingredients
        #for i,ingred in enumerate(ingredients):
        #    print ingred,ingred_freq[i]
 
    def are_neighbors_outside_self(self,ring1,ring2):
        #check if two rings are neighbors with at least one of their common T-atoms not part of the self ring
        are_neighbors = False
        t_indices = ring1.indices[::2]
        nt_indices = len(t_indices)
        for i in xrange(nt_indices):
            first = t_indices[i]
            second = t_indices[(i+1)%nt_indices] # pair of atoms, they are neighbors in the ring
            if first in ring2.indices and second in ring2.indices and not (first in self.indices and second in self.indices):
                are_neighbors = True
                break
        return are_neighbors
    def compare_neighbors(self,neighbors1,neighbors2):
        are_identical = False
        deque1=deque(neighbors1)
        deque1rev=deque(neighbors1[::-1])
        for i in xrange(len(neighbors1)):
            if list(deque1) == neighbors2 or list(deque1rev) == neighbors2:
                are_identical = True
                break
            deque1.rotate(1)
            deque1rev.rotate(1)
        return are_identical
    def set_connected_neighbors(self,neighbors):
        nneighbors = len(neighbors)
        connected_neighbors = [neighbors[0]]
        neighbors.remove(neighbors[0])
        for i in xrange(nneighbors-1):
            foundaneighbor = False
            for j,ring in enumerate(neighbors):
                if self.are_neighbors_outside_self(connected_neighbors[i],ring):
                    foundaneighbor = True
                    break
            if not foundaneighbor:
                break
            connected_neighbors.append(ring)
            neighbors.remove(ring)
        return neighbors, connected_neighbors
    def set_neighbors(self,rg):
        #self.neighbors1 and self.neighbors2 are two lists of neighboring rings at opposite sides of the self ring
        print "self",self.indices
        neighbors = []
        t_indices = self.indices[::2]   # t_indices are consecutive T-atoms in the ring
        nt_indices = len(t_indices)
        for i in xrange(nt_indices):
            first = t_indices[i]
            second = t_indices[(i+1)%nt_indices] # pair of atoms, they are neighbors in the ring
            for j,ring in enumerate(rg.list_rings):
                if ring.indices == self.indices:
                    continue
                if first in ring.indices and second in ring.indices:
                    #r = ring.copy()
                    neighbors.append(ring)
        neighbors = list(set(neighbors)) #remove duplicate rings
        print "neighbors",[ring.indices for ring in neighbors]
        print "nneighbors",len(neighbors)
        neighbors, self.neighbors1 = self.set_connected_neighbors(neighbors)
        neighbors, self.neighbors2 = self.set_connected_neighbors(neighbors)
        assert len(neighbors) == 0
        print "neighbors1",[ring.indices for ring in self.neighbors1]
        print "neighbors2",[ring.indices for ring in self.neighbors2]
        self.neighbors1_sizes = [len(ring.indices) for ring in self.neighbors1]
        self.neighbors2_sizes = [len(ring.indices) for ring in self.neighbors2]
        print "neighbors1_sizes",self.neighbors1_sizes
        print "neighbors2_sizes",self.neighbors2_sizes
        self.neighbors_identical = self.compare_neighbors(self.neighbors1_sizes,self.neighbors2_sizes)
        print "neighbors_identical",self.neighbors_identical


class RingGroup(object):
    def __init__(self,list_indices,pos,atomtypes,unitcell=None):
        list_rings = []
        for i,indices in enumerate(list_indices):
            #print "ring",i
            ring = Ring(indices,pos,atomtypes,unitcell=unitcell)
            list_rings.append(ring)
#            # too large?
#            center = np.sum(ring.pos,0)/len(ring.indices)
#            dist2 = np.sum((ring.pos-center)**2,1)
#            if np.sum(dist2 > 8.**2) > 0:
#                print "alert",dist2
#            else:
#                list_rings.append(ring)

# for testing:
#            # too large?
#            dist2 = np.sum((ring.newpos-ring.center)**2,1)
#            if np.sum(dist2 > 8.**2) > 0:
#                print "alert",dist2
#            else:
#                list_rings.append(ring)

        self.list_rings = list_rings
        self.list_indices = [ring.indices for ring in list_rings]
        self.list_pos = [ring.pos for ring in list_rings]
        self.list_atomtypes = [ring.atomtypes for ring in list_rings]
        self.nring = len(self.list_rings)
    def set_ellips(self):
        for i in range(self.nring):
            self.list_rings[i].set_ellips()
    def set_passing(self):
        for i in range(self.nring):
            self.list_rings[i].set_passing()
    def analyze_ellips(self):
        # histogram theta
        theta = [ring.theta for ring in self.list_rings]
        import matplotlib.pyplot as plt
        plt.figure()
        plt.hist(theta,bins=50)
        plt.savefig("histogram.theta.png")
        # histogram ax1
        ax1 = [ring.ax1 for ring in self.list_rings]
        ax2 = [ring.ax2 for ring in self.list_rings]
        import matplotlib.pyplot as plt
        plt.figure()
        plt.hist(ax1,bins=100)
        #plt.hist(ax2,bins=100)    
        plt.savefig("histogram.ax.png")
        # correlation
        plt.figure()
        plt.plot(ax1,ax2)
        plt.savefig("fig_corr.ax12.png")
        plt.figure()
        plt.plot(ax1,theta)
        plt.savefig("fig_corr.ax1.theta.png")

    def analyze_composition(self,ingredients=None):
        self.ingred = ingredients
        for i in range(self.nring):
            self.list_rings[i].analyze_composition(ingredients=ingredients)
        self._analyze_composition()

    def _analyze_composition(self,):
        # assert all rings have same ingredients list
        # TODO
        # first ring
        list_comps = [[f for f in self.list_rings[0].ingred_freq]]
        list_rings_percomp = [[0]]
        list_compsrings = [0]
        # loop over other rings
        for i in range(1,self.nring):
            orig = [False for freq in list_comps]
            for j,freq in enumerate(list_comps):
                for k in range(len(freq)):
                    if freq[k] != self.list_rings[i].ingred_freq[k]:
                        orig[j] = True
            orig1 = np.product(orig)
            if orig1:   # everywhere True
                list_comps.append([f for f in self.list_rings[i].ingred_freq])
                list_rings_percomp.append([i])  # add separately
                list_compsrings.append(len(list_comps)-1)
            else:       # at least one False (and should be just one False)
                assert np.sum(orig)==len(orig)-1
                index = orig.index(False)
                list_rings_percomp[index].append(i)  # add to existing composition
                list_compsrings.append(index)
        self.list_comps = list_comps
        self.list_rings_percomp = list_rings_percomp
        self.list_compsrings = list_compsrings
    def print_composition(self,):
        print "-"*5
        print "Composition per ring"
        if self.ingred is None:
            print "ingredients: determined ring by ring"
            for i in range(self.nring):
                print " ".join(str(a) for a in self.list_rings[i].ingred)," ".join(str(a) for a in self.list_rings[i].ingred_freq)          
        else:
            print "Composition per ring"   
            print "ingredients:"," ".join(self.ingred)
            for i in range(self.nring):
                print "ring %i: "%i," ".join(str(a) for a in self.list_rings[i].ingred_freq)
        print "-"*5
        print "Rings ordened per composition"
        if self.ingred is None:
            print "not done as not all rings might have identical ingredient list"
        else:
            print " ".join(self.ingred)
            for i,freq in enumerate(self.list_comps):
                print "composition"," ".join(str(a) for a in freq)," for rings ", " ".join(str(a) for a in self.list_rings_percomp[i])
        print "-"*5
        print "List of ring compositions"
        if self.ingred is None:
            print "not done as not all rings might have identical ingredient list"
        else:
            print self.list_compsrings


class Channel(object):
    def __init__(self,ring1,ring2,):
        self.ring1 = ring1
        self.ring2 = ring2
        d2 = np.sum((ring1.center-ring2.center)**2)
        d = np.sqrt(d2)
        vec = ring1.center-ring2.center
        vec = vec/np.linalg.norm(vec)
        vec = vec*abs(vec[0])/vec[0]
        self.d = d
        self.vec = vec


class ChannelGroup(object):
    def __init__(self,rg):
        self.pairs = []
        self.list_channel = []
        nring = rg.nring
        for i in range(nring):
            for j in range(i+1,nring):
                channel = Channel(rg.list_rings[i],rg.list_rings[j])
                print "i,j,vec",i,j,channel.vec
                if channel.d < 10.:
                    self.list_channel.append(channel)
                    self.pairs.append([i,j])
    def plot_channel(self,fn_hist,fn_corrvec,fn_arrows):
        list_channel = self.list_channel
        import matplotlib.pyplot as plt
        # plot vecs
        plt.figure()
        #vecs = np.array((3,len(list_channel)))
        for i in range(3):
            vecs = [chan.vec[i] for chan in list_channel]
            plt.hist(vecs,bins=100)
        plt.savefig(fn_hist)
        # plot corr vecs
        plt.figure()
        v0 = [chan.vec[0] for chan in list_channel]
        v1 = [chan.vec[1] for chan in list_channel]
        v2 = [chan.vec[2] for chan in list_channel]
        plt.plot(v0,v1,"o")
        plt.plot(v0,v2,"o")
        plt.plot(v1,v2,"o")
        plt.savefig(fn_corrvec)
        # cluster
        import matplotlib.pyplot as plt
        import mayavi.mlab as mlab
        #from mlab import quiver3d
        mlab.options.backend = 'envisage'         # one way to save visualization
        #f = mlab.figure()
        mlab.figure()
        x = [np.zeros(v.shape) for v in v0]
        mlab.quiver3d(x,x,x,v0,v1,v2)
        mlab.savefig(fn_arrows)

    def make_tcl_channel(self,fn_tcl,fn_config,):
        #Define our own version of arrow in vmd
        tcl = r"""proc vmd_draw_arrow {mol start end} {
            # an arrow is made of a cylinder and a cone
            set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
            graphics $mol cylinder $start $middle radius 0.1
            graphics $mol cone $middle $end radius 0.2
        }"""
        #Plot the molecule itself
        #tcl += """\nmol new /home/an/work/germansastre/zeoTsites/results_sam/sapo34/config.xyz\n"""
        tcl += "\nmol new %s\n"%(fn_config)
        #Set the charge of all molecules
        #for i in xrange(sys.natom):
        #    tcl += """\tset sel [atomselect 0 "index %d"]\n\tatomselect%d set charge %.3f\n""" % (i,i,charges[i])
        #Update settings for plotting atoms
        tcl += r"""mol delrep 0 top
        mol representation CPK 0.5 0.3 10 10
        #mol color Charge
        mol selection {all}
        mol material Opaque
        mol addrep top"""
        #Set arrows color
        tcl += """\ndraw color green\n"""
        #Write an arrow for all dipoles
        for i in xrange(len(self.list_channel)):
            start = self.list_channel[i].ring1.center
            end = self.list_channel[i].ring2.center  #+list_channel[i].vec
            tcl += """\tdraw arrow {%.2f %.2f %.2f} {%.2f %.2f %.2f}\n""" % (start[0], start[1], start[2] , end[0], end[1], end[2] )
        #for i in xrange(rg.nring):
        #  for j in xrange(len(rg.list_rings[i].indices)):
        #    start = rg.list_rings[i].pos[j,:]
        #    end = rg.list_rings[i].center  #+np.ones((3),float)*0.2  #+list_channel[i].vec
        #    tcl += """\tdraw arrow {%.2f %.2f %.2f} {%.2f %.2f %.2f}\n""" % (start[0], start[1], start[2] , end[0], end[1], end[2] )
        #Write tcl script to file
        with open(fn_tcl,'w') as f:
            f.write(tcl)


class Passing(object):
    """Object to gather all the information,
    such that transitions can be counted afterwards"""
    def __init__(self,atoms,time,dist,ksi,signchange,shift=0):
        self.atoms = atoms
        self.time = time + shift
        self.dist = dist
        self.ksi = ksi
        self.signchange = signchange
        # extra
        self.transitions = np.sum(abs(signchange))
        self.netto = np.sum(signchange)
    


def make_tcl_ring(ring,fn_tcl,fn_config,):
    # Define our own version of arrow in vmd
    # use as follows:
    #     vmd -e fn_tcl config.xyz
    tcl = r"""proc vmd_draw_arrow {mol start end} {
        # an arrow is made of a cylinder and a cone
        set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
        graphics $mol cylinder $start $middle radius 0.1
        graphics $mol cone $middle $end radius 0.2
    }"""
    #Plot the molecule itself
    #tcl += """\nmol new /home/an/work/germansastre/zeoTsites/results_sam/sapo34/config.xyz\n"""
    tcl += "\nmol new %s\n"%(fn_config)
    #Set the charge of all molecules
    #for i in xrange(sys.natom):
    #    tcl += """\tset sel [atomselect 0 "index %d"]\n\tatomselect%d set charge %.3f\n""" % (i,i,charges[i])
    #Update settings for plotting atoms
    tcl += r"""mol delrep 0 top
    mol representation CPK 0.5 0.3 10 10
    #mol color Charge
    mol selection {all}
    mol material Opaque
    mol addrep top"""
    #Set arrows color
    tcl += """\ndraw color green\n"""
    #Write an arrow for all dipoles
    for i in xrange(len(ring.indices)):
        start = ring.pos[i,:]
        end = ring.pos[i,:] + np.ones(3,float)*0.3
        tcl += """\tdraw arrow {%.2f %.2f %.2f} {%.2f %.2f %.2f}\n""" % (start[0], start[1], start[2] , end[0], end[1], end[2] )
    for i in xrange(len(ring.indices)):
        start = ring.newpos[i,:]
        end = ring.newpos[i,:] + np.ones(3,float)*0.3
        tcl += """\tdraw arrow {%.2f %.2f %.2f} {%.2f %.2f %.2f}\n""" % (start[0], start[1], start[2] , end[0], end[1], end[2] )
    start = ring.center
    end = start + np.ones(3,float)*0.3
    tcl += """\tdraw arrow {%.2f %.2f %.2f} {%.2f %.2f %.2f}\n""" % (start[0], start[1], start[2] , end[0], end[1], end[2] )
    start = ring.center
    end = start + ring.plane
    tcl += """\tdraw arrow {%.2f %.2f %.2f} {%.2f %.2f %.2f}\n""" % (start[0], start[1], start[2] , end[0], end[1], end[2] )
    #for i in xrange(rg.nring):
    #  for j in xrange(len(rg.list_rings[i].indices)):
    #    start = rg.list_rings[i].pos[j,:]
    #    end = rg.list_rings[i].center  #+np.ones((3),float)*0.2  #+list_channel[i].vec
    #    tcl += """\tdraw arrow {%.2f %.2f %.2f} {%.2f %.2f %.2f}\n""" % (start[0], start[1], start[2] , end[0], end[1], end[2] )
    #Write tcl script to file
    with open(fn_tcl,'w') as f:
        f.write(tcl)


########################################
### functions
########################################
def passing_average(list_passing,average="runs"):
    nfile = len(list_passing)
    tottime = np.sum([len(pas.time) for pas in list_passing])
    #print "nfile",nfile
    #print "tottime",tottime
    allatoms = np.zeros(tottime,int)
    alltime = np.zeros(tottime,int)
    alldist = np.zeros(tottime,float)
    allksi = np.zeros(tottime,float)
    allsignchange = np.zeros(tottime,int)
    if average == "oneheap":
        i = 0
        for n in range(nfile):
            pas = list_passing[n]
            ntime = len(pas.time)
            allatoms[i:i+ntime] = pas.atoms[:]
            alltime[i:i+ntime] = pas.time[:]
            alldist[i:i+ntime] = pas.dist[:]
            allksi[i:i+ntime] = pas.ksi[:]
            allsignchange[i:i+ntime] = pas.signchange[:]
            i += ntime
        assert i == tottime
        return Passing(allatoms,alltime,alldist,allksi,allsignchange)

def reposition_pbc(pos,unitcell):
    """Translate atoms such that ring is coherent, not cut by pbc.
       The function returns new positions that are
       adapted to form a coherent ring. In practice, a ring is gradually
       built up by adding each of the ring atoms one by one, and
       by choosing the image of the atom that is closest to the
       center of mass of the ring. The center of mass is updated.
    """
    newpos = np.zeros(pos.shape,float)
    newpos[:] = pos[:]   # copy coordinates
    natom = len(pos)
    notused = range(natom)
    used = []
    #print "used",used
    #print "notused",notused

    masses = np.ones(natom)   # equal weight for each atom
    com = np.dot(masses,pos)/np.sum(masses)
    for i in range(natom):   # repeat untill all atoms have been used
        pos_diff = np.take(newpos,notused,0) - com  # distance to com
        reciproc = np.linalg.inv(unitcell.transpose())  # (unitcell^T)^(-1)
        direct = np.dot(pos_diff,reciproc)
        direct_minimage = direct - np.round(direct)   # minimum image
        pos_diff_minimage = np.dot(direct_minimage,unitcell.transpose())

        # distance
        dist2 = np.sum(pos_diff_minimage**2,1)
        #print "dist2",dist2
        closest = np.argsort(dist2)
        #print "closest",closest
        closest = closest[0]
        closest_at = notused[closest]
        #print "closest,closest_at",closest,closest_at

        # update
        used.append(closest_at)
        notused.remove(closest_at)
        #print "used",used
        #print "notused",notused

        # adapt newpos, com
        newpos[closest_at,:] = pos_diff_minimage[closest,:]+com   # that shortest vector piece
        pos_used = np.take(newpos,used,0)
        masses_used = np.take(masses,used,0)
        com = np.dot(masses_used,pos_used)/np.sum(masses_used)
    return newpos

def get_centers_planes(allcoorframework,rg,atomtypes,unitcell=None):
    # create time series of center of rings 
    # allcoorframework  --  ntime x natom x 3
    # rg  --  RingGroup instance
    # atomtypes  --  atomtypes of all frameworkatoms
    ntime = len(allcoorframework)
    centers = np.zeros((ntime,rg.nring,3))
    planes = np.zeros((ntime,rg.nring,3))
    #print "ntime,natom,3 (allcoorframework)",allcoorframework.shape
    #print "ntime,nring,3 (centers)",centers.shape
    for t in range(ntime):
        for i in range(rg.nring):
            ring = Ring(rg.list_indices[i],allcoorframework[t,:,:],atomtypes,unitcell=unitcell)
            centers[t,i,:] = ring.center[:]
            planes[t,i,:] = ring.plane[:]
    return centers,planes

#without ring object
def get_centers_planes_bis(allcoorframework,list_ring_indices,atomtypes,unitcell=None):
    # create time series of center of rings 
    # allcoorframework  --  ntime x natom x 3
    # list_ring_indices  --  list of ring indices lists
    # atomtypes  --  atomtypes of all frameworkatoms
    ntime = len(allcoorframework)
    nring = len(list_ring_indices)
    centers = np.zeros((ntime,nring,3))
    planes = np.zeros((ntime,nring,3))
    #print "ntime,natom,3 (allcoorframework)",allcoorframework.shape
    #print "ntime,nring,3 (centers)",centers.shape
    for t in range(ntime):
        for i in range(nring):
            ring = Ring(list_ring_indices[i],allcoorframework[t,:,:],atomtypes,unitcell=unitcell)
            centers[t,i,:] = ring.center[:]
            planes[t,i,:] = ring.plane[:]
    return centers,planes

def shortest_vector(delta,unitcell):
    reciprocal = np.linalg.inv(unitcell).T
    direct = np.dot(delta, reciprocal)
    direct_minimage = direct - np.round(direct)    # minimum image
    return np.dot(direct_minimage,unitcell.T)  # change cartesian
    #direct_minimage = np.floor(direct + 0.5)  # fractional
    #return delta - np.dot(direct_minimage, unitcell.T)  # change cartesian

def check_location(ring,pos,unitcell):
    delta = pos-ring.center
    shortest = shortest_vector(delta,unitcell)
    dist2 = np.sum(shortest**2,axis=-1)
    dist = np.sqrt(dist2)
    ksi = np.dot(shortest,ring.plane)
    h = np.sqrt(dist2-ksi**2)
    #print "ring.plane",ring.plane.shape
    #print "pos",pos.shape
    #print "dist",dist.shape
    #print "ksi",ksi.shape
    return dist,ksi,h

def check_location_bis(center,plane,pos,unitcell):
    delta = pos-center
    shortest = shortest_vector(delta,unitcell)
    dist2 = np.sum(shortest**2,axis=-1)
    dist = np.sqrt(dist2)
    #ksi = np.dot(shortest,plane)
    ksi = np.sum(shortest*plane,axis=-1)
    #h = np.sqrt(dist2-ksi**2)
    #print "in check_location_bis"
    #print "shortest",shortest.shape
    #print "ring.plane",ring.plane.shape
    #print "pos",pos.shape
    #print "dist",dist.shape
    #print "ksi",ksi.shape
    return dist,ksi

# code from MOLMOD package 
if False:
    def to_fractional(self, cartesian):
        """Convert Cartesian to fractional coordinates

           Argument:
            | ``cartesian``  --  Can be a numpy array with shape (3, ) or with shape
                                 (N, 3).

           The return value has the same shape as the argument. This function is
           the inverse of to_cartesian.
        """
        return numpy.dot(cartesian, self.reciprocal)

    def to_cartesian(self, fractional):
        """Converts fractional to Cartesian coordinates

           Argument:
            | ``fractional``  --  Can be a numpy array with shape (3, ) or with shape
                                  (N, 3).

           The return value has the same shape as the argument. This function is
           the inverse of to_fractional.
        """
        return numpy.dot(fractional, self.matrix.transpose())

    def shortest_vector(self, delta):
        """Compute the relative vector under periodic boundary conditions.

           Argument:
            | ``delta``  --  the relative vector between two points

           The return value is not necessarily the shortest possible vector,
           but instead is the vector with fractional coordinates in the range
           [-0.5,0.5[. This is most of the times the shortest vector between
           the two points, but not always. (See commented test.) It is always
           the shortest vector for orthorombic cells.
        """
        fractional = self.to_fractional(delta)
        fractional = numpy.floor(fractional + 0.5)
        return delta - self.to_cartesian(fractional)

    def reciprocal(self):
        """The reciprocal of the unit cell

           In case of a three-dimensional periodic system, this is trivially the
           transpose of the inverse of the cell matrix. This means that each
           column of the matrix corresponds to a reciprocal cell vector. In case
           of lower-dimensional periodicity, the inactive columns are zero, and
           the active columns span the same sub space as the original cell
           vectors.
        """
        U, S, Vt = numpy.linalg.svd(self.matrix*self.active)
        Sinv = numpy.zeros(S.shape, float)
        for i in xrange(3):
            if abs(S[i]) < self.eps:
                Sinv[i] = 0.0
            else:
                Sinv[i] = 1.0/S[i]
        return numpy.dot(U*Sinv, Vt)*self.active


########################################
### Functions: reading
########################################

# reading functions for files generated by ZeoTSites/Sam Moors

def read_rings(fn_rings):
    print "Reading file...", fn_rings
    # read rings
    list_rings_indices = []
    f = file(fn_rings)
    for line in f:
        words = line.split()
        list_rings_indices.append([int(word) for word in words])
    f.close()
    #print "rings"
    #print list_rings_indices
    return list_rings_indices

def read_coords(fn_xyz):
    print "Reading file...", fn_xyz
    # read coordinates
    f = file(fn_xyz)
    for line in f:
        natom = int(line)
        break
    f.next()
    pos = np.zeros((natom,3),float)
    atomtypes = []
    i = 0
    for line in f:
        words = line.split()
        assert len(words)==4
        pos[i,:] = np.array([float(word) for word in words[1:]])
        atomtypes.append(words[0])
        i+=1
    assert len(atomtypes)==natom
    #print "natom",natom
    f.close()
    return natom,atomtypes,pos

def read_unitcell(fn_unitcell):
    print "Reading file...", fn_unitcell
    # read unitcell vectors, one vector per line
    unitcell = np.zeros((3,3),float)  # one vector per column
    i = 0
    f = file(fn_unitcell)
    for line in f:
        words = line.split()
        assert len(words)==3
        unitcell[:,i] = [float(word) for word in words]
        i+=1
    f.close()
    assert i==3
    return unitcell

########################################
# reading functions for data-files, own format

def read_dist(filename):
    data = []
    f = file(filename)
    for line in f:
        if not line.startswith("#"):
            words = line.split()
            data.append([float(word) for word in words])
    f.close()
    return np.array(data)

def read_ksi_truncated(filename):
    atoms = []
    time = []
    dist = []
    ksi = []
    signchange = []
    f = file(filename)
    for line in f:
        if line.startswith("at") or line.startswith("#at"):
            words = line.split()
            at = int(words[1])
        else:
            atoms.append(at)
            words = line.split()
            time.append(int(words[0]))
            dist.append(float(words[1]))
            ksi.append(float(words[2]))
            signchange.append(int(words[3]))
    f.close()
    return np.array(atoms),np.array(time),np.array(dist),np.array(ksi),np.array(signchange)


########################################
### Functions: writing
########################################

def write_ksi(fn_ksi,ksi):
    f = file(fn_ksi, "w+")
    nstep = len(ksi)
    #print >> f, "#---doing ring%i"%i
    for t in range(nstep):
        line = " ".join(map(str,ksi[t,:].tolist()))
        print >> f, line
    f.close()

def write_ksi_select(fn_ksi,dist_ma,ksi,changedsign):
    # ksi -- nstep x at
    # assume masked
    f = file(fn_ksi, "w+")

    nstep,natom = ksi.shape
    for at in range(natom):
        print >> f, "#at",at
        #a = ksi[:,at][~dist.mask[:,at]]
        #b = dist[:,at][~dist.mask[:,at]]
        #c = changedsign[:,at][~dist.mask[:,at]]
        times = np.arange(nstep)[~dist_ma.mask[:,at]]
        for i,t in enumerate(times[:-1]):
            print >> f, "%d  %10.4f  %10.4f %i" %(t,dist_ma[t,at],ksi[t,at],changedsign[t,at])
    f.close()

