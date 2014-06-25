#!/usr/bin/env python


import numpy as np
import transformations as trans
from collections import deque
from math import atan2
from ring_sasa import ring_sasa


def ringpucker(a1,align=False):
    #see: A General Definition of Ring Puckering Coordinates, Kremer and Pople, JACS 1975
    #a1 is a vector of atom coordinates
    natom=len(a1) #number of atoms in the ring
    assert natom >= 4
    #fit coordinates to xy-plane with geometric center on the origin
    if align == True:
        a1,theta=align_xy(a1)
    #select the z-component of the aligned coordinates:
    z=a1[:,2]
    if natom%2==1:
        nm=(natom-1)/2 #number of possible m values
        q_phi=np.zeros([nm,2])
    else:
        nm=natom/2
        q_phi=np.zeros([nm+1,2])
    q_cos_phi=np.zeros(nm)
    q_sin_phi=np.zeros(nm)
    norm=np.sqrt(2.0/natom)
    for m in np.arange(2,nm):
        for j in np.arange(natom):
            q_cos_phi[m]=q_cos_phi[m]+z[j]*np.cos(2*np.pi*m*j/natom)
            q_sin_phi[m]=q_sin_phi[m]+z[j]*np.sin(2*np.pi*m*j/natom)
        q_cos_phi[m]=q_cos_phi[m]*norm #normalization
        q_sin_phi[m]=q_sin_phi[m]*(-1)*norm #normalization
        q_phi[m,1]=atan2(q_sin_phi[m],q_cos_phi[m])
        q_phi[m,0]=q_cos_phi[m]/np.cos(q_phi[m,1])
    if natom%2==0:
        m=natom/2
        qn2=0
        for j in np.arange(natom):
            qn2=qn2+(-1)**j*z[j]
        qn2=qn2*norm/2
        q_phi[m,0]=qn2
    q_phi[0,0]=sum(z*z) #total puckering amplitude
    return q_phi

######### Functions ELLIPS #########

def align_xy(pos):
    #align a set of ring coordinates to the xy-plane, with the center on the origin
    #1)shift the coordinates to to have the center on the origin
    #2)rotate the coordinates so that the plane is in the xy-plane
  
    natom=len(pos) #number of atoms in the ring
    #shift coordinates so that the center is on the origin:
    center=np.average(pos,axis=0)
    pos0=pos-center
    #calculate the mean plane:
    R1=sum(pos0[i,:]*np.sin(2*np.pi*i/natom) for i in xrange(natom))
    R2=sum(pos0[i,:]*np.cos(2*np.pi*i/natom) for i in xrange(natom))
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
      pos0=np.dot(pos0, M[:3,:3].T)
  #     for i in xrange(natom):
  #       pos[i]=dot(M[:3,:3],pos[i])
  #       pos[i]=dot(M,pos[i])
    return pos0,theta,plane,center


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
        pos_align,theta,plane,center=align_xy(pos)
    x = pos_align[:,0]
    y = pos_align[:,1]
    try:
        a = fitEllipse(x,y)
    except:
        a = np.arange(1.,7.)/np.linalg.norm(np.arange(1.,7.))
    return a,pos_align,theta,plane,center
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
    def __init__(self,indices,pos,atomtypes,unitcell=None,maxradius=6.,orient=None):
        self.indices = indices
        self.pos = np.take(pos,indices,0)    # natom-in-ring x 3
        self.atomtypes = np.take(atomtypes,indices)
        self.unitcell = unitcell
        # some extra ring properties
        if unitcell is None:
            self.newpos = np.take(pos,indices,0)  # just a copy of self.pos
        else:
            # too large?
            center = np.sum(self.pos,0)/float(len(self.indices))
            dist2 = np.sum((self.pos-center)**2,1)
            if np.sum(dist2 > maxradius**2) > 0:  # if any of the atoms is too far from center
                self.newpos = reposition_pbc(self.pos,unitcell)  # make ring coherent
            else:
                self.newpos = np.take(pos,indices,0)  # just a copy of self.pos
        posalign,theta,plane,center = align_xy(self.newpos)
        self.posalign = posalign  # aligned positions, more or less in xy-plane
        self.theta = theta        # angle of plane with ...
        self.plane = plane        # normal vector on ring plane
        self.center = center      # center of the ring
        # adapt vector orient to largest direction or to given direction
        if orient is None:
            index = self.get_orient()
        else: index = orient
        if self.plane[index] < 0:
            self.plane *= -1.
        self.orient = index
        # print in xyz format:
        #print len(posalign)
        #print "posalign"
        #for i in range(len(posalign)): print "C %8.3f %8.3f %8.3f"%(posalign[i,0],posalign[i,1],posalign[i,2])
    def get_orient(self):
        args = np.argsort(abs(self.plane))   # the largest component
        return args[-1]
    def set_ellips(self):
        a = fitEllipse(self.posalign[:,0],self.posalign[:,1])
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
    def set_dmin(self):
        #calculates the minimum diagonal O-O distance
        no_atom = len(self.indices)/2
        o_pos1 = self.posalign[1:no_atom:2]
        o_pos2 = self.posalign[no_atom+1::2]
        difvecs = o_pos1 - o_pos2
        difvecs = difvecs.T
        self.diameters = np.sqrt(difvecs[0]**2 + difvecs[1]**2 + difvecs[2]**2)
#         print "diameters",self.diameters
        self.dmin = np.min(self.diameters)
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
        from collections import deque
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
    def set_neighbors(self,rg,doprint=False):
        #self.neighbors1 and self.neighbors2 are two lists of neighboring rings at opposite sides of the self ring
        print "setting neighbors"
        if doprint: print "self",self.indices
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
        if doprint:
            print "neighbors",[ring.indices for ring in neighbors]
            print "nneighbors",len(neighbors)
        neighbors, self.neighbors1 = self.set_connected_neighbors(neighbors)
        neighbors, self.neighbors2 = self.set_connected_neighbors(neighbors)
        assert len(neighbors) == 0
        self.neighbors1_sizes = [len(ring.indices) for ring in self.neighbors1]
        self.neighbors2_sizes = [len(ring.indices) for ring in self.neighbors2]
        self.neighbors_identical = self.compare_neighbors(self.neighbors1_sizes,self.neighbors2_sizes)
        if doprint:
            print "neighbors1",[ring.indices for ring in self.neighbors1]
            print "neighbors2",[ring.indices for ring in self.neighbors2]
            print "neighbors1_sizes",self.neighbors1_sizes
            print "neighbors2_sizes",self.neighbors2_sizes
            print "neighbors_identical",self.neighbors_identical
    def set_ringpucker(self):
        self.q_phi = ringpucker(self.posalign)
    def set_radii(self, dictradii):
        #dictradii is a dictionary, for example: {'Si': 2.95, 'Al': 2.86, 'P': 2.97, 'O2': 2.66, 'O1': 2.81}
        self.radii = np.array([dictradii[i] for i in self.atomtypes])
    def set_sasa(self, solventradius):
        #calculate the solvent accessible area of the ring window,
        #based on the projected 2D surface of the ring atoms
        self.sasapolygon,self.sasa = ring_sasa(self.posalign[:,:2], self.radii, solventradius)


class RingGroup(object):
    def __init__(self,list_indices,pos,atomtypes,unitcell=None):
        list_rings = []
        for i,indices in enumerate(list_indices):
            #print "ring",i
            ring = Ring(indices,pos,atomtypes,unitcell=unitcell,orient=None)
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
        self.list_orient = [ring.orient for ring in list_rings]
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

        # set compositon of individual rings
        for i in range(self.nring):
            self.list_rings[i].analyze_composition(ingredients=ingredients)

        # analyze composition of all rings
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
        signchange3 = []
        c = 0
        for i in range(len(signchange)):
            if c == len(signchange)-1:
                signchange3.append(0)
                c += 1
            elif signchange[c]==0:
                signchange3.append(0)
                c += 1
            elif signchange[c]==1:
                if signchange[c+1]==-1:
                    signchange3.append(0)
                    signchange3.append(0)
                    c += 2
                elif signchange[c+1]==0:
                    signchange3.append(1)
                    signchange3.append(0)
                    c += 2
                elif signchange[c+1]==1:
                    signchange3.append(1) #raise ValueError
                    c += 1
            elif signchange[c]==-1:
                #if signchange[c+1]==-1:
                #    print "time",c
                #    print "ksi(time-1,time,time+1,time+2)",ksi[c-1:c+3]
                #    print "dist(time-1,time,time+1,time+2)",dist[c-1:c+3]
                #    print "signchange(time-1,time,time+1,time+2)",signchange[c-1:c+3]
                #    #print signchange[c],ksi[c]
                #    #print signchange[c+1],ksi[c+1]
                #    raise ValueError
                if signchange[c+1]==1:
                    signchange3.append(0)
                    signchange3.append(0)
                    c += 2
                elif signchange[c+1]==0:
                    signchange3.append(-1)
                    signchange3.append(0)
                    c += 2
                elif signchange[c+1]==-1:
                    signchange3.append(-1) #raise ValueError
                    c += 1
            if c >= len(signchange): break

        signchange3 = np.array(signchange3)
        self.transitions3 = np.sum(abs(signchange3))
        self.netto3 = np.sum(signchange3)


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
            ring = Ring(rg.list_indices[i],allcoorframework[t,:,:],atomtypes,unitcell=unitcell,orient=rg.list_orient[i])
            centers[t,i,:] = ring.center[:]
            planes[t,i,:] = ring.plane[:]
    return centers,planes

#without ring object
#def get_centers_planes_bis(allcoorframework,list_ring_indices,atomtypes,unitcell=None):
#    # create time series of center of rings 
#    # allcoorframework  --  ntime x natom x 3
#    # list_ring_indices  --  list of ring indices lists
#    # atomtypes  --  atomtypes of all frameworkatoms
#    ntime = len(allcoorframework)
#    nring = len(list_ring_indices)
#    centers = np.zeros((ntime,nring,3))
#    planes = np.zeros((ntime,nring,3))
#    #print "ntime,natom,3 (allcoorframework)",allcoorframework.shape
#    #print "ntime,nring,3 (centers)",centers.shape
#    for t in range(ntime):
#        for i in range(nring):
#            ring = Ring(list_ring_indices[i],allcoorframework[t,:,:],atomtypes,unitcell=unitcell,orient=list_orient[i])
#            centers[t,i,:] = ring.center[:]
#            planes[t,i,:] = ring.plane[:]
#    return centers,planes

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

def write_ksi_select_masked(fn_ksi,dist_ma,ksi,changedsign):
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
            print >> f, "%-6d  %10.4f  %10.4f %i" %(t,dist_ma[t,at],ksi[t,at],changedsign[t,at])
    f.close()

def write_ksi_select(fn_ksi,mask,dist,ksi,changedsign):
    # dist,ksi,mask -- arrays of size nstep x at
    # mask  --  array that gives True if the element should be masked (left out) 
    f = file(fn_ksi, "w+")

    nstep,natom = ksi.shape
    for at in range(natom):
        print >> f, "#at",at
        times = np.arange(nstep)[~mask[:,at]]   # gather the non-masked elements
        for i,t in enumerate(times[:-1]):
            print >> f, "%-6d  %10.4f  %10.4f %i" %(t,dist[t,at],ksi[t,at],changedsign[t,at])
    f.close()

########################################
### Functions: analyze ring passages
########################################

def create_summarytransitions(datadir,rg,runs=[1]):
    for i in range(rg.nring):
        for run in runs:
            # add one Passing instance for each run (to each of the rings)
            fn_ksi = "%s/ksi.%i.run%i.dat"%(datadir,i,run)
            atoms,time,dist,ksi,signchange = read_ksi_truncated(fn_ksi)
            shift = run*10000

            pas = Passing(atoms,time,dist,ksi,signchange,shift=shift)
            rg.list_rings[i].passing.append(pas)
    # rg will be adapted

def write_summarytransitions(logfile,rg):
    f = file(logfile,"w+")
    print >> f,"#ring run composition   transitions netto   transitions-corr netto-corr"
    for i in range(rg.nring):
        c = rg.list_compsrings[i]
        for run in range(len(rg.list_rings[i].passing)):
            print >> f, "%-3i %3i %3i   %6d %6d   %6d %6d" %(
                     i, c, run, rg.list_rings[i].passing[run].transitions, \
                     rg.list_rings[i].passing[run].netto, rg.list_rings[i].passing[run].transitions3,\
                     rg.list_rings[i].passing[run].netto3 )
    f.close()
    print "file written...",logfile

def write_averagetransitions(logfile,rg):
    if isinstance(logfile,str):
        f = file(logfile,"w+")
    else: f = logfile
    print >> f,"#ring composition   transitions netto   transitions-corr netto-corr"
    for i in range(rg.nring):
        c = rg.list_compsrings[i]
        pas = passing_average(rg.list_rings[i].passing, average="oneheap")
        print >> f,"%-3i %3i   %6d %6d   %6d %6d" %(
                 i,c,pas.transitions,pas.netto,pas.transitions3,pas.netto3 )
    if isinstance(logfile,str):
        f.close()
    print "file written...",logfile

#create_summarytransitions(datadir,rg,runs=np.arange(1,21))
#logfile = "%s.%s.transitions.dat" % (system,temp)
#write_summarytransitions(logfile,rg)
#logfile = "%s.%s.transitions.average.dat" % (system,temp)
#write_averagetransitions(logfile,rg)

#################################################################################

import matplotlib.pyplot as plt
colors = ['blue','green','red','black','grey','orange','magenta']

def plot_Fprofiles(basename,rg,):
    print "-"*20
    print "#plotting"
    print "#ring  composition transitions netto"
    
    plt.figure(1)
    plt.figure(2)
    for i in range(rg.nring):
        c = rg.list_compsrings[i]
        label = " ".join(str(a) for a in rg.list_comps[c]) if i==rg.list_rings_percomp[c][0]  else ""
        pas = passing_average(rg.list_rings[i].passing, average="oneheap")
        hist,edges = np.histogram(pas.ksi.ravel(),bins=50)
        # plot
        plt.figure(1)
        plt.subplot(3,1,1)
        plt.plot((edges[1:]+edges[:-1])/2.,hist,color=colors[c],label=label)
        plt.subplot(3,1,2)
        plt.plot((edges[1:]+edges[:-1])/2.,-np.log(hist)-min(-np.log(hist)),color=colors[c],)
        plt.subplot(3,1,3)
        plt.plot((edges[1:]+edges[:-1])/2.,-np.log(hist),color=colors[c],)
        plt.figure(2)
        plt.plot((edges[1:]+edges[:-1])/2.,-np.log(hist)-min(-np.log(hist)),color=colors[c],)
        print "ring",i,c,pas.transitions,pas.netto

    plt.figure(1)
    plt.subplot(3,1,1)
    # assume rg.ingred is set to ingredients, not to None
    plt.title(" ".join(rg.ingred))
    plt.legend()
    #plt.legend([" ".join(str(a) for a in comp) for comp in rg.list_comps])
    
    plt.subplot(3,1,1)
    plt.ylabel("hist(ksi)")
    plt.xlim(-4,4)
    plt.subplot(3,1,2)
    plt.ylabel("F(ksi) in [kBT]")
    plt.xlim(-4,4)
    plt.ylim(0,10)
    plt.grid()
    plt.subplot(3,1,3)
    plt.xlabel("ksi in [A]")
    plt.ylabel("F(ksi) in [kBT]")
    plt.xlim(-4,4)
    plt.grid()

    plt.savefig(basename+".png")

    plt.figure(2)
    plt.title(" ".join(rg.ingred))
    plt.legend()
    plt.grid()
    plt.xlim(-4,4)
    plt.xlabel("ksi in [A]")
    plt.ylabel("F(ksi) in [kBT]")

    plt.savefig(basename+".noshift.png")

def plot_Fprofiles_ringtypeidentical(basename,rg,):
    print "-"*20
    print "#plotting"
    print "#ring  composition transitions netto"

    for i in range(rg.nring):
        print "ring",i,"neighbors_identical",rg.list_rings[i].neighbors_identical

    plt.figure()
    for i in range(rg.nring):
        c = rg.list_rings[i].neighbors_identical
        # plot
        pas = passing_average(rg.list_rings[i].passing, average="oneheap")
        hist,edges = np.histogram(pas.ksi.ravel(),bins=np.arange(-4,4.,8./50.))
        plt.subplot(2,1,1)
        plt.plot((edges[1:]+edges[:-1])/2.,hist,color=colors[c])
        plt.subplot(2,1,2)
        plt.plot((edges[1:]+edges[:-1])/2.,-np.log(hist)-min(-np.log(hist)),color=colors[c])
        print "ring",i,c,pas.transitions,pas.netto

    plt.subplot(2,1,1)
    plt.title("neighbors_identical: False(blue) / True(green)")

    plt.subplot(2,1,1)
    plt.ylabel("hist(ksi)")
    plt.xlim(-4,4)
    plt.ylim(0,15000)
    plt.grid()
    plt.subplot(2,1,2)
    plt.xlabel("ksi in [A]")
    plt.ylabel("F(ksi) in [kBT]")
    plt.xlim(-4,4)
    plt.ylim(0,10)
    plt.grid()
    plt.savefig("%s.png"%(basename))

def plot_Fprofiles_perringtype(fignamebase,rg,):
    print "-"*20
    print "#plotting"
    print "#ring  composition transitions netto"

    for c,crings in enumerate(rg.list_rings_percomp):
        print "Doing comp",c,"for rings",crings

        plt.figure()
        for i,ringnb in enumerate(crings):
            # plot
            pas = passing_average(rg.list_rings[ringnb].passing, average="oneheap")
            hist,edges = np.histogram(pas.ksi.ravel(),bins=np.arange(-4,4.,8./50.))
            plt.subplot(2,1,1)
            plt.plot((edges[1:]+edges[:-1])/2.,hist,color=colors[c])
            plt.subplot(2,1,2)
            plt.plot((edges[1:]+edges[:-1])/2.,-np.log(hist)-min(-np.log(hist)),color=colors[c])
            print "ring",ringnb,c,pas.transitions,pas.netto
    
        plt.subplot(2,1,1)
        # assume rg.ingred is set to ingredients, not to None
        plt.title(" ".join(rg.ingred)+" "+" ".join(str(a) for a in rg.list_comps[c]))
    
        plt.subplot(2,1,1)
        plt.ylabel("hist(ksi)")
        plt.xlim(-4,4)
        plt.ylim(0,1000)
        plt.grid()
        plt.subplot(2,1,2)
        plt.xlabel("ksi in [A]")
        plt.ylabel("F(ksi) in [kBT]")
        plt.xlim(-4,4)
        plt.ylim(0,10)
        plt.grid()
        plt.savefig("%s.comp%i.png"%(fignamebase,c))


def write_Fprofiles(basename,rg,shift=False):
    print "-"*20
    print "writing Fprofiles"

    f = file(basename+".max.dat","w+")
    print >> f, "#ring Fmax"
    for i in range(rg.nring):
        pas = passing_average(rg.list_rings[i].passing, average="oneheap")
        hist,edges = np.histogram(pas.ksi.ravel(),bins=np.arange(-4,4.,8./50.))
        if shift:
            fun = -np.log(hist)-min(-np.log(hist))
        else:
            fun = -np.log(hist)
        filename = basename+".ring.%i.dat" %(i)
        write_histogram(filename,fun,edges)
        print >> f, i,max(fun)
        print i,max(fun)
    print "file written...",basename+".max.dat"
    f.close()

def write_histogram(filename,fun,edges):
    assert len(fun)+1==len(edges)
    middle = (edges[:-1]+edges[1:])/2.
    nbin = len(middle)
    f = file(filename,"w+")
    print >> f, "#bin fun edge-left edge-right middle"
    for i in range(nbin):
        print >> f, "%-4i %f %f %f %f" %(i,fun[i],edges[i],edges[i+1],middle[i],)
    print "file written...",filename
    f.close()

def write_xyz(filename,positions,atomtypes,labels=None):
    natom = len(positions)
    assert natom==len(atomtypes)
    f = file(filename,"w+")
    print >> f, natom
    print >> f,"title"
    for i in range(natom):
      if labels is None:
        print >> f, atomtypes[i], " ".join(str(val) for val in positions[i,:]) 
      else:
        print >> f, atomtypes[i], " ".join(str(val) for val in positions[i,:]), labels[i]
    f.close()
    print "file written...",filename


