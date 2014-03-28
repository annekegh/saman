#!/bin/python

"""Analyze ksi-coordinate
AG, 19 March, 2014"""

import numpy as np
# for matplotlib figure generation...
import matplotlib
matplotlib.use('Agg')


###################################################################################

from saman.ringfunc import read_coords,read_rings
from saman.ringfunc import read_ksi_truncated, passing_average
from saman.ringfunc import Passing, RingGroup
import sys,glob

system = sys.argv[1]
temp = sys.argv[2]
datadir = sys.argv[3]
#list_datadir = glob.glob("../data-channels/%s_%sK/%s_*ethylene_%sK/" %(system,temp,system,temp))
list_datadir = glob.glob("%s/*/"%datadir)
assert len(list_datadir) == 1
datadir = list_datadir[0]

# ring definitions
#somedir = "../channels/ringdefinitions/"
somedir = sys.argv[4]
fn_xyz = "%s/config.xyz" %(somedir,)
fn_rings = "%s/ringatoms" %(somedir,)
#fn_rings_atomtypes = "../%s/%s/atomtypes" %(somedir,system)


## set datadir
#system = "sapo34"
#temp = "450"
#temp = "300"
#import glob
#list_datadir = glob.glob("../data-channels/%s_%sK/%s_*ethylene_%sK/" %(system,temp,system,temp))
#list_datadir = glob.glob("/home/an/%s_*ethylene_%sK/" %(system,temp))
#assert len(list_datadir) == 1
#datadir = list_datadir[0]
##datadir = "../data-channels/%s_450K/%s_26ethylene_450K/" %(system,system)
#
## ring definitions
#fn_xyz = "../zeoTsites/results_sam/%s/config.xyz" %system
#fn_rings = "../zeoTsites/results_sam/%s/ringatoms" %system
#fn_rings_atomtypes = "../zeoTsites/results_sam/%s/atomtypes" %system

#################################################################################

natom,atomtypes,pos = read_coords(fn_xyz)
list_rings_indices = read_rings(fn_rings)

rg = RingGroup(list_rings_indices,pos,atomtypes)
rg.set_ellips()
rg.set_passing()

#################################################################################
print "*"*20
print "Analyze composition"
print "*"*20
#ingredients = None
ingredients = ["Al","O1","O2","P","Si",]
rg.analyze_composition(ingredients=ingredients)
rg.print_composition()

print "*"*20



#################################################################################

logfile = "%s.%s.log.dat" % (system,temp)
f = file(logfile,"w+")
print >> f,"#ring run transitions netto"


nfile = 20   # loop over rings

for i in range(rg.nring):
    for n in np.arange(1,nfile+1):
        fn_ksi = "%s/ksi.%i.run%i.dat"%(datadir,i,n)
        atoms,time,dist,ksi,signchange = read_ksi_truncated(fn_ksi)
        shift = n*10000

        pas = Passing(atoms,time,dist,ksi,signchange,shift=shift)
        rg.list_rings[i].passing.append(pas)
        
        #print >> f, "ring",i,"file",n,"transitions",np.sum(abs(signchange)),np.sum(signchange)
        print >> f, i, n, pas.transitions, pas.netto
f.close()


print "-"*20
print "#ring  composition transitions netto"
for i in range(rg.nring):
    c = rg.list_compsrings[i]
    # plot
    pas = passing_average(rg.list_rings[i].passing, average="oneheap")
    print "ring",i,c,pas.transitions,pas.netto
print "-"*20

#################################################################################
import matplotlib.pyplot as plt

colors = ['blue','green','red','black','grey','orange']
plt.figure()
for i in range(rg.nring):
    c = rg.list_compsrings[i]
    # plot
    pas = passing_average(rg.list_rings[i].passing, average="oneheap")
    hist,edges = np.histogram(pas.ksi.ravel(),bins=50)
    plt.subplot(2,1,1)
    plt.plot((edges[1:]+edges[:-1])/2.,hist,color=colors[c])
    plt.subplot(2,1,2)
    plt.plot((edges[1:]+edges[:-1])/2.,-np.log(hist)-min(-np.log(hist)),color=colors[c])

plt.subplot(2,1,1)
plt.title(" ".join(ingredients))
plt.legend([" ".join(str(a) for a in comp) for comp in rg.list_comps])

plt.subplot(2,1,1)
plt.xlabel("ksi in [A]")
plt.ylabel("hist(ksi)")
plt.xlim(-4,4)
plt.subplot(2,1,2)
plt.xlabel("ksi in [A]")
plt.ylabel("F(ksi) in [kBT]")
plt.xlim(-4,4)
plt.savefig("histogram.ksi_ma.png")

  #  import numpy.ma as ma
  #  dist_ma = ma.masked_greater(dist,5.)
  #  ksi_ma = ma.masked_array(ksi,mask=dist_ma.mask)
  #  
  #  #print "masked dist:",sum(dist_ma.mask)  # this is the number of "True" values = the number of masked values
  #  #print "masked ksi:",sum(ksi_ma.mask)  # this is the number of "True" values = the number of masked values
  #  
  #  # detect ring passage
  #  sign = (ksi>0)       # True if ksi>0, False if ksi<0
  #  signchange = np.array(sign[1:,:],int)-np.array(sign[:-1,:],int)
  #  dksi_ma = ma.masked_array(signchange,mask=dist_ma.mask[:-1])
  #  hist,edges = np.histogram(ksi_ma[~ksi_ma.mask].ravel(),bins=50)

