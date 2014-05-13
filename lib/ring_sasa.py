#!/usr/bin/env python

import numpy as np

def area_polygon(a):
    #it is assumed that the vertices in a are ordered clockwise or counterclockwise,
    #and that the polygon is not intersecting
    ar = np.roll(a, -1, axis = 0)
    area  = 0.5 * abs(np.sum(a[i,0]*ar[i,1]-ar[i,0]*a[i,1] for i in xrange(len(a))))
    return area

def intersect_2circles(a,b,ra,rb):
    #see: http://2000clicks.com/mathhelp/GeometryConicSectionCircleIntersection.aspx
    #a and b are arrays representing the coordinates of the centers of 2 circles
    #K is the area of the triangle formed by the centers of the two circles and one of their points of intersection;
    #d is the distance between the circles centers;
    #ra, rb are the circles' radii;
    #(xa,ya) and (xb,yb) are the circles' centers
    #example: c, d = intersect_2circles(a,b,ra,rb)
    xa = a[0]
    ya = a[1]
    xb = b[0]
    yb = b[1]
    d2 = (xa - xb)**2 + (ya - yb)**2
    if d2 > (ra + rb)**2 : #no intersection
        return np.empty((2,2)) * np.nan
    K = 0.25 * np.sqrt(((ra + rb)**2 - d2)*(d2 - (ra - rb)**2))
    x1 = 0.5*(xb + xa) + 0.5*(xb - xa)*(ra**2 - rb**2)/d2 + 2*(yb - ya)*K/d2
    x2 = 0.5*(xb + xa) + 0.5*(xb - xa)*(ra**2 - rb**2)/d2 - 2*(yb - ya)*K/d2
    y1 = 0.5*(yb + ya) + 0.5*(yb - ya)*(ra**2 - rb**2)/d2 - 2*(xb - xa)*K/d2
    y2 = 0.5*(yb + ya) + 0.5*(yb - ya)*(ra**2 - rb**2)/d2 + 2*(xb - xa)*K/d2
    return np.array([[x1,y1],[x2,y2]])

def select_shortest(a):
    #select from a set of pairs of 2D positions the position that has the shortest vector (closest to the origin)
    #a must have shape (N,2,2) : N pairs, 2 positions per pair, 2 coordinates (x and y) per position
    dist2 = np.sum(a**2,axis=2)
    ordered = np.argsort(dist2,axis=1)
    return a[ordered==0]

def get_intersect(a,sasaradii,natom,step):
    #calculate for a set of atoms (a) the intersections with neighbors located at (step) distance
    a_roll = np.roll(a, -step, axis = 0)
    sasaradii_roll = np.roll(sasaradii, -step)
    return np.array([intersect_2circles(a[i],a_roll[i],sasaradii[i],sasaradii_roll[i]) for i in xrange(natom)])

def remove_between_next_neighbors(a,sasaradii,natom):
    #check if next-neighbor intersections overlap with the 1 atom in between
    #if that is not the case remove that 1 atom
    #for each pair of intersections select the one that is closest to the ring center
    intersect = get_intersect(a,sasaradii,natom,2)
    newintersect = select_shortest(intersect)
#     print "newintersect",newintersect
    sasaradii_roll = np.roll(sasaradii, -1)
    a_roll = np.roll(a, -1, axis = 0)
    toremove = np.array([False]*natom)
    for i in xrange(natom):
        if np.sum((newintersect[i]-a_roll[i])**2) >= sasaradii_roll[i]**2:
            remove_index1 = (i+1)%natom
            toremove[(i+1)%natom] = True
#             print "i,remove_index1",i,remove_index1
#     print "atoms remaining (of the",natom,")",[i for i in xrange(natom) if toremove[i] == False]
    if np.sum(toremove) > 0:
        a = np.array([a[i] for i in xrange(natom) if toremove[i] == False])
        sasaradii = np.array([sasaradii[i] for i in xrange(natom) if toremove[i] == False])
        natom = len(a)
    return a,sasaradii,natom

def ring_sasa(a,radii,solventradius):
    #calculate that solvent accessible surface area in a 2D ring
    #it is assumed that the ringatoms are:
        #ordered clockwise or counter-clockwise
        #aligend with the xy plane
        #centered around the origin

    natom = len(a)
    sasaradii = radii + solventradius

    #step1: check if any diagonal atom pair distance is shorter than the sum of their sasaradii
    #in that case, the area is zero
#     print "step 1"
    diameters2 = np.sum((a[:natom/2] - a[natom/2:])**2,axis=1)
    sasasum2 = (sasaradii[:natom/2] + sasaradii[natom/2:])**2
    sasadiameters2 = diameters2 - sasasum2
    if np.min(sasadiameters2) <= 0.:
        return np.array([]),0.0


    #step 2:check if next-neighbor intersections overlap with the 1 atom in between
    #if that is not the case remove that 1 atom
    #do this step twice
#     print "step 2"
    a,sasaradii,natom = remove_between_next_neighbors(a,sasaradii,natom)
#     print "natom",natom
    a,sasaradii,natom = remove_between_next_neighbors(a,sasaradii,natom)
#     print "natom",natom
    if natom < 3:
        return np.array([]),0.0

    #step 3: for each pair of (remaining) neighbor rings calculate the intersection point closest to the ring center
#     print "step 3"
    intersect = get_intersect(a,sasaradii,natom,1)
    newintersect = select_shortest(intersect)
#     print "newintersect",newintersect

    #step 4:calculate the sasa area of the polygon formed by the intersection points as vertices
#     print "step 4"
    area = area_polygon(newintersect)

    return newintersect, area

if __name__ == "__main__":
    from yaff import *
    import sys
    from saman.ringfunc import write_xyz
    import time
    coord=sys.argv[1]
    sys1=System.from_file(coord)
    a1=sys1.pos/angstrom
    natom = len(a1)
    radii = np.array([1.4]*natom)
    solventradius = 1.5

    start = time.clock()
    points,area = ring_sasa(a1[:,0:2],radii,solventradius)
    end = time.clock()

    npoint = len(points)
    if npoint == 0:
        print "no points, exiting."
        sys.exit()
    atomtypes = np.array([12]*npoint)
    zeros = np.zeros(npoint).reshape(npoint,1)
    points3D = np.append(points,zeros,axis=1)
    print "points3D",points3D
    print "area",area
    print "time",end-start
    write_xyz(coord + "points.xyz",points3D,atomtypes)

