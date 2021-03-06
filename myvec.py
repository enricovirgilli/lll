"""Vector algebra"""
import sys
from math import sqrt, acos, sin, cos, atan2
from pylab import array, dot, norm
################################################################################
# Vectors
#
def normalize(v):
    #Returns the normal vector parallel to a n-dimensional vector
    epsilon = 1e-100
    return array(v)/(norm(v)+epsilon)

def directorcos(v, w):
    """Directo cosine between 2 n-dimensional vectors."""
    return dot(v, w)/(norm(v)*norm(w))

def anglebetween(v, w):
    """Angle between 2 n-dimensional vectors."""
    dcos=directorcos(v, w)
    if dcos>1: dcos=1.
    if dcos<-1: dcos=-1.
    return acos(dcos)

def makenormalvec(phi, theta):
    """Returns a normal vector in cartesian coordinates, given the two angular coordinates
    """
    x=cos(phi)*sin(theta)
    y=sin(phi)*sin(theta)
    z=cos(theta)
    return array((x,y,z))

def spherical2cartesian(rho, phi, theta):
    """Returns a vector in cartesian coordinates, given the three spherical ones
    """
    return rho*makenormalvec(phi, theta)

def cartesian2spherical(*v):
    """Returns a vector in cartesian coordinates, given the three spherical ones
    """
    if len(v)==1: v=v[0]
    r=norm(v)
    phi=atan2(v[1], v[0])
    if len(v)<3: return r, phi
    return r, phi, v[2]/r

def rotx(v, alpha):
    '''Rotation of angle alpha around the x axys
    '''
    v_=array(v)
    c, s = cos(alpha), sin(alpha)
    line0=(1., 0., 0.)
    line1=(0.,  c,  s)
    line2=(0., -s,  c)
    return (sum(line0*v_), sum(line1*v_), sum(line2*v_))
#
def roty(v, alpha):
    '''Rotation of angle alpha around the y axys
    '''
    v_=array(v)
    c, s = cos(alpha), -sin(alpha)
    line0=( c, 0., -s)
    line1=(0., 1., 0.)
    line2=( s, 0.,  c)
    return (sum(line0*v_), sum(line1*v_), sum(line2*v_))
#
def rotz(v, alpha):
    '''Rotation of angle alpha around the z axys
    '''
    v_=array(v)
    c, s = cos(alpha), sin(alpha)
    line0=( c,  s, 0.)
    line1=(-s,  c, 0.)
    line2=(0., 0., 1.)
    return (sum(line0*v_), sum(line1*v_), sum(line2*v_))
#
def rot(v, rotations):
    '''Rotates the crystal of angles theta and phi
    '''
    rotdict={'x':rotx, 'y':roty, 'z':rotz}
    for rotation in rotations:
        axis, angle = rotation
        v = rotdict[axis](v, angle)
    return v
#
if __name__=='__main__':
    for x in arange(-10., 10.):
        for y in arange(-10., 10.):
            if atan2(y,x)-azimuth((x,y))!=0:
                print x,y,atan2(y,x)/pi-azimuth((x,y))/pi
    pass
