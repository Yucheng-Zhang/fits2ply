import numpy as np
import healpy as hp


def get_ra_dec(theta, phi):
    '''Get RA, DEC [degree] from theta, phi [radians] used in Healpy.'''
    rot = hp.Rotator(coord=['G', 'C'])
    theta_equ, phi_equ = rot(theta, phi)
    dec, ra = 90. - np.rad2deg(theta_equ), np.rad2deg(phi_equ)
    # move RA in [-180,0) to [180,360)
    ra = np.where(ra < 0., ra + 360., ra)

    return ra, dec


def get_theta_phi(ra, dec):
    '''Get theta, phi [radians] used in Healpy from RA, DEC [degree].'''
    rot = hp.Rotator(coord=['C', 'G'])
    theta_equ, phi_equ = np.deg2rad(90.-dec), np.deg2rad(ra)
    theta_gal, phi_gal = rot(theta_equ, phi_equ)
    # move phi in [-pi,0) to [pi,2pi)
    phi_gal = np.where(phi_gal < 0., phi_gal + 2*np.pi, phi_gal)

    return theta_gal, phi_gal
