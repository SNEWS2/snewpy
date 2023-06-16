#!/usr/bin/env python3

if __name__ == "__main__":

    import numpy as np
    from astropy import units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz

    # skycoordinates of neutrino source
    Betelgeuse = SkyCoord.from_name('Betelgeuse') 
    
    # neutrino detector
    SuperK = EarthLocation(lat=36.425722*u.deg, lon=137.310306*u.deg, height=389*u.m)
    #SuperK = EarthLocation.of_site('Super-Kamiokande')
    utcoffset = +9*u.hour  

    # when the supernova occured
    time = Time('2021-5-26 23:14:00') - utcoffset

    # altaz of supernovae at detector
    SNaltaz = Betelgeuse.transform_to(AltAz(obstime=time,location=SuperK)) 

    # load the mdule that does the Earth-matter effect calculation 
    import Sqa3Earth

    # class to accumulate input data for calcultion
    ID = Sqa3Earth.InputDataSqa3Earth()

    # assign data fields
    ID.altitude = SNaltaz.alt.deg
    ID.azimuth = SNaltaz.az.deg

    ID.outputfilenamestem = "out/Sqa3Earth:PREM"

    ID.densityprofile = "PREM.rho.dat"
    ID.electronfraction = "PREM.Ye.dat"

    ID.NE = 296
    ID.Emin = 1	
    ID.Emax = 60

    ID.deltam_21 = 7.59E-005
    ID.deltam_32 = 2.43E-003
    ID.theta12 = 34.4
    ID.theta13 = 9E-000	
    ID.theta23 = 45
    ID.deltaCP = 0	

    ID.accuracy = 1.01E-009
    ID.stepcounterlimit = 1000

    ID.outputflag = False

    # do the calculation. The return is a four dimensional array of transition probabilities nu_alpha -> nu_i: 
    # index order is matter/antimatter, energy, i, alpha
    Pfm = Sqa3Earth.RunSqa3Earth(ID)

    print(Pfm[0][0][0][0])

