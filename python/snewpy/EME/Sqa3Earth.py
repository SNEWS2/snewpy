#!/usr/bin/env python3

if __name__ == "__main__":

    import numpy as np
    from astropy import units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz

    # skycoordinates of neutrino source
    Betelgeuse = SkyCoord.from_name('Betelgeuse') 
    
    # neutrino detector
    SuperK = EarthLocation.of_site('SuperK')

    # time when the supernova occured 
    # the first time option means the neutrinos traveled through the Earth, the second means they did not
    time = Time('2021-5-26 14:14:00')
    #time = Time('2021-5-26 14:14:00') - 12*u.hour    

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

    # if set to True the Sqa3Earth module will output files in the 'out' directory
    # The stepcounterlimit controls how often output is written. The larger the number, the less often it happens. 
    ID.outputflag = False 
    ID.stepcounterlimit = 1



    # do the calculation. The return is a four dimensional array of transition probabilities nu_alpha -> nu_i: 
    # index order is matter/antimatter, energy, i, alpha
    Pfm = Sqa3Earth.RunSqa3Earth(ID)

    print("finished")

