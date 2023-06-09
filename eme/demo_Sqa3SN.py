#!/usr/bin/env python3

if __name__ == "__main__":

    import numpy as np
    from astropy import units as u
    # load the mdule that does the Earth-matter effect calculation 
    import Sqa3SN

    ID = Sqa3SN.InputDataSqa3SN()

    ID.outputfilenamestem = "./out/Sqa3SN:"

    ID.rmin = 1e7
    ID.rmax = 1e12

    ID.densityprofile = "rhoSN.dat"
    ID.electronfraction = "YeSN.dat"

    ID.NE = 500
    ID.Emin = 0.2	
    ID.Emax = 100

    ID.deltam_21 = 7.59E-005
    ID.deltam_32 = 2.43E-003
    ID.theta12 = 34.4
    ID.theta13 = 9E-000	
    ID.theta23 = 45
    ID.deltaCP = 0	

    ID.accuracy = 1.01E-009
    ID.stepcounterlimit = 10

    ID.outputflag = True

    # do the calculation, The return is a four dimensional array of transition probabilities nu_alpha -> nu_i: 
    # index order is matter/antimatter, energy, i, alpha

    Pmf = Sqa3SN.RunSqa3SN(ID)

    # the return is a four dimensional array of transition probabilities nu_alpha -> nu_i: 
    # index order is matter/antimatter, energy, i, alpha

