#! /usr/bin/env python
from __future__ import print_function
import pyfplo.wanniertools as wt


# ===================================================================

if __name__ == '__main__':


    wdc=wt.WanDefCreator(rcutoff=27,wftol=0.0005,coeffformat='bin',
                         wfgriddirections=[[1,0,0],[0,1,0],[0,0,1]],
                         wfgridsubdiv=[1,1,1],savespininfo=False)

    emin=-13
    emax= 0
    delower=1
    deupper=5  

    wdc.add(wt.MultipleOrbitalWandef('Si',[1,2],['3sb','3pb'],
                                  emin=emin,emax=emax,
                                  delower=delower,deupper=deupper))
    


    wdc.writeFile('=.wandef')

    

