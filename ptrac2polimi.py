from mcnptools import Ptrac
from sys import stdout
import numpy as np


mt2ntyn = { -1:1  ,\
            -2:2  ,\
            -4:5  ,\
             2:-99,\
           102:0  ,\
           107:0   }
for i in range(51,92): mt2ntyn[i] = -1



masses = {1001:0.99916733,\
          6000:11.9078563} 

qvals = {6000:\
              {107:-5.70205}}



def info(evt):
    """Return event information."""

    ipt = int(evt.Get(Ptrac.PARTICLE))
    rxn = int(evt.Get(Ptrac.RXN))
    za  = int(evt.Get(Ptrac.ZAID))
    dep = evt.Get(Ptrac.ENERGY)
    x   = evt.Get(Ptrac.X)
    y   = evt.Get(Ptrac.Y)
    z   = evt.Get(Ptrac.Z)
    tme = evt.Get(Ptrac.TIME)
    wgt = evt.Get(Ptrac.WEIGHT)

    return ipt,rxn,za,dep,x,y,z,tme,wgt


def record(nps,npar,ipt,rxn,za,cel,dep,tme,x,y,z,wgt,ngen,nsca,ncode,erg):
    """Return collision record string"""

    strng  = "{:10d}{:6d}{:3d}{:5d}{:6d}{:11d}".format(nps,npar,ipt,rxn,za,cel)
    strng += "{:11.6f}{:17.2f}{:9.2f}{:8.2f}{:8.2f}".format(dep,tme,x,y,z)
    strng += "{:11.3e}{:5d}{:6d}{:5d}{:10.3e}\n".format(wgt,ngen,nsca,ncode,erg)

    return strng


def convert(ptracFilename, polimiFilename="polimi", numCoinc=1, \
            cells=None, nEthresh=0.0, pEthresh=0.0, eEthresh=0.0, \
            verbose=False):
    """Converts an MCNP PTRAC file into an MCNPX/PoliMi file
       =====================================================
       Input: ptracFilename  != required, input filename
              polimiFilename != optional, output filename
              numCoinc       != optional, default all events
              cells          != optional, default all cells
              npeEthresh     != optional, default all energy
              verbose        != optional, information edits
       =====================================================
    """

    eThresh = [nEthresh,pEthresh,eEthresh]

    try:
        PT = Ptrac(ptracFilename,Ptrac.BIN_PTRAC)
    except RuntimeError:
        try:
            PT = Ptrac(ptracFilename,Ptrac.ASC_PTRAC)
        except RuntimeError:
            print("MCNPTOOLS doesn't recognize input file, {}: check format".format(ptracFilename))
            return

    f = open(polimiFilename,'w')
    hists = PT.ReadHistories(20049)
    banks = []
    terms = []
    colls = []
    while hists:
        for hist in hists:
            nps     = hist.GetNPS()
            numEvts = hist.GetNumEvents()
            evts = [hist.GetEvent(e) for e in range(numEvts)]

            strng = ""
            unqcl = []
            npar = 0
            for e in range(numEvts):
                evt = evts[e]
                cel = int(evt.Get(PT.CELL))
                detect = (cells is None) or (cel in cells)
                # Source Events
                if evt.Type() == PT.SRC: npar += 1
                # Bank Events
                if evt.Type() == PT.BNK:
                    if evt.BankType() == 0:
                        npar += 1
                    else:
                        if detect:
                            if verbose: stdout.write("BANK - {} {} {}\n".format(nps.NPS(),e,evt.BankType()))
                            banks.append(evt.BankType())
                # Collision Events
                if evt.Type() == PT.COL:
                    if detect:
                        erg = evts[e-1].Get(PT.ENERGY)
                        ipt,rxn,za,dep,x,y,z,tme,wgt = info(evt)
                        # Neutron collisions
                        if ipt == 1:
                            # Elastic scatter
                            if rxn == 2:
                                dep = erg-dep
                            # Inelastic scatter
                            elif (51 <= rxn) and (rxn <=91):
                                d1 = np.array([evts[e-1].Get(PT.U),evts[e-1].Get(PT.V),evts[e-1].Get(PT.W)])
                                d2 = np.array([evt.Get(PT.U),evt.Get(PT.V),evt.Get(PT.W)])
                                mu = np.dot(d1,d2) / (np.dot(d1,d1)*np.dot(d2,d2))**0.5
                                try:
                                    dep = (erg + dep - 2*mu*(erg*dep)**0.5)/masses[za]
                                except KeyError:
                                    stdout.write("WARNING: Inelastic Collision Missing Mass, {} {}".format(nps.NPS(),za))
                                    dep = erg
                            # Radiative capture
                            elif rxn == 102:
                                try:
                                    dep = erg/(masses[za]+1)
                                except KeyError:
                                    stdout.write("WARNING: Capture Collision Missing Mass, {} {}".format(nps.NPS(),za))
                                    dep = erg
                            # Alpha production
                            elif rxn == 107:
                                try:
                                    dep = erg + qvals[za][rxn]
                                except KeyError:
                                    stdout.write("WARNING: Alpha Production Missing Q-value, {} {} {}".format(nps.NPS(),za,rxn))
                                    dep = erg
                            # Other
                            else:
                                stdout.write("WARNING: Other Reaction Not Handled, {} {} {} {}\n".format(nps.NPS(),za,erg,rxn))
                                dep = erg-dep
                                colls.append(rxn)
                        # Photon collisions
                        if ipt == 2:
                            za = int(za/1000)
                            # Coherent/Incoherent scatter
                            if rxn == -1 or rxn == -2:
                                dep = erg-dep
                            # Pair production
                            elif rxn == -4:
                                dep = erg-1.022
                            # Other
                            else:
                                stdout.write("WARNING: Other Reaction Not Handled, {} {} {} {}\n".format(nps.NPS(),za,erg,rxn))
                                dep = erg-dep
                                colls.append(rxn)

                        # Converting MT to PoliMi NTYN
                        try:
                            rxn = mt2ntyn[rxn]
                        except KeyError:
                            stdout.write("WARNING: MT to NTYN Conversion Failed for rxn={}, nps={}\n".format(rxn,nps.NPS()))
                            pass

                        if dep > eThresh[ipt-1]:
                            unqcl.append(cel)
                            strng += record(nps.NPS(),npar,ipt,rxn,za,cel,dep,tme,x,y,z,wgt,0,0,0,erg)

                        if verbose: stdout.write("COLL - {} {} {}\n".format(nps.NPS(),e,rxn))
                # Termination Events
                if evt.Type() == PT.TER:
                    if detect:
                        if verbose: stdout.write("TERM - {} {} {}\n".format(nps.NPS(),e,evt.Get(PT.TERMINATION_TYPE)))
                        terms.append(evt.Get(PT.TERMINATION_TYPE))
            # Record Collision Event
            if len(set(unqcl)) >= numCoinc:
                if verbose: stdout.write(strng)
                f.write(strng)
        hists = False #PT.ReadHistories(100)

    f.close()
    stdout.write("Collision Types: {}\nBank Events Types: {}\nTermination Types: {}\n".\
                 format(set(colls),set(banks),set(terms)))

#==============================================================================
# Execute the following if running script from command line: python ptrac2polimi.py
if __name__ == "__main__":

    from sys import argv
    import os.path

    if len(argv)>1:
        ptracFilename = argv[1]
    else:
        ptracFilename = "ptrac"

    if os.path.isfile(ptracFilename):
        stdout.write("Converting "+ptracFilename+"...\n")
#        convert(ptracFilename)
        convert(ptracFilename,nEthresh=0.01,pEthresh=0.01,numCoinc=2, \
                cells=[3,14,25,36,47,58,69,80,91,102,113,124,135,146])#,verbose=True)
    else:
        stdout.write("Cannot find file named "+ptracFilename+", check input args.\n")


