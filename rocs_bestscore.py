#!/usr/bin/env python
# 
# Author: Tiange Cui (tiange_cui@vrtx.com)
#
# Copyright (C) 2016 Vertex. All rights reserved.
#
# To run: python omega_rocs.py input output
#

import re
import os
import sys
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeshape import *
    
def omegagen(smile):
    omega = OEOmega()
    omega.SetMaxConfs(25)
    omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(False)
    mol = OEMol()
    OESmilesToMol(mol, smile)    
    ok_omega = omega(mol)
    return mol, ok_omega    

def rocs_bestscore(mol1,mol2):
    #omega = OEOmega()
    #omega.SetMaxConfs(25)
    #omega.SetStrictStereo(False)
    #if (omega(mol1)==True) and (omega(mol2)==True):
    best = OEBestOverlay()
    best.SetRefMol(mol1)
    best.SetColorForceField(OEColorFFType_ImplicitMillsDean)
    best.SetColorOptimize(True)
    scoreiter = OEBestOverlayScoreIter()
    OESortOverlayScores(scoreiter, best.Overlay(mol2), OEHighestTanimotoCombo())
           
    for score in scoreiter:
        #print "FitConfIdx: %-4d RefConfIdx: %-4d Tanimoto: %.2f" % (score.fitconfidx, score.refconfidx, score.GetTanimotoCombo())
        bestscore = score.GetTanimotoCombo()
        break
    #else:
    #    bestscore = -1.00
    return bestscore

def main():
    outputfile = ".".join((sys.argv[2], "txt"))
    if os.path.exists(outputfile):
        os.remove(outputfile)
    input = open(sys.argv[1], 'r')
    for i, line in enumerate(input):
        pairs=re.split(r' ', line)
        smi1 = pairs[0]
        smi2 = pairs[2]
        chemid1 = pairs[1]
        chemid2 = pairs[3]
        mol1, ok1 = omegagen(smi1)
        mol2, ok2 = omegagen(smi2)
        if (ok1==True) and (ok2==True):
            score = rocs_bestscore(mol1,mol2)
        else:
            score = -1.00
            
        with open(outputfile, 'a') as output:
            print "writing...Line %d", i
            newline = []
            newline = " ".join((smi1, chemid1, smi2, chemid2, pairs[4], pairs[5], pairs[6], '%0.2f'%(score), pairs[7]))
            output.write(newline)

if __name__ == "__main__":
    main()