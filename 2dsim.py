#!/usr/bin/env python
# 
# Author: Tiange Cui (tiange_cui@vrtx.com)
#
# Copyright (C) 2016 Vertex. All rights reserved.
#
# To run: python omega_rocs.py input output output_directory
#
# Function: Add 2D Tanimoto similarity score to the FTrees/ROCS result (.sdf file).
#
import os
import sys
import fnmatch
from openeye.oechem import *
from openeye.oegraphsim import *

filename = sys.argv[1]
output_rocs = sys.argv[3]
regex = r'{0}*.*'.format(filename)

for file in os.listdir(output_rocs):
    if fnmatch.fnmatch(file, regex):
        filename = file

if len(sys.argv) != 4:
    OEThrow.Usage("%s <inputfile> <outputfile> <output_directory>" % sys.argv[0])

ifs = oemolistream(filename)
if not ifs.open(filename):
    OEThrow.Fatal("Unable to open %s for reading" % filename)

ofs = oemolostream()
if not ofs.open(sys.argv[2]):
    OEThrow.Fatal("Unable to open %s for writing" % sys.argv[2])
if ofs.GetFormat() != OEFormat_SDF:
    OEThrow.Fatal("%s output file has to be an SDF file" % sys.argv[2])

qfp = OEFingerPrint()
mo = OEGraphMol()
OEReadMolecule(ifs,mo) # read first mol (reference)
OEMakeFP(qfp, mo, OEFPType_MACCS166)
mo.SetData("Tanimoto", qfp)

ifs = oemolistream(filename)
if not ifs.open(filename):
    OEThrow.Fatal("Unable to open %s for reading" % filename)
	
fp = OEFingerPrint()
for mol in ifs.GetOEGraphMols():
    OEMakeFP(fp, mol, OEFPType_MACCS166)
    tanimotoscore = OETanimoto(qfp, fp)
    fptypestr = "Tanimoto similarity"
    fpdata = str('{:05.3f}'.format(tanimotoscore))
    OESetSDData(mol, fptypestr, fpdata)
    OEWriteMolecule(ofs, mol)