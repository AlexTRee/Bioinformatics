#!/bin/env python
# 
# Author: Tiange Cui (tiange_cui@vrtx.com)
#
# Copyright (C) 2016 Vertex. All rights reserved.
#
# This file is part of VertexScriptLibrary.
#
"""
ftrees-fs.py - Perform fuzzy similarity search to identify similar compounds using FTrees-FS (Fragment Space).

Usage:
    ftrees.py [-s <minSimilarityThreshold>] [-n <maxNofResults> ] [-i <infile>] [-l <libfile>] [-d <outputDir>]
    ftrees.py [-s <minSimilarityThreshold>] [-n <maxNofResults> ] [-i <infile>] [-l <libfile>] [-c <coresToUse>] [-q <query3D>] [-b <bestHitsNo>] [-d <outputDir>]
    ftrees.py -h | --help | -v | --version | -e | --examples

Category:
    2D searching and fingerprints

Description:
    Perform substructure and similarity searching against given
    fragment space using specified MDL, MOL or SMILE query file by the
	FTrees-FS search. OMEGA is then used for generating  multi-conformer 
    structure databases. ROCS is also used for identifying potentially 
    active compounds by shape comparison. 2D Tanimoto similarity score 
    is also added to the final result.

Options:
    -c, --coresToUse <coresToUse>  [default: 15]
        Number of CPU cores to use for running OMEGA and ROCS. [1 to 50]
    -q, --query3D <query3D>
        Provide query file with 3D coordinates in '.sdf' or '.mol2' format.
        Note: If query3D is provided, the script will use OMEGA program to generate 
        multi-conformer OEBinary file, then use ROCS program to overlay collection 
        of multi-conformers onto the query (reference) molecule.
    -b, --bestHitsNo <bestHitsNo>  [default: 500] 
        Maximum number of hits to return in ROCS results.
    -s, --minSimilarityThreshold <minSimilarityThreshold>  [default: 1]
        Minimun similarity threshold. Anything above this minimum similarity score
        will be reported. [0.0 to 1.0]
    -n, --maxNofResults <maxNofResults>  [default: 100]
        Maximum number of hits to return in FTrees-FS search. [1 to 100000]
    -i, --infile <infile>  [default: stdin]
        The query file in .smi format only. 
        The recommended file format is:
        SMILEs/SMARTs QueryName (seperated by one whitespace)
		eg.
		C1C[C@H](N(C1)C(=O)C(F)(F)F)C(=O)Cl 491672
    -l, --libfile <libfile>
        The library file you want to search against. (only .fsf type is allowed)
    -d, --outputDir <outputDir> [default: /cluster/home/cuit/results/FTrees/]
        Location to store the result files.
    -e, --examples
        Print examples retrieved from the examples section.
    -h, --help
        Print this help message.
    -v, --version
        Print version number.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.
        
Examples:
    To perform FTrees-FS similarity search using molecule query.smi against the defined library 
	library.smi for 500 similar molecule hits above 0.85 similarity threshold, including 
    base matching images in .pdf format, the query file with 3D coordinates is query3D.sdf. 
    Use 15 cores to run OMEGA and ROCS, the number of hits for ROCS is set to 1000, type:

        $ ftrees-fs.py -s 0.85 -n 500 -i query.smi -l lib.fsf -c 15 -q query3D.sdf -b 1000 -d /path/for/results

Author:
    Tiange Cui (tiange_cui@vrtx.com)

See also:
    ftrees.py
    2dsim.py
    
Copyright:
    Copyright (C) 2016 Vertex. All rights reserved.

    This file is part of VertexScriptLibrary.

"""
import os
import re
import sys
from docopt import docopt

try:
    import vertexutil
except ImportError, errMsg:
    sys.stderr.write("\nImportError: %s\n" % errMsg)
    sys.stderr.write("\nIt looks like your VertexScriptLibrary (VSL) environment is not set up correctly.\n")
    sys.stderr.write("Please execute the following command and try again: \"module load vsl\"\n\n")
    sys.exit(1)

__version__ = "1.0"

scriptName = os.path.basename(sys.argv[0])
optionsInfo = {}
optionsInfo = docopt(__doc__, version=__version__)

####################################################################
########## Here defines the location for the result files.##########
####################################################################
outputDir = optionsInfo["--outputDir"]
if (outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    if outputDir[-1] == "/":
        outputDir = outputDir[:-1]
    output_ftrees_fs= outputDir+"/FTrees-FS/"   # FTrees-FS results
    output_rocs = outputDir+"/rocs/"            # ROCS results with 2D Tanimoto similarity scores
    output_omegarocs = outputDir+"/omegarocs/"  # ROCS results without 2D Tanimoto similarity scores
else:    
    output_ftrees_fs="/cluster/home/cuit/results/FTrees-FS/"    # FTrees-FS results
    output_rocs = "/cluster/home/cuit/results/rocs/"            # ROCS results with 2D Tanimoto similarity scores
    output_omegarocs = "/cluster/home/cuit/results/omegarocs/"  # ROCS results without 2D Tanimoto similarity scores


        
def main():
    """Start script execution"""
    
    vertexutil.info("\n%s: Starting...\n" % scriptName)

    (wallClockTime, processorTime) = vertexutil.getWallClockAndProcessorTime()

    # Process command line arguments and options...
    processOptions()

    if not os.path.exists(output_ftrees_fs):
        os.makedirs(output_ftrees_fs)
    if not os.path.exists(output_omegarocs):
        os.makedirs(output_omegarocs)
    if not os.path.exists(output_rocs):
        os.makedirs(output_rocs)
    cmd = "rm -r "+output_ftrees_fs+"*"
    os.system(cmd)
    cmd = "rm -r "+output_omegarocs+"*"
    os.system(cmd)
    cmd = "rm -r "+output_rocs+"*"
    os.system(cmd)

    # Perform appropriate actions required by the script...
    processScriptActions()
    
    vertexutil.info("\n%s: Done...\n" % scriptName)
    vertexutil.info("Total time: %s\n" % vertexutil.getFormattedElapsedTime(wallClockTime, processorTime))
    vertexutil.info("The FTrees results can be found at: %s\n" % output_ftrees_fs)
    vertexutil.info("The ROCS results can be found at: %s\n" % output_rocs)
    
def processScriptActions():
    """FTrees3 search using SMILES query against selected library"""
    
    infile = optionsInfo["--infile"]
    libfile = optionsInfo["--libfile"]
    minSimilarityThreshold = float(optionsInfo["--minSimilarityThreshold"])
    maxNofResults = int(optionsInfo["--maxNofResults"])
    coresToUse = int(optionsInfo["--coresToUse"])
    query3D = optionsInfo["--query3D"]
    bestHitsNo = int(optionsInfo["--bestHitsNo"])
    outputDir = optionsInfo["--outputDir"]
    
    setSim="SET FFS_TARGET_SIMILARITY "+str(minSimilarityThreshold)
    timestamp=os.system("date \"+%H%M%S\"")
    outputfile = output_ftrees_fs+"ftrees-fs.bat"
    absqueryname = infile.split(".")[0]         # Get query file name with absolute path
    queryname=os.path.basename(absqueryname)    # Get query file name without absolute path
    abslibname = libfile.split(".")[0]          # Get library file name with absolute path
    libname=os.path.basename(abslibname)        # Get library file name without absolute path
    outfiles=output_ftrees_fs+queryname+"_results_$(query_nr)"

    with open(outputfile, 'w') as output:
        content = """MOL
CONVERT \""""+infile+"""\" % % \""""+output_ftrees_fs+queryname+""".fdf\" n
END   
"""+setSim+"""

FTREES
  DELETE 0 %
  FFSPACE
    DELETE  
    SELINIT !*
    READ \""""+libname+""".fsf\" n y
  END
  READ 0 \""""+output_ftrees_fs+queryname+""".fdf\" % %
  SETVAR $(last_query) $(LAST_SLOT)
  FFSEARCH
    FOR_EACH $(query_nr) FROMTO 0 $(last_query)
      MCREATE $(query_nr) % %
      SELOUTP """+outfiles+""" o y
      MSEARCH """+str(maxNofResults)+""" 20 y
      SELOUTP screen
      2MOLS """+outfiles+""".mol2 y 1-"""+str(maxNofResults)+"""
      SELOUTP """+outfiles+""".smi o y
      2SMILES 1-"""+str(maxNofResults)+""" n n y n s
      SELOUTP screen
    END_FOR
  END
END"""
        output.write(content)
    
    # Perform FTrees-FS search
    cmd = "ftrees -b "+outputfile
    os.system(cmd)
    
    # If query3D is provided, multi-conformer overlay will be performed as well as FTrees search.
    if (query3D):
        for file in os.listdir(output_ftrees_fs):
            ext = re.search(r'(\w+)\.(\w+)', file)
            if (ext):
                ext = ext.group(2)
                if re.match("^(smi)$", ext):
                    ##########################################################################################################################
                    ########## Here defines the parameters for OMEGA and ROCS, if you need to change the default ones, revise below ##########
                    ##########################################################################################################################
                    cmd = "omega2 -in "+output_ftrees_fs+file+" -out "+output_omegarocs+file+".oeb -strict false -flipper true -mpi_np "+str(coresToUse)
                    os.system(cmd)
                    cmd = "rocs -query "+query3D+" -dbase "+output_omegarocs+file+".oeb -oformat sdf -prefix "+output_omegarocs+file+" -mpi_np "+str(coresToUse)+" -besthits "+str(bestHitsNo)
                    os.system(cmd)
                    
                    # Only keeps .sdf file as ROCS results, remove the rest.
                    cmd = "find "+output_omegarocs+" -type f ! -name \"*.sdf\" -exec rm -rf {} \;"
                    os.system(cmd)
                    
                    # Use 2dsim.py to add 2D Tanimoto similarity score to the FTrees search result (.sdf file)
                    for file in os.listdir(output_omegarocs):                      
                        cmd = "python /cluster/home/cuit/scripts/2dsim.py "+output_omegarocs+file+" "+output_rocs+"rocs_"+file+" "+output_omegarocs
                        os.system(cmd)
    
def processOptions():
    """Process command line options and arguments"""

    vertexutil.info("Processing options...")

    global optionsInfo
    optionsInfo = docopt(__doc__, version=__version__)

    # Set current working directory to specified directory...
    #workingDir = optionsInfo["--workingdir"]
    #if workingDir:
    #    os.chdir(workingDir)

    # Handle examples option...
    if "--examples" in optionsInfo and optionsInfo["--examples"]:
        vertexutil.info(vertexutil.getExamplesText(__doc__))
        sys.exit(0)

    validateOptions()
    
def validateOptions():
    """Validate options and their values"""
    
    minSimilarityThreshold = float(optionsInfo["--minSimilarityThreshold"])
    if not 0 < minSimilarityThreshold <= 1:
        vertexutil.error("The value specified, %s, for option \"-s --minSimilarityThreshold\" is not valid. It must be a number between [0 to 1].\n" % minSimilarityThreshold)
        
    maxNofResults = int(optionsInfo["--maxNofResults"])
    if not 0 <= maxNofResults < 100000:
        vertexutil.error("The value specified, %s, for option \"-n --maxNofResults\" is not valid. It must be a positive integer between [0 to 100,000].\n" % maxNofResults)
    
    # Make sure infile has a valid SMILES file extension...
    #infile = optionsInfo["--infile"]
    #if infile != "stdin":
    #    if not vertexutil.checkFileExt(infile, "smi can"):
    #        vertexutil.error("The input file, %s, is not a SMILES file.\n" % infile)

    #    if not os.path.exists(infile):
    #        vertexutil.error("The input file, %s, doesn't exist.\n" % infile)
    
    # Make sure library file was provided...
    libfile = optionsInfo["--libfile"]
    if not os.path.exists(libfile):
        vertexutil.error("The library file, %s, doesn't exist.\n" % libfile)
    
    ##################################################################################
    ########## Here defines the default range of CPU core to use (1 to 50). ##########
    ##################################################################################
    coresToUse = int(optionsInfo["--coresToUse"])
    if not 1 <= coresToUse <= 50:
        vertexutil.error("The value specified, %s, for option \"-s --coresToUse\" is not valid. It must be a positive integer between [1 to 50].\n" % coresToUse)
    
    # Make sure query file with 3D coordinates has a valid file extension...
    query3D = optionsInfo["--query3D"]
    if (query3D):
        ext = query3D.split('.', 1)[1]
        if not re.match("^(oeb|oeb.gz|sdf|mol|sdf.gz|mol.gz|mol2|mol2.gz|pdb|ent|pdb.gz|ent.gz|mmod|mmod.gz)$", ext):
            vertexutil.error("The 3D query file, %s, is not a OEBinary/SDF/MOL2/PDB/MacroModel file, it has to contain 3D coordinates.\n" % query3D)
            
    # Make sure number of hits for ROCS is valid...
    bestHitsNo = int(optionsInfo["--bestHitsNo"])
    if not 0 <= bestHitsNo:
        vertexutil.error("The value specified, %s, for option \"-s --bestHitsNo\" is not valid. It must be a positive integer.\n" % bestHitsNo)
        
if __name__ == "__main__":
    main()