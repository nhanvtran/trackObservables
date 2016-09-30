import os,sys
from optparse import OptionParser
import math

pwd=os.environ['PWD']

# get options
parser = OptionParser(description='Submit condor jobs for anaSubstructure calls.')
parser.add_option('--indir',      action='store', dest='indir',      default="",    help='Input LHE sample directory')
parser.add_option('--outdir',     action='store', dest='outdir',     default="./",  help='Location of output directory')
parser.add_option('--outdirname', action='store', dest='outdirname', default="anafull",  help='Location of output directory')
parser.add_option('--evPerJob',   action='store', dest='evPerJob',   default=10000, help='Number of events to run over in each job')
parser.add_option('--maxEvents',  action='store', dest='maxEvents',  default=50000, help='Number of events in each input sample')
parser.add_option('--anaSubLoc',  action='store', dest='anaSubLoc',  default="%s/anaSubstructure"%(pwd),            help='Location of anaSubstructure')
parser.add_option('--fastJetLoc', action='store', dest='fastJetLoc', default="%s/fastjet/fastjet-install/"%(pwd),   help='Location of fastjet-install')

(options, args) = parser.parse_args()
cmssw_base = os.environ['CMSSW_BASE']

# check that input directory is specified
if options.indir == "":
    print "ERROR: no --indir specified. Exiting..."
    quit()

# check that output directory exists
options.outdir = os.path.abspath(options.outdir)
if not os.path.exists(options.outdir):
    os.system('mkdir -p ' + options.outdir)

# get ls of samples
files = os.listdir("/eos/uscms/%s"%options.indir)

from math import ceil
nFilesPerLHE=int(ceil(float(options.maxEvents)/float(options.evPerJob)))

# prepare all configs
for i,j in [(i,j) for i in range(len(files)) for j in range(nFilesPerLHE)]:
    if files[i].split('.')[-1] != "lhe" : continue
    if "WW" not in files[i].split('.')[0] : continue
    current_name = options.outdir + '/' + files[i].split('.')[0] + "_%i"%j
    current_conf = open(current_name + '.condor','w')
    current_shel = open(current_name + '.sh','w')

    minEv = int(options.evPerJob)*j
    maxEv = int(options.evPerJob)*(j+1)-1

    # write condor config
    conf_tmpl = open('../condor/CondorConf.tmpl.condor')
    for line in conf_tmpl:
        if 'OUTPUT_PATH' in line: line = line.replace('OUTPUT_PATH', options.outdir)
        if 'OUTDIRFOLD'  in line: line = line.replace('OUTDIRFOLD',  options.outdirname)
        if 'INDIR'       in line: line = line.replace('INDIR',       "root://cmseos.fnal.gov///%s"%options.indir)
        if 'FILE'        in line: line = line.replace('FILE',        files[i].split('.')[0])
        if 'NAME'        in line: line = line.replace('NAME',        files[i].split('.')[0] + "_%i"%j)
        if 'CMSSWBASE'   in line: line = line.replace('CMSSWBASE',   cmssw_base)
        if 'ANASUBLOC'   in line: line = line.replace('ANASUBLOC',   options.anaSubLoc)
        if 'FASTJETLOC'  in line: line = line.replace('FASTJETLOC',  options.fastJetLoc)
        if 'MINEV'       in line: line = line.replace('MINEV',       minEv)
        if 'MAXEV'       in line: line = line.replace('MAXEV',       maxEv)
        if 'TAG'         in line: line = line.replace('TAG',         "%i"%j)

        current_conf.write(line)

    conf_tmpl.close()

    # write shell script
    shel_tmpl = open('../condor/CondorShel.tmpl.sh')
    for line in shel_tmpl:
        if 'OUTPUT_PATH' in line: line = line.replace('OUTPUT_PATH', options.outdir)
        if 'OUTDIRFOLD'  in line: line = line.replace('OUTDIRFOLD',  options.outdirname)
        if 'INDIR'       in line: line = line.replace('INDIR',       "root://cmseos.fnal.gov///%s"%options.indir)
        if 'FILE'        in line: line = line.replace('FILE',        files[i].split('.')[0])
        if 'NAME'        in line: line = line.replace('NAME',        files[i].split('.')[0] + "_%i"%j)
        if 'CMSSWBASE'   in line: line = line.replace('CMSSWBASE',   cmssw_base)
        if 'ANASUBLOC'   in line: line = line.replace('ANASUBLOC',   options.anaSubLoc)
        if 'FASTJETLOC'  in line: line = line.replace('FASTJETLOC',  options.fastJetLoc)
        if 'MINEV'       in line: line = line.replace('MINEV',       "%i"%minEv)
        if 'MAXEV'       in line: line = line.replace('MAXEV',       "%i"%maxEv)
        if 'TAG'         in line: line = line.replace('TAG',         "%i"%j)

        current_shel.write(line)

    shel_tmpl.close()

    current_conf.close()
    current_shel.close()

# run jobs
for i,j in [(i,j) for i in range(len(files)) for j in range(nFilesPerLHE)]:
    if files[i].split('.')[-1] != "lhe" : continue
    if "WW" not in files[i].split('.')[0] : continue
    current_name = files[i].split('.')[0] + "_%i"%j
    os.chdir(options.outdir + '/')
    os.system('condor_submit ' + current_name + '.condor')

