#!/usr/bin/python

"""
metagenomics Pan Genome Analysis pipeline (metaPGAP) tool for generating core genome from multiple bacterial strains.
"""

__author__ = 'Mitul Patel'
__copyright__ = "Copyright 2016-17"
__license__ = "GPL"
__version__ = "v1.0"
__maintainer__ = "Mitul Patel"
__email__ = "mitul428@gmail.com"
__status__ = "Production"

# metaPGAP ##################################################################
#
# Author: Mitul Patel
#
# LICENSE ######################################################################
#
#    Copyright (C) 2016 Mitul Patel
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

################################################################################
# Imports

bExitForImportFailed=0
try:
    import sys
    from time import strftime
    import argparse
    from optparse import OptionParser, OptionGroup
    import copy
    import subprocess
    import os
    import shutil
    import time
    import fnmatch
    import platform
    import urllib
    from shutil import copyfile
    import zipfile
    import glob
    import csv
    from collections import OrderedDict
except Exception, e:
    print 'Basic imports failed!'
    print e
    bExitForImportFailed=1

################################################################################
# Color messages

class Highlighter:
    def __init__(self):
        self._msgTypes={'INF':'\033[0m',
                'IMP':'\033[1;32m',
                'DEV':'\033[1;34m',
                'ERR':'\033[1;31m',
                'WRN':'\033[1;33m'}
        self._reset='\033[0m'
        self._default='INF'
    def ColorMsg(self,msg,msgLevel='INF'):
        try:
            s=self._msgTypes[msgLevel]+msg+self._reset
        except:s=s=self._msgTypes[self._default]+msg+self._reset
        return s

def ColorOutput(msg,msgLevel='INF'):
    o=Highlighter()
    return o.ColorMsg(msg,msgLevel)

################################################################################
# Options

def getOptions():
    '''Retrieve the options passed from the command line'''

    usage = "usage: python metaPAGAP.py"
    version="metaPGAP "+__version__
    description=("metagenomics Pan Genome Analysis pipeline (metaPGAP) tool for generating core genome from multiple bacterial strains."+
                 "For bug reports, suggestions or questions mail to Mitul Patel: mitul428@gmail.com")
    parser = OptionParser(usage,version=version,description=description)
        # Parse the options
    return parser.parse_args()

def CheckRequirements():
    #debug
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Checking software requirements\n')

    # Check for BioPython
    try:import Bio
    except:
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' ERROR: BioPython is missing!\n','ERR'))
        return 1
    '''
    # Check for Prokka, get_homologues, Mafft, AMAS, RaXML
    lex=['prokka', 'get_homologues', 'mafft', 'AMAS', 'RaXML']
    for ex in lex:
        p = subprocess.Popen('which '+str(ex),shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        out=p.communicate()

        if out[0]=='':
            sys.stderr.write(strftime("%H:%M:%S")+
                  ColorOutput(' ERROR: '+ex+' executable is missing or not reachable from this location!\n','ERR'))
            return 1
    '''
    return 0

def PrintRequirements():
    sys.stderr.write('metaPGAP requirements:\n')
    #sys.stderr.write('\tInputs: -c contigs -r reference(s)\n')
    sys.stderr.write('\tSoftwares:\n')
    sys.stderr.write('\t\tPython (python.org)\n')
    sys.stderr.write('\t\tBioPython (1.5 or higher) with NumPy (biopython.org / numpy.scipy.org)\n')
    sys.stderr.write('\t\tProkka (https://github.com/tseemann/prokka)\n')
    sys.stderr.write('\t\tget_homologues (https://github.com/eead-csic-compbio/get_homologues)\n')
    sys.stderr.write('\t\tMafft (http://mafft.cbrc.jp/alignment/software/)\n')
    sys.stderr.write('\t\tPerl (Perl.org)\n')
    sys.stderr.write('\t\tAMAS (https://github.com/marekborowiec/AMAS)\n')
    sys.stderr.write('\t\tRaXML (https://github.com/stamatak/standard-RAxML)\n')
    sys.stderr.write('\t\twget (python.org)\n')
    sys.stderr.write('\tAll software should be accessible from the command line\n')
    sys.stderr.write('\t\tAdd the executables path to the PATH environmental variable\n')
    sys.stderr.write('\t\tor create a symbolic link in usr/bin or /usr/local/bin '+
                    '(as root: \"ln -s /usr/bin/EXECUTABLE /PATH/TO/EXECUTABLE\")\n')

def getmetaPGAPpath():
    current_path = os.path.abspath(sys.argv[0])
    current_path = os.path.abspath(os.path.join(current_path, ".."))
    tmp_path = current_path
    metaPGAP_path = current_path
    last_tmp_path = tmp_path
    it_count = 0
    while len(tmp_path) > 1 and it_count < 20:
        it_count += 1
        current_dir = os.path.basename(tmp_path)
        current_path = os.path.abspath(os.path.join(tmp_path, ".."))
        sys.stdout.flush()
        
        if current_dir.startswith("metaPGAP"):
            metaPGAP_path = tmp_path + "/"
            break
        elif os.path.exists(current_path + "/metaPGAP/"):
            metaPGAP_path = current_path + "/metaPGAP/"
            break
        else:
            tmp_path = current_path
        if last_tmp_path == tmp_path or it_count == 20:
            print("Cannot determine metaPGAP path")
            exit()
    return metaPGAP_path;

def directories():
    step = 1
    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Creating directories .....\n','IMP'))

    dirs = ["data", "prokka_out", "pan_genome", "phylogeny"]

    for d in dirs:
        if os.path.isdir(d):
            pass
        else:
            os.mkdir(d)
    return step

def fetchDATA():
    step = 2
    os.chdir(os.path.abspath("data/"))
    file_path=os.path.abspath(os.curdir)

    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Downloading input files .....\n','IMP'))

    t0_time = time.time()
    if len(os.listdir(file_path)) > 0:
        unzip("master.zip", file_path)
        #print "already exist"
        pass
        #if (urllib.urlretrieve("https://jjcloud.app.box.com/s/rw0v7r6thtv6efsfs2ekp8gontaw7quf", "test-contigs.zip")):
    elif (urllib.urlretrieve("https://github.com/mitul-patel/data/archive/master.zip", "master.zip")): 
        unzip("master.zip", file_path)
       #urllib.urlretrieve('https://jjcloud.app.box.com/s/rw0v7r6thtv6efsfs2ekp8gontaw7quf', 'test-contigs.zip')
    t1_time = time.time()
    print(t1_time - t0_time, "Download")
    #urllib.urlretrieve("https://github.com/marekborowiec/AMAS/archive/master.zip", "master.zip")
    os.chdir("../")
    return step

def unzip(zipFilePath, destDir):

    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Extracting files for zip folder .....\n','IMP'))

    zfile = zipfile.ZipFile(zipFilePath)
    for name in zfile.namelist():
        (dirName, fileName) = os.path.split(name)
        # Check if the directory exisits
        newDir = destDir + '/' + dirName
        if not os.path.exists(newDir):
            os.mkdir(newDir)
        if not fileName == '':
            # file
            fd = open(destDir + '/' + name, 'wb')
            fd.write(zfile.read(name))
            fd.close()
    zfile.close()

def prokka(platform):
    step = 3

    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Annotating genome with Prokka .....\n','IMP'))

    if platform.lower() == "darwin":
        prokka_path=os.path.abspath("binaries/prokka/darwin/bin/")
    elif platform.lower() == "linux":
        prokka_path=os.path.abspath("binaries/prokka/linux/bin/")
    
    data_path=os.path.abspath("data/data-master/")
    out_path=os.path.abspath("prokka_out/")
    files = os.listdir(data_path)
    
    params = ["--kingdom", "--outdir", "--prefix", "--locustag", ""]
    args = ["--kingdom Bacteria"]
    #passed_params = createPara(params, parameters)
    
    t0_time = time.time()

    if len(os.listdir(out_path)) == len(files):
        pass
    else:     
        for i in files:
            sample = i.split('-')[0]                
            cmd = '%s/prokka --kingdom Bacteria --outdir %s/prokka_%s --prefix %s --locustag %s %s/%s' % (prokka_path, out_path, sample, sample, sample, data_path, i)
            if os.system(cmd):exit()

    t1_time = time.time()
    print(t1_time - t0_time, "Prokka")

    return step

def copyGBK():
    step = 4
    prokka_res=os.path.abspath("prokka_out/")
    
    if os.path.isdir("prokka_allGBK"):
        pass
    else:
        os.mkdir("prokka_allGBK")
        #os.chdir("../")
        #gbkPath = os.path.abspath("prokka_allGBK/")
        for root, dirs, files in os.walk(prokka_res):
            for file in files:
                if file.endswith('.gbk'):
                    shutil.copy2(os.path.join(root, file), 'prokka_allGBK')    
    return step

def get_panGenome():
    step = 5

    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Pan genome analysis. Creating core genome .....\n','IMP'))

    get_homologues_path=os.path.abspath("binaries/get_homologues-macosx-20161027/")
    data_path=os.path.abspath("prokka_allGBK/")
    out_path=os.path.abspath("prokka_allGBK_homologues/")
    files = os.listdir(data_path)
    
    BDBH = 'perl %s/get_homologues.pl -d %s -n 10' % (get_homologues_path, data_path)
    OMCL = 'perl %s/get_homologues.pl -d %s -M -n 10' % (get_homologues_path, data_path)
    COG = 'perl %s/get_homologues.pl -d %s -G -n 10' % (get_homologues_path, data_path)
    
    #folders = ["algBDBH", "algCOG", "algOMCL"]

    t0_time = time.time()
    if os.system(BDBH):exit()
    if os.system(OMCL):exit()
    if os.system(COG):exit()
    
    thelist = os.listdir(out_path)
    for fo in thelist:
        if fo.find("alltaxa_alg") != -1:
            samp = fo.split('.')[0]
            intersec = 'perl %s/compare_clusters.pl -o %s/intersection -d %s/%s,%s/%s,%s/%s -n' % (get_homologues_path, out_path, out_path, samp, out_path, samp, out_path, samp)
            if os.system(intersec):exit()

    t1_time = time.time()
    print(t1_time - t0_time, "Pan Genome")

    intersec_files_path=os.path.abspath("prokka_allGBK_homologues/intersection")
    renameHeaders(intersec_files_path)

    return step

def renameHeaders(intersecFiles):
    intersec_path=intersecFiles
    if os.path.isdir("tmp"):
        shutil.rmtree('tmp')
        os.mkdir("tmp")
        if os.path.isdir(intersec_path):
            for root, dirs, files in os.walk(intersec_path): 
                for file in files:               
                    if file.endswith('.fna'):
                        shutil.copy2(os.path.join(root, file), 'tmp')
        else:
            exit()
    else:
        os.mkdir("tmp")
        if os.path.isdir(intersec_path):
            for root, dirs, files in os.walk(intersec_path): 
                for file in files:               
                    if file.endswith('.fna'):
                        shutil.copy2(os.path.join(root, file), 'tmp')

    tmp_path=os.path.abspath("tmp/")
    pan_path=os.path.abspath("pan_genome/")
    for root, dirs, files in os.walk(tmp_path):
            for file in files:
                if file.endswith('.fna'):
                    fl = os.path.join(root, file)
                    cmd="sed 's/^>ID:/>/' %s | sed 's/\_.*//' > %s/%s" % (fl,pan_path,file)
                    if os.system(cmd):exit()
    shutil.rmtree('tmp')


def mafft(platform):
    step = 6
    
    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Multiple sequence alignemnt of core genes .....\n','IMP'))

    if platform.lower() == "darwin":
        mafft_path=os.path.abspath("binaries/mafft/darwin/")
    elif platform.lower() == "linux":
        mafft_path=os.path.abspath("binaries/mafft/linux/bin/")

    data_path=os.path.abspath("pan_genome/")
    out_path=os.path.abspath("phylogeny/")
    
    t0_time = time.time()
    for root, dirs, files in os.walk(data_path):
        for file in files:
            if file.endswith('.fna'):
                fl = os.path.join(root, file)
                fname=file.split('.')[0] 
                cmd= '%s/mafft --thread 10 --auto --quiet %s > %s/%s.aln' % (mafft_path,fl,out_path,fname)
                if os.system(cmd):exit()
    
    t1_time = time.time()
    print(t1_time - t0_time, "Mafft")
    
    return step

def AMAS():
    step = 7

    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Concatenating alignments .....\n','IMP'))

    AMAS_path=os.path.abspath("binaries/AMAS-master/amas/")
    data_path=os.path.abspath("phylogeny/")
    out_path=os.path.abspath("phylogeny/")
        
    os.chdir(out_path)
    t0_time = time.time()
    cmd="python %s/AMAS.py concat -f fasta -d aa -i %s/*aln --out-format phylip --part-format raxml" % (AMAS_path,data_path)
    if os.system(cmd):exit()
    t1_time = time.time()
    print(t1_time - t0_time, "AMAS")
    os.chdir("../")

    return step

def RaXML(platform):
    step = 8

    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Phylogenetic analysis of core genome using RaXML .....\n','IMP'))

    if platform.lower() == "darwin":
        RaXML_path=os.path.abspath("binaries/standard-RAxML-master/darwin/")
    elif platform.lower() == "linux":
        RaXML_path=os.path.abspath("binaries/standard-RAxML-master/linux/")

    data_path=os.path.abspath("phylogeny/")
    out_path=os.path.abspath("phylogeny/")
        
    os.chdir(out_path)
    t0_time = time.time()
    cmd="%s/raxmlHPC -m PROTGAMMABLOSUM62 -p 12345 -x 12345 -s %s/concatenated.out -n CORE -f a -N 50" % (RaXML_path,data_path)
    if os.system(cmd):exit()
    t1_time = time.time()
    print(t1_time - t0_time, "RaXML")
    os.chdir("../")

    return step

def nwTree(platform):
    step = 9
    
    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Building phylogenetic tree using Newick tools .....\n','IMP'))

    if platform.lower() == "darwin":
        nwTree_path=os.path.abspath("binaries/nwTree/darwin/")
    elif platform.lower() == "linux":
        neTree_path=os.path.abspath("binaries/nwTree/linux/bin/")

    data_path=os.path.abspath("phylogeny/")
    out_path=os.path.abspath("phylogeny/")
    
    CORE_file="RAxML_bestTree.CORE"
    t0_time = time.time()
    for root, dirs, files in os.walk(data_path):
        if CORE_file in files:
            fl = os.path.join(root, CORE_file)
            fname=CORE_file.split('.')[0] 
            cmd= '%s/nw_display -s -S -b opacity:0 %s > %s/%s.svg' % (nwTree_path,fl,out_path,fname)
            #pdf= 'inkscape -f %s/%s.svg -A %s.pdf'% (fname,out_path,fname)
            if os.system(cmd):exit()
            #if os.system(pdf):exit()
      
    t1_time = time.time()
    print(t1_time - t0_time, "nwTree")
    
    return step

sys.path += ["binaries/", "../", "../../", "../../../"]
sys.path += ["data/", "../", "../../", "../../../"]

#sys.path += "binaries/prokka/darwin/bin/"
#sys.path += "binaries/prokka/darwin/bin/"
#sys.path += "binaries/mafft/darwin/bin/"
#sys.path += "binaries/mafft/linux/bin/"
#sys.path += "binaries/standard-RAxML-master/darwin/"
#sys.path += "binaries/standard-RAxML-master/linux/"
#sys.path += "binaries/nwTree/darwin/"
#sys.path += "binaries/nwTree/linux/bin/"

run_path = os.path.abspath(os.curdir)
metaPGAP_path = getmetaPGAPpath()
os.chdir(metaPGAP_path)

def main():
    os.chdir(metaPGAP_path)
    plat=platform.system()
    step = 1

    if bExitForImportFailed:
        pass
    elif CheckRequirements():
        PrintRequirements()
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
                            ColorOutput(' Stopping metaPGAP\n','WRN'))
        mylog.WriteLog('INF','Stopping metaPGAP')
        sys.exit(1)
    else:
        # Message
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
            ColorOutput(' Starting metaPGAP\n','IMP'))
        for step in range(1,10):
            if step == 1:
                stepCOM=directories()
            elif step == 2:
                stepCOM=fetchDATA()
            elif step == 3:
                stepCOM=prokka(plat)      
            elif step == 4:
                stepCOM=copyGBK()      
            elif step == 5:
                stepCOM=get_panGenome()      
            elif step == 6:
                stepCOM=mafft(plat)
            elif step == 7:
                stepCOM=AMAS()  
            elif step == 8:
                stepCOM=RaXML(plat)
            elif step == 9:
                stepCOM=nwTree(plat)
    # Message
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
        ColorOutput(' Finished metaPGAP\n','IMP'))
    sys.exit(1)

if __name__ == '__main__':
    main()
