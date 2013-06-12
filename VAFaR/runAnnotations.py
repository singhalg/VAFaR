#-------------------------------------------------------------------------------
# Name:        VariantAnnotation.py
# Purpose:     This script annotates the variants present in the vcf file using snpEff and SnpSift.
#
# Author:      Gaurav Singhal
# Contact:     gsinghal@wustl.edu
#
# Created:     13/09/2012
# Copyright:   (c) Gaurav Singhal 2012
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------
#!/usr/bin/python

"""


"""

from subprocess import Popen, PIPE, STDOUT
import sys
import pickle
import optparse

"""
@refGenome: reference genome database version of the snpEff.
More details can be found here
<http://snpeff.sourceforge.net/download.html#databases>

"""

def snpEffSnpSift(cases, controls, refGenome='GRCh37.69', SnpSift_ref_db="./dbNSFP2.0b3.txt",  verbose=True):
    print "This script is going to invoke snpEff and SnpSift, and will require \
    4 GB system memory."
    all_files = cases + controls

    for eachFile in all_files:


        vac_file = eachFile
        outFileName = eachFile[:-4] + '_snpEff.vcf'
        summaryFile = eachFile[:-4] + '_summary.html'

        snpEffCmd = 'java -Xmx4g -jar ./snpEff.jar eff  -s ' + summaryFile + '  ' + refGenome + ' ' + legal_vcf + ' > ' + outFileName
##        cmd = 'java -jar SnpSift.jar dbnsfp -v dbNSFP2.0b3.txt ../' + fileName    + ' > '     + outFileName
        if verbose:
            print snpEffCmd
        job1 = Popen(snpEffCmd,  shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
##        print snpEffCmd
        job1.wait()


##        print snpEffCmd

        snpSiftOutput = outFileName[:-4] + '_SnpSift.vcf'

        SIFTcmd = 'java -jar ./SnpSift.jar dbnsfp -v  ' + SnpSift_ref_db + '  ' + outFileName    + '  >  '     + snpSiftOutput
        if verbose:
            print SIFTcmd
        job2 = Popen(SIFTcmd,  shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        job2.wait()


        gwasOutput = snpSiftOutput[:-4] + '_GWAS.vcf'

        GWAScmd = 'java -jar ./SnpSift.jar gwasCat ./gwascatalog.txt ' + snpSiftOutput    + '  >  '     + gwasOutput
        if verbose:
            print GWAScmd
        job3 = Popen(GWAScmd,  shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        job3.wait()

    if verbose:
        print "vcf annotation complete. Now intiating variant reading and filtration"


##def make_legal_vcf(illegal_vcf):
##    data = open(illegal_vcf, 'rU').readlines()
##    fhout_name = illegal_vcf[:-4] + '_legal.vcf'
##    fhout = open(fhout_name, 'w')
##    for line in data:
##        if line[0] == '#':
##            fhout.write(line)
##        else:
##            outline = 'chr' + line
##            fhout.write(outline)
##
##    fhout.close()
##    return fhout_name


def main():
    desc="""This script is a wrapper around snpEff and SnpSift program commonly
    used for annotating variants (SNVs).
    """
##
##    parser = optparse.OptionParser(description=desc)
##    parser.add_option('-snpEff', help='This option specifies that snpEff annotations will be used', dest='snpEff', dafault=True)
##    parser.add_option('-snpEff_ref_genome', help='This option specifies which reference genome is used for annotation. eg, for humans, the current build is GRCh37.69', dest='snpEff_ref_genome')
##    parser.add_option('-SnpSift', help='This option specifies that SnpSift annotations will be used', dest='SnpSift', default=True)
##    parser.add_option('-SnpSift_ref_db', help='This option specifies the reference database of annotations used by SnpSift', dest='SnpSift_ref_db')
##    parser.add_option('-GWAS', help='This option specifies that SnpSift annotations will be used', dest='GWAS' , default=True)
##    (opts, args) = parser.parse_args()
##    if opts.snpEff==False and opts.SnpSift==False and opts.GWAS==False:
##        print "All options are set to False, hece no annotations performed. Check Usage with script_name.py -h "
##    if opts.snpEff==True and opts.snpEff_ref_genome==None:
##        print "Reference Genome needs to be downloaded and specified. Check documentatation of snpEff"
##    if opts.SnpSift==True and opts.SnpSift_ref_db==None:
##        print "Reference Genome needs to be downloaded and specified. Check documentatation of SnpSift"

    snpEffSnpSift(cases, controls)



if __name__ == '__main__':
    main()
