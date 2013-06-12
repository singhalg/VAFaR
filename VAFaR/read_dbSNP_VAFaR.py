#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     10/09/2012
# Copyright:   (c) Gaurav Singhal 2012
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------
import sys, pickle
from sets import Set

'''
reads the vcf file having dbSNP variants, saves them to a pickle. This script needs large amount of memory, approx 10x the size of common_all.vcf.
'''

"""
The read_dbSNP method has been edited to read 'common_all.vcf' downloaded from
dbSNP website. This vcf file contains all common SNPs that are polymorphic in
at least one among the 1000 Genomes project or
any of the following handles: CSHL-Hapmap, EGP_SNPS, NHLBI-ESP, PGA-UW-FHCRC.
For more details about 'common_all.vcf', please check
 <ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00--README.txt>.
common_all.vcf.gz can be found
at <ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all.vcf.gz>

Saves a pickle of the dict object consisting of common variants derived
 from common_all.vcf

@filename: common_all.vcf, Required
@build: dbSNP build specification, Default = 135
@maf: minor allele frequency threshold, Default = 0.03
@verbose: Prints status of script if True. Deault = False
@returns: Null


"""
def read_dbSNP(fileName, build=135, maf=0.03, verbose=False):
    if verbose:
        print 'Running read_dbSNP with the following arguments: ', mfileName,
        build, maf, verbose
    fh = open(fileName, 'rU')
    gmaf = float(maf) # global minor allele freq threshold. This can be set to
                      # a float between 0 and 1. If gmaf is set to 0.03, then
                      # any variant which has gmaf>=0.03 will be saved in the
                      # dict object

    #dbSNP build id: If buildID is set to 132, then all variants in dbSNP builds
    # 132 and earlier would be included in the dict object
    buildID = int(build)
    variantsVCF = {}
    """



    """
    for line in fh:
        if line[0] != '#': #indicates that it is contains a variant in vcf file

            processVCF(line, variantsVCF, buildID, gmaf)

        else:
            pass
    fh.close()

##    variantPos = sorted(variantsVCF.keys())
##    for apos in variantPos:
##        info = variantsVCF[apos]
##        if ('GMAF' in info) and (float(info['GMAF'][0]) >=0.01):
##            commonVariantsVCF[apos] = info
##
##
##        elif 'G5A' in info:
##            commonVariantsVCF[apos] = info

##    print '# of all variants in dbSNP file = ', len(fh)

    outfileName = 'dbSNP_build_' + str(build) + '_maf_'+ str(maf)+'.pkl'
    #  eg: outfileName = dbSNP_build_135_maf_0.03.pkl
    if verbose:
        print '# of variants saved in the pickle ', outfileName, len(variantsVCF.keys())
    fhPickledbSNP = open(outfileName, 'w')
    pickle.dump(variantsVCF, fhPickledbSNP)
    fhPickledbSNP.close()
    return outfileName

"""
Parses every line of the vcf file in dbSNP vcf file and adds the variant to the
 dict if it meets the specifications of build id and maf
@line: line of vcf file
@variantDict: variant dict object. This script adds variants to this dict. Each
 variant is a key value pair, where variant location is the key and variant
 information is the pair.
@build: dbSNP Build (int)
@maf: threshold minor allele freq

"""

def processVCF(line, variantDict, build, maf):
    flds = line.split('\t')
    key = flds[0].strip() + '#' + flds[1].strip() # key = chrm#coordinate or 1$10034234  if the variant is located on chr1 and base 10034234
    value = {}
    value['REF'] = flds[3].strip()
    value['ALT'] = flds[4].strip().split(',')
    tags = Set(['G5A', 'G5'])
    for each in flds[7].strip().split(';'):
        info = each.partition('=')
        infoVal = info[2].split(',')
        if info[0] in tags:
            value[info[0]] = 1
        elif info[0] =='GMAF':
            value[info[0]] = infoVal
        elif info[0] == 'dbSNPBuildID':
            value[info[0]] = infoVal
##

    if ('GMAF' in value) and (float(value['GMAF'][0]) >=maf) and ('dbSNPBuildID' in value) and (int(value['dbSNPBuildID'][0]) <= build) :
        variantDict[key] = value


    del value # try to minimize the memory footprint for our objects




def read_clin_sign_snps(fileName):
    fh = open(fileName, 'rU')

    variantsVCF = {}
##    commonVariantsVCF = {}

    for line in fh:
        if line[0] != '#':

            processVCFClin(line, variantsVCF)

        else:
            pass
    fh.close()

    fhout = fileName + '.pkl'
    fhPickledbSNP = open(fhout, 'w')
    pickle.dump(variantsVCF, fhPickledbSNP)
    fhPickledbSNP.close()
    return fhout



def processVCFClin(line, variantDict):
    flds = line.split('\t')
    key = flds[0].strip() + '#' + flds[1].strip()
    value = {}
    value['REF'] = flds[3].strip()
    value['ALT'] = flds[4].strip().split(',')
    tags = Set(['G5A', 'G5', 'PM', 'TPA', 'PMC', 'NSF', 'NSM', 'NSN', 'ASS', 'DSS', 'U3', 'U5', 'INT', 'R3', 'R5' ])

    fields = Set(['GMAF', 'CLNSIG', 'CLNDBN', 'GENEINFO'])

    for each in flds[7].strip().split(';'):
        info = each.partition('=')

        if info[0] in tags:
            value[info[0]] = 1
        elif info[0] in fields:
            value[info[0]] = info[2]



    variantDict[key] = value


    del value # making sure that if the value isn't stored inside variantDict, it is removed from memory. The dbSNP file is very big and we want to
        # try to minimize the memory footprint for our objects/variables






def main():

    read_dbSNP('common_all.vcf', 'dbSNP_common_build135_maf3.pkl', '134', 0.01 )
##    read_clin_sign_snps('common_and_clinical_20130226.vcf')

if __name__ == '__main__':
    main()
