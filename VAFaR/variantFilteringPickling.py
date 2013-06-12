#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     14/09/2012
# Copyright:   (c) Gaurav Singhal 2012
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------

import pickle
import sys
import numpy as np
from sets import Set
from mutation_recurrence import mutationRecurrence, returnRecurrentVariants
'''
Modufying the script to create 3 pickles. One which contains all the variants, one which contains only deleterious variants, and one
which contains only rare delerious variants. Also, I collect all the variant detection stats and report them in a csv file.


vcf2pickle.py was written to convert files from vcf format like 'PT_10_brain_1200360.vcf' into 2 pickles, one containing all variants and other containing only
deleterious variants. This vcf file is the output of snpEff

Here I am modifying this script. This new script (vcf2PickleV2.py) takes in a vcf file which is the output of snpEff and SnpSift (has predictions from both)
and saves 2 pickles, one containing all variants and the other containing only rare deleterious (MISSENSE and NONSENSE) variants. Thus, the second pickle
does not have any variant that appears in internal controls or external set of common variants (dbSNP, 1000 Genomes, EVS etc)

Additional tags:

'SIFT', 'GERP', 'PolyPhen2', 'score', 'logOdds'

'''

def VCFtoPickle(cases, controls, dbSNP_common, dbSNP_common_clinical, recurrence_level):
    cases_internal_common = internal_recurrence(cases, dbSNP_common, dbSNP_common_clinical, recurrene_level, 'cases_pre_filtered')
    controls_internal_common = internal_recurrence(controls, dbSNP_common, dbSNP_common_clinical, recurrene_level, 'controls_pre_filtered')

    filter_pickle_VCF(cases, dbSNP, dbSNP_common_clinical, cases_internal_common, list_name_cases)
    filter_pickle_VCF(controls, dbSNP, dbSNP_common_clinical, controls_internal_common, list_name_controls)




"""

This method is
@dbSNP_
@dbSNP_common_clinical

"""

def internal_recurrence(list_of_vcf_filenames, dbSNP_common, dbSNP_common_clinical, recurrene_level, list_name):

    interal_control_set = Set([])

    filter_pickle_VCF(list_of_vcf_filenames, dbSNP_common, dbSNP_common_clinical, interal_control_set, list_name)
    internal_common = mutationRecurrence(list_of_vcf_filenames, recurrene_level)
    return internal_common





"""
This method takes in a list of filenames and saves 3 pickles for each filename.
One pickle containing all variants, one containing only non-synonymous variants,
 and one containing only non-synonymous germline/somatic variants.
"""
def filter_pickle_VCF(list_of_vcf_filenames, dbSNP_common, dbSNP_common_clinical, internal_control, list_name):

    fhout_file = list_name +'.csv'
    fhout = open(fhout_file, 'w')
    outline = 'TISSUE_SOURCE, ALL_VARIANTS, DELETERIOUS_VARIANTS, RARE_DEL_VARIANTS' + '\n'
    fhout.write(outline)

    dbSNP_all = open(dbSNP_common, 'rU')
    dbSNP = pickle.load(dbSNP_all)
    dbSNP_pickle



    dbSNP_clinical_fh = open(dbSNP_common_clinical, 'rU')
    dbSNP_common_clinical = pickle.load(dbSNP_clinical_fh)


    dbSNP_common_clinical_variants = returnVariants(dbSNP_common_clinical)
    dbSNP_common_variants = returnVariants(dbSNP)

    dbSNP_common_nonClinical = dbSNP_common_variants - dbSNP_common_clinical_variants

    all_common = dbSNP_common_nonClinical | internal_control

    for eachVCF in list_of_vcf_filenames:




        variantVCF, deleteriousVCF, deleteriousRareVCF, variantCounts = readVCF(eachVCF, all_common)

        allVCFPickleName = eachVCF[:-4]+'.pkl'
        fhAllVCF = open(allVCFPickleName, 'w')
        pickle.dump(variantVCF, fhAllVCF)
        fhAllVCF.close()

        deleteriousVCFpickleName = vcfFile[:-4]+'_del.pkl'
        fhDelVCF = open(deleteriousVCFpickleName, 'w')
        pickle.dump(deleteriousVCF, fhDelVCF)
        fhDelVCF.close()

        deleteriousRareVCFPickle = vcfFile[:-4]+'_rare_del.pkl'
        fhDelRare = open(deleteriousRareVCFPickle, 'w')
        pickle.dump(deleteriousRareVCF, fhDelRare)


        print vcfFile, ' processed !'
        print variantCounts
        outline = tissue_identifier + ',' + str(variantCounts[0]) + ',' + str(variantCounts[1]) + ',' + str(variantCounts[2]) + '\n'

        fhout.write(outline)
        del variantVCF, deleteriousVCF  , deleteriousRareVCF # finally purging the objects from memory
    fhout.close()









def returnVariants(dbSNP_dict):
    variant_set = Set()
    for eachLoc in dbSNP_dict:
        ref_allele = dbSNP_dict[eachLoc]['REF']
        if len(ref_allele)>1:
            pass

        else:
            alt_allele =dbSNP_dict[eachLoc]['ALT']
            for eachAlt in alt_allele:

                key = eachLoc + '#' + ref_allele + '#' + eachAlt
                variant_set.add(key)
    return variant_set


def readVCF(fileName, commonSnpSet):
    fh = open(fileName, 'rU')
    data = fh.readlines()
    variantsVCF = {}
    deleteriousVariantsVCF = {}
    deleteriousRareVariantsVCF = {}

##    varSet = Set([])

    for line in data:
        if line[0] != '#':
            flds = line.split('\t')
            if len(flds[7].split(';'))>2:

                processVCFline(line, variantsVCF)



    variantPos = sorted(variantsVCF.keys())
    for apos in variantPos:
##        if apos not in commonSnpSet:  # filtering out all those variants which are present in allCommon (not rare)

        info = variantsVCF[apos]
        if 'Effect_Impact' in info:  # changing from Functional_Class to Effect_Impact
            func_class = info['Effect_Impact']

            if ('MODERATE' in func_class) or ('HIGH' in func_class):

                deleteriousVariantsVCF[apos] = info

    deleteriousVariantPos = deleteriousVariantsVCF.keys()
    for eachPos in deleteriousVariantPos:

        ref_allele = deleteriousVariantsVCF[eachPos]['REF']
        alt_allele = deleteriousVariantsVCF[eachPos]['ALT']

        for eachAlt in alt_allele:
            variant = eachPos + '#' + ref_allele + '#' + eachAlt

            if variant not in commonSnpSet:
                deleteriousRareVariantsVCF[eachPos] = deleteriousVariantsVCF[eachPos]


    variantCounts = [len(variantPos), len(deleteriousVariantPos), len(deleteriousRareVariantsVCF.keys())]
    return variantsVCF, deleteriousVariantsVCF, deleteriousRareVariantsVCF, variantCounts


'''
This method has been especially edited to parse the snpEff output vcf.

'''
def processVCFline(line, variantDict):
    flds = line.split('\t')
    key = flds[0].strip() + '#' + flds[1].strip()
    value = {}
    value['REF'] = flds[3].strip()
    value['ALT'] =   flds[4].strip().split(',') # ALT is a list
    for each in flds[7].strip().split(';'):
        info = each.partition('=')
        infoVal = info[2].split(',')
        value[info[0]] = infoVal  # DP, AF, PVAL, EFF are all lists

##    variantDict[key] = value
    snpEff = value['EFF']
##    snpEff = flds[7].strip().split(';')[2][4:].split(',')
##    print snpEff

##    valueBack = variantDict[key]
####    if 'EFF' in valueBack:
##    snpEff = valueBack['EFF']
##    else:
##        print line
##        sys.exit()
##        except:
##            print line
##            print flds[7]
##            sys.exit()
##    impacts = []
    '''
'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change|
Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript |
Exon [ | ERRORS | WARNINGS ] )'

    '''

    for eachEffect in snpEff:
        flds = eachEffect.partition('(')[2].split('|')
        effect = eachEffect.partition('(')[0]

        addValuetoDict('Effect', effect, value) # value is the dictionary here, this dictionary will be added as a value to variantDict later
        addValuetoDict('Effect_Impact', flds[0], value)
        addValuetoDict('Functional_Class', flds[1],value)
##        addValuetoDict('Codon_Change',flds[2],value)
##        addValuetoDict('Amino_Acid_Change',flds[3],value)
        addValuetoDict('Amino_Acid_Length',flds[4],value)
        addValuetoDict('Gene_Name',flds[5],value)
##        addValuetoDict('Gene_BioType',flds[6],value)
        addValuetoDict('Coding',flds[7],value)
        addValuetoDict('Transcript',flds[8],value)
        addValuetoDict('Exon',flds[9].partition(')')[0] ,value)

#calculating the SIFT score
    scores = []

    if 'dbnsfpSIFT_score' in value:
##        try:
        siftScores = list2FloatNParray(value['dbnsfpSIFT_score'])
##        for eachSift in value['dbnsfpSIFT_score']:
##            try:
##                fSift = float(eachSift)
##                siftScores.append(fSift)
##            except:
##                pass

        if len(siftScores)>0:

            SiftNew = 1 - (np.array(siftScores)).min()
            if SiftNew >=0.95:
                value['SIFT'] = SiftNew
                scores.append(SiftNew)

            else:
                value['SIFT'] = SiftNew/2
                newSc = SiftNew/2
                scores.append(newSc) #-------

##        except:
##            print value['dbnsfpSIFT_score']
##        SiftNew = 1 - siftScores.min()
##
##
##        if SiftNew >=0.95:
##            value['SIFT'] = SiftNew
##            scores.append(SiftNew)
##
##        else:
##            value['SIFT'] = SiftNew/2
##            newSc = SiftNew/2
##            scores.append(newSc) #-------

#calculating the GERP score

    if ('dbnsfpGERP++_RS' in value) and ('dbnsfpGERP++_NR' in value):
        gerp_rs = list2FloatNParray(value['dbnsfpGERP++_RS'])
##        try:
##            gerp_rs = np.array([float(n) for n  in value['dbnsfpGERP++_RS']])
##            gerp_rs_max = gerp_rs.max()
##        except :
##            print line

        gerp_nr = list2FloatNParray(value['dbnsfpGERP++_NR'])
        if (len(gerp_rs)>0) and (len(gerp_nr)>0):
            gerp_nr_max = gerp_nr.max()
            gerp_rs_max = gerp_rs.max()
            if gerp_rs_max>0:
                gerpScore = gerp_rs_max/gerp_nr_max
            else:
                gerpScore = 0
            if gerpScore >=1:
                value['GERP'] = 1
                scores.append(1)
            else:
                value['GERP'] = gerpScore
                scores.append(gerpScore)

#calculating polyPhen2 score

    if 'dbnsfpPolyphen2_HVAR_pred' in value:
        polyphen2 = value['dbnsfpPolyphen2_HVAR_pred']
        polyPhenScores = []
        for eachPred in polyphen2:
            if eachPred =='D':
                polyPhenScores.append(0.9)
            elif eachPred =='P':
                polyPhenScores.append(0.5)
            elif eachPred =='B':
                polyPhenScores.append(0)
        if len(polyPhenScores)>0:
            ppScore = np.array(polyPhenScores).max()
            value['PolyPhen2'] = ppScore
            scores.append(ppScore)

#calculating logOdds
    if 'dbnsfp29way_logOdds' in value:
        logOdds = list2FloatNParray(value['dbnsfp29way_logOdds'])
##        logOdds = np.array([float(n) for n in value['dbnsfp29way_logOdds']]).max()
        if len(logOdds)>0:

            LOScore  = (logOdds.max())/14
            if LOScore>=1:

                value['logOdds'] = 1
                scores.append(1)

            else:
                value['logOdds'] = LOScore
                scores.append(LOScore)
    if len(scores)>0:

##        print scores

        combinedScore = (np.array(scores)).sum() / float(len(scores))
        value['score'] = combinedScore






    variantDict[key] = value

def addValuetoDict(akey, newValue, dictionary):
    if akey in dictionary:
        value = dictionary[akey]
        if len(newValue) > 0:
            value.add(newValue)
            dictionary[akey] = value
    else:
        value = Set([])
        if len(newValue) > 0:
            value.add(newValue)
            dictionary[akey] = value

'''
Converts a list of numbers into floats, and returns them as an numpy array. Also takes care if some elements are strings,
it removes the elements which are strings and just returns the float numbers.
'''

def list2FloatNParray(alist):
    floatArray = []
    for each in alist:
        try:

            fScore = float(each)
            floatArray.append(fScore)
        except:
            pass

    return np.array(floatArray)

def main():
    dbSNp_pickle = sys.argv[1]
    VCFtoPickle(dbSNp_pickle)




if __name__ == '__main__':
    main()
