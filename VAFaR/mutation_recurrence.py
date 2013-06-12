#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     07/06/2013
# Copyright:   (c) Gaurav Singhal 2013
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------

import sys
from sets import Set


def mutationRecurrence(vcf_file_list, recurrence_level):

    mutation_patient_dict = {}

    for each in vcf_file_list:

        pickleFile = each+'_snpEff_SnpSift_GWAS_del.pkl'

        mutation_data = pickle.load(open(pickleFile, 'rU'))

        for eachMutation in mutation_data.keys():
            if eachMutation in mutation_patient_dict:
                count = mutation_patient_dict[eachMutation][0]
                count+=1
                mutation_patient_dict[eachMutation] = [count, mutation_data[eachMutation]]
            else:
                mutation_patient_dict[eachMutation] = [1, mutation_data[eachMutation]]




    recurrentVariants = returnRecurrentVariants(mutation_patient_dict, recurrence_level)

##    pickle.dump(mutation_patient_dict, open('mutation_patient_control_dict.pkl', 'w'))
##
##    pickle.dump(recurrentVariants, open('recurrent_mutations_control.pkl', 'w'))

    return recurrentVariants




def returnRecurrentVariants(mutation_patient_dict, recurrence_level):

    variant_set = Set()
    for eachLoc in mutation_patient_dict:
        ref_allele = mutation_patient_dict[eachLoc][1]['REF']
        alt_allele = mutation_patient_dict[eachLoc][1]['ALT']
        if (mutation_patient_dict[eachLoc][0] > recurrence_level) and (len(ref_allele)<=1):
            for eachAlt in alt_allele:

                key = eachLoc + '#' + ref_allele + '#' + eachAlt
                variant_set.add(key)


    return variant_set

def main():
    mutationRecurrence(vcf_file_list, recurrence_level)

if __name__ == '__main__':
    main()
