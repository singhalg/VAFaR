#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     14/05/2013
# Copyright:   (c) Gaurav Singhal 2013
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------

import sys, pickle
from sets import Set

def compound_heterozygosity():
    gene_patient_dict = pickle.load(open('gene_patient_data_v_rare_deleterious.pkl', 'rU'))
    outfile = open('gene_patient_compound_heterozygous.csv', 'w')
    genes = []

    patient_gene_dict = {}

    for eachGene in gene_patient_dict.keys():
        patient_dict = gene_patient_dict[eachGene]
        for eachPatient in patient_dict:
            mutation_list = patient_dict[eachPatient]
            if len(mutation_list)>1:
                genes.append(eachGene)
                for eachMutation in mutation_list:
                    if eachPatient in patient_gene_dict:



                        gene_mutation_dict = patient_gene_dict[eachPatient]
                        if eachGene in gene_mutation_dict:
                            gene_mutation_list = gene_mutation_dict[eachGene]
                            gene_mutation_list.append(eachMutation)
                            gene_mutation_dict[eachGene] = gene_mutation_list
                            patient_gene_dict[eachPatient] = gene_mutation_dict
                        else:
                            gene_mutation_dict[eachGene] = [eachMutation]
                            patient_gene_dict[eachPatient] = gene_mutation_dict

                    else:
                        gene_mutation_dict = {}
                        gene_mutation_dict[eachGene] =[eachMutation]

                        patient_gene_dict[eachPatient] = gene_mutation_dict


##
##                gene_mutation_dict = patient_gene_dict[eachGene]
##
##                print 'Gene = ', eachGene
##                genes.append(eachGene)
##                print 'Patient = ', eachPatient
##                print '# of mutations in this patient = ', len(mutation_list)
##                for eachMutation in mutation_list:
##                    print eachMutation['EFF']
##                    print eachMutation['REF'], ',ALT= ', eachMutation['ALT'], ', score = ', eachMutation['score']



    print 'ALL GENES = ', Set(genes)

    pickle.dump(patient_gene_dict, open('patient_gene_dict_comp_heterozygous.pkl', 'w'))

def compound_heterozygosity_2_csv():
    fhout = open('compound_heterozygous_genes.csv', 'w')
    outline_header = 'Patient_Name, Gene_Name, Score, Depth, Functional_Class, EFF' + '\n'
    fhout.write(outline_header)

    data = pickle.load(open('patient_gene_dict_comp_heterozygous.pkl', 'rU'))
    for eachPatient in data:
        patient_gene_dict = data[eachPatient]
        for eachGene in patient_gene_dict:
            mutation_list = patient_gene_dict[eachGene]
            for eachMutation in mutation_list:
                outline = eachPatient + ','+ eachGene + ',' + str(eachMutation['score']) +',' + '#'.join(eachMutation['DP']) + ',' +\
                set_2_string(eachMutation['Effect']) + ','+'#'.join(eachMutation['EFF']) + '\n'
                fhout.write(outline)

    fhout.close()


def set_2_string(aset):

    return '#'.join(list(aset))


def main():
##    compound_heterozygosity()
    compound_heterozygosity_2_csv()

if __name__ == '__main__':
    main()
