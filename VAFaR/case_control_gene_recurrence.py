#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     09/05/2013
# Copyright:   (c) Gaurav Singhal 2013
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------

"""
This script takes the case-gene-by-patient data and control-gene-by-patient-data and compiles a csv file.

"""


import sys
import pickle
from sets import Set

def case_control_gene_recurrence(composite_outfile, case_csv, control_csv, cases_list_name, controls_list_name):

    outFile = open(composite_outfile, 'w')
    geneByPatient_binary_control(outFile, case_csv, control_csv, cases_list_name, controls_list_name)


def geneByPatient_binary_control(outFile, case_csv, control_csv, cases_list_name, controls_list_name):

    fileList = open(case_csv, 'rU').readlines()
    control_patient_list = []

    for each in fileList:

        control_patient_list.append(each.split(',')[0].strip())


    fileList = open(control_csv, 'rU').readlines()

    cases_patient_list = []

    for each in fileList:

        cases_patient_list.append(each.split(',')[0].strip()[:-4])

    cases_gene_patient_dict_pkl = cases_list_name +  '_gene_patient_dict.pkl'
    controls_gene_patient_dict_pkl = controls_list_name + '_gene_patient_dict.pkl'

    case_gene_patient_dict = pickle.load(open(cases_gene_patient_dict_pkl, 'rU'))
    control_gene_patient_dict = pickle.load(open(controls_gene_patient_dict_pkl, 'rU'))

##    outFile = open('Gene_by_Patient_score_v_rare_deleterious_binary.csv', 'w')

# now writing the gene mutation status in a csv format, genes are columns, patients are rows.

    all_genes = list(Set(control_gene_patient_dict.keys()) | Set(case_gene_patient_dict.keys()))

    outline = ','
    for i in range(len(cases_patient_list)):
        outline+='1,'
    for i in range(len(control_patient_list)):
        outline+='0,'

    outline+='\n'
    outFile.write(outline)

    for eachGene in all_genes:
        #for each patient in cases:
            # 1 if present; 0 if absent
        #for each patient in controls:
            # 1 if present; 0 if absent
        #go to next line, next gene
        outline=eachGene+','
        if eachGene in case_gene_patient_dict:
            patient_dict = case_gene_patient_dict[eachGene]


            for eachCase in cases_patient_list:
                if eachCase in patient_dict:
                    #print 1
                    outline+='1,'
                else:
                    outline+='0,'
        else:
            for eachCase in cases_patient_list:
                outline+='0,'

        if eachGene in control_gene_patient_dict:

            patient_dict = control_gene_patient_dict[eachGene]
            for eachControl in control_patient_list:
                if eachControl in patient_dict:
                    outline+='1,'
                else:
                    outline+='0,'
        else:
            for eachControl in control_patient_list:
                outline+='0,'

        outline+='\n'
        outFile.write(outline)


    outFile.close()



def main():
    case_control_gene_recurrence(composite_outfile, case_csv, control_csv, cases_list_name, controls_list_name)

if __name__ == '__main__':
    main()
