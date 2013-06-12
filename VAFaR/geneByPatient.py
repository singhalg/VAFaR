#-------------------------------------------------------------------------------
# Name:        geneRecurrence
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     07/09/2012
# Copyright:   (c) Gaurav Singhal 2012
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------

from sets import Set
import pickle, sys


def geneByPatient(file_list, list_name):
##    fileList = open('control_samples.csv', 'rU').readlines()
    gene_patient_dict = {}

    for each in file_list:

        bamFile = each.split(',')[0].strip()

        pt_name = bamFile

        pickleFile = bamFile+'_snpEff_SnpSift_GWAS_rare_del.pkl'

        mutation_data = pickle.load(open(pickleFile, 'rU'))

        for eachMutation in mutation_data.keys():
            if 'Gene_Name' in mutation_data[eachMutation]:
                genes = mutation_data[eachMutation]['Gene_Name']
                for eachGene in genes:
                    if eachGene in gene_patient_dict:
                        gene_dict_val = gene_patient_dict[eachGene]
                        if pt_name in gene_dict_val:
                            val = gene_dict_val[pt_name]
                            val.append(mutation_data[eachMutation])
                            gene_dict_val[pt_name] = val
                            gene_patient_dict[eachGene] = gene_dict_val
                        else:
##                            val = {}
                            gene_dict_val[pt_name] =[mutation_data[eachMutation]]
                            gene_patient_dict[eachGene] = gene_dict_val


                    else:
                        val = {}
                        val[pt_name] = [mutation_data[eachMutation]]

                        gene_patient_dict[eachGene] = val

    pickle_file_name = list_name + '_gene_patient_dict.pkl'

    fhPickle = open(pickle_file_name, 'w')
    pickle.dump(gene_patient_dict, fhPickle)
    fhPickle.close()

    return gene_patient_dict


def genePatientCountCSV(file_list, list_name):

    gene_patient_dict = geneByPatient(file_list, list_name)

    gene_patient_count = []
    for eachGene in gene_patient_dict.keys():
        patient_dict = gene_patient_dict[eachGene]
        mutation_count = 0
        for eachPT in patient_dict:
            mutation_count+= len(patient_dict[eachPT])

        gene_patient_count.append([eachGene, patient_dict.keys(), len(patient_dict.keys()), mutation_count])

    #saving as pickle
    pickle_name = list_name + '_gene_patient_count.pkl'
    pickle.dump(gene_patient_count, open(pickle_name, 'w'))

    #writing to csv
    csv_name = list_name + '_gene_patient_count.csv'
    fhout = open(csv_name, 'w')
    outline = 'Gene, Patients, Patient_Count, Mutation_count\n'
    fhout.write(outline)

    sorted_gene_patient_count = sorted(gene_patient_count, key=sort_function, reverse=True)

    for eachGene in sorted_gene_patient_count:



        pt_list = ''

        for eachElement in eachGene[1]:
            pt_list+= str(eachElement) + '_'
        outline = eachGene[0] + ',' + pt_list + ',' + str(eachGene[2]) + ',' + str(eachGene[3]) + '\n'


        fhout.write(outline)
    fhout.close()
    return csv_name


def geneByPatientScoreCSV(patient_list, list_name):


    gene_patient_dict_pickle = list_name + '_gene_patient_dict.pkl'



    gene_patient_dict = pickle.load(open(gene_patient_dict_pickle, 'rU'))
    outfile_name = list_name + '_scores.csv'

    outFile = open(outfile_name, 'w')
    outline = ','
    for eachPT in patient_list:
        outline+=eachPT+','
    outline+='\n'
    outFile.write(outline)

    for eachGene in gene_patient_dict.keys():
        outline = eachGene+','
        pt_dict = gene_patient_dict[eachGene]

        for eachPatient in patient_list:
            if eachPatient in pt_dict:
                scores = []
                mutation_data = pt_dict[eachPatient]
                for eachMutation in mutation_data:
                    if 'score' in eachMutation:
                        scores.append(float(eachMutation['score']))
                comp_score=sum(scores)/len(scores)
                outline+=str(comp_score) + ','


            else:
                outline+=','
        outline+='\n'
        outFile.write(outline)

    outFile.close()
    return outfile_name


def geneByPatientEffectCSV(patient_list, list_name):


    gene_patient_dict_pickle = list_name + '_gene_patient_dict.pkl'
    gene_patient_dict = pickle.load(open(gene_patient_dict_pickle, 'rU'))

    outfile_name = list_name + '_effect.csv'


    outFile = open(outfile_name, 'w')
    outline = ','
    for eachPT in patient_list:
        outline+=eachPT+','
    outline+='\n'
    outFile.write(outline)

    for eachGene in gene_patient_dict.keys():
        outline = eachGene+','
        pt_dict = gene_patient_dict[eachGene]

        for eachPatient in patient_list:
            if eachPatient in pt_dict:
                scores = []
                mutation_data = pt_dict[eachPatient]
                for eachMutation in mutation_data:
                    if ('Effect_Impact' in eachMutation) and ('HIGH' in eachMutation['Effect_Impact']) :
##                        scores.append(float(eachMutation['score']))
##                comp_score=sum(scores)/len(scores)
                        effect = ''
                        for eachEffect in eachMutation['Effect']:
                            effect+=eachEffect+'#'
                        outline+=effect + ','


            else:
                outline+=','
        outline+='\n'
        outFile.write(outline)

    outFile.close()
    return outfile_name




def geneByPatientBinaryCSV(patient_list, list_name):

    gene_patient_dict_pickle = list_name + '_gene_patient_dict.pkl'
    gene_patient_dict = pickle.load(open(gene_patient_dict_pickle, 'rU'))

    outfile_name = list_name + '_binary.csv'


    outFile = open(outfile_name, 'w')
    outline = ','
    for eachPT in patient_list:
        outline+=eachPT+','
    outline+='\n'
    outFile.write(outline)

    for eachGene in gene_patient_dict.keys():
        outline = eachGene+','
        pt_dict = gene_patient_dict[eachGene]

        for eachPatient in patient_list:
            if eachPatient in pt_dict:
##                scores = []
##                mutation_data = pt_dict[eachPatient]
##                for eachMutation in mutation_data:
##                    if 'score' in eachMutation:
##                        scores.append(float(eachMutation['score']))
##                comp_score=sum(scores)/len(scores)
                outline+= '1,'


            else:
                outline+='0,'
        outline+='\n'
        outFile.write(outline)

    outFile.close()
    return outfile_name





def sort_function(alist):
    int(alist[3])




def geneRecurrenceCancer():

    files = [   ['PT_10_lung_1200384.vcf','PT_10_brain_1200360.vcf'],
                ['PT_11_lung_1200363.vcf', 'PT_11_brain_1200385.vcf'],
                ['PT_12_lung_1200381.vcf', 'PT_12_brain_1200368.vcf'],
                ['PT_3_lung_1123221.vcf','PT_3_brain_1123220.vcf'],
                ['PT_1_lung_1123233.vcf', 'PT_1_brain_1123232.vcf'],
                ['PT_6_lung_1123227.vcf', 'PT_6_brain_1123226.vcf'],
                ['PT_9_lung_1200382.vcf',  'PT_9_brain_1200383.vcf'],
                ['PT_13_lung_1202946.vcf', 'PT_13_brain_1200370.vcf'],
                [ 'PT_2_lung_1123231.vcf', 'PT_2_brain_1123230.vcf']  ]


    geneDict = {}



    for eachPatient in files:
        primPick = eachPatient[0][:-4] +'_snpEff_del.pkl'
        metPick = eachPatient[1][:-4] +'_snpEff_del.pkl'


        fhPickleInPrim = open(primPick, 'rU')
        primPatientData = pickle.load(fhPickleInPrim)
        fhPickleInPrim.close()

        fhPickleInMet = open(metPick, 'rU')
        metPatientData = pickle.load(fhPickleInMet)
        fhPickleInMet.close()

        end = primPick.find('_lung')
        patientName = primPick[:end]


        primPositions = primPatientData.keys()
        metPositions = metPatientData.keys()

        ''' pick those genes which have
        DP>=20,
        AF>=0.1
        Present in both prim and met for the same allele

        {Gene_name: [{PT_name:recurrence_score},{func_class:mut_info}, {mut_loc: mut_loc_val_prim}, {mut_loc:, mut_loc_val_met}]}
            '''

def mutationRecurrence():

    files = [   ['PT_10_lung_1200384.vcf','PT_10_brain_1200360.vcf'],
                ['PT_11_lung_1200363.vcf', 'PT_11_brain_1200385.vcf'],
                ['PT_12_lung_1200381.vcf', 'PT_12_brain_1200368.vcf'],
                ['PT_3_lung_1123221.vcf','PT_3_brain_1123220.vcf'],
                ['PT_1_lung_1123233.vcf', 'PT_1_brain_1123232.vcf'],
                ['PT_6_lung_1123227.vcf', 'PT_6_brain_1123226.vcf'],
                ['PT_9_lung_1200382.vcf',  'PT_9_brain_1200383.vcf'],
                ['PT_13_lung_1202946.vcf', 'PT_13_brain_1200370.vcf'],
                [ 'PT_2_lung_1123231.vcf', 'PT_2_brain_1123230.vcf']  ]


    # considering mutations only in primary tumors.
    allDeleterious = {}  # here we are considering only deleterious mutations and only in primary tumors.


    for eachTumor in files:
        prim = eachTumor[0]
        delMutPickleName = prim[:-4] + '_snpEff_del.pkl'
        fh = open(delMutPickleName, 'rU')
        delMutationsPrimary = pickle.load(fh)
        fh.close()

        end = delMutPickleName.find('_lung')
        patientName = delMutPickleName[:end]

        delMutLocations = delMutationsPrimary.keys()
        for eachLoc in delMutLocations:
            AF =  delMutationsPrimary[eachLoc]['AF']
            DP = delMutationsPrimary[eachLoc]['DP']
            ALT =  delMutationsPrimary[eachLoc]['ALT']
            for i in range(len(AF)):
                if (float(AF[i]) >= 0.1) and (int(DP[0])>=20):

                    if eachLoc not in allDeleterious:

                            newVal = [[patientName, ALT[i], float(AF[i]), DP[0]]]
                            allDeleterious[eachLoc] = newVal

                    else:
                        val = allDeleterious[eachLoc]

                        val.append([patientName, ALT[i], float(AF[i]), DP[0]])
                        allDeleterious[eachLoc] = val
##                        except:
##                            print delMutationsPrimary[eachLoc]
##                            sys.exit()

        del delMutationsPrimary, delMutLocations

    allDeleteriousPickle = open('All_deleterious.pkl', 'w')
    pickle.dump(allDeleterious, allDeleteriousPickle)
    allDeleteriousPickle.close()


def topRecurrentLocations():

    commonEVS = pickle.load(open('commonEVS.pkl', 'rU'))
    agilentCommon = pickle.load(open('agilent_V3_common.pkl', 'rU'))
    IlluminaCommon = pickle.load(open('Illumina_common.pkl', 'rU'))
    Illumin_agilent = pickle.load(open('Illumina_Agilent_common.pkl', 'rU'))

    allCommon = commonEVS | agilentCommon | IlluminaCommon | Illumin_agilent


    fh = open('All_deleterious.pkl', 'rU')



    allDeleterious = pickle.load(fh)
    fh.close()
    allDelList = []
    locations = allDeleterious.keys()
    for eachLoc in locations:
        if eachLoc not in allCommon:

            allDelList.append([eachLoc, allDeleterious[eachLoc]])

    allDelListSorted = sorted(allDelList, reverse=True, key = myFun)

    fhCommonDelList = open('common_del_list.pkl', 'w')
    pickle.dump(allDelListSorted, fhCommonDelList)

    count = 0
    for eachMostRecurrent in allDelListSorted[:1000]:
        if (int(eachMostRecurrent[1][0][3])>100) :
            print eachMostRecurrent
            print '# of patients having this mutation = ', len(eachMostRecurrent[1])



def myFun(alist):
    return len(alist[1])



def main():
##    VCFtoPickle()
##    mutationRecurrence()
##    topRecurrentLocations()
    geneRecurrence()
    geneByPatient()
    geneByPatient_effect()
    geneByPatient_binary()
if __name__ == '__main__':
    main()
