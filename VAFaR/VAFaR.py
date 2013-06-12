#-------------------------------------------------------------------------------
# Name:        VAFaR
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     01/06/2013
# Copyright:   (c) Gaurav Singhal 2013
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------
import sys
import runAnnotations
import read_dbSNP_VAFaR
import variantFilteringPickling
import extractControlVCF
import geneByPatient
import case_control_gene_recurrence

"""
This script brings together the functionality of all the core scripts. The
sequential calling of methods from the core scripts constitues an analysis pipeline.
These methods can be called individually, or run directly via the core scripts.

@cases: list of vcf files from cases (eg, tumor tissue, genetic disease patiets)
@controls: list of vcf files, controls (healthy human/sample or tissue)
@dbSNP_common_all_vcf_file: common_all.vcf (download the latest file from dbSNP ftp server)
@dbSNP_common_clinical_vcf_filename: common_clinical vcf file from dbSNP ftp server
@cases_list_name:The gene by patient dict pickle for cases will be saved with this
name, cases_list_name_gene_patient_dict.pkl.
@cases_patient_list: List of patient names. These names need to be unique sample names.
@composite_csv: File name with which the final composite gene by patient data (from
both cases and controls) need to be saved in.

"""


def VAFaR(cases, controls, dbSNP_common_all_vcf_file, dbSNP_common_clinical_vcf_filename, cases_list_name, controls_list_name, cases_patient_list, controls_patient_list, composite_csv):


    # get control vcfs (optional, use if you do not have control vcfs of your own
    # and want to use vcf files extracted from NIEHS 95 sample variant data
    # Note that phenotype information is NOT available for NIEHS's 95 samples
    controls = extractControlVCF.vcf_extractor()

    # run vcf annotations
    # refGenome refers to the database version of snpEff
    # for details of refGenome see http://snpeff.sourceforge.net/download.html#databases
    #SnpSift_ref_db refers to the location of the dbNSFP database. For details,
    # see http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0b3.readme.txt

    runAnnotations.snpEffSnpSift(cases, controls, refGenome='GRCh37.69',
                                SnpSift_ref_db="./dbNSFP2.0b3.txt",
                                verbose=True)


    # run dbSNP's common_all.vcf and specify build and maf, and construct a
    # dict object
    # if the specified build is greater than the build of the common_all.vcf
    # then all variants would be included in dict object
    dbSNP_common = read_dbSNP_VAFaR.read_dbSNP(dbSNP_common_all_vcf_file, build=135, maf=0.03, verbose=True)

    #download the latest vcf file having common snps with clinical significance.
    # common clinical vcf file names look like "common_and_clinical_20130226.vcf"
    dbSNP_common_clinical = read_dbSNP_VAFaR.read_clin_sign_snps(dbSNP_common_clinical_vcf_filename)

    #recurrence_level indicates that a variant which occurs twice or more within
    # the set of cases or controls will be considered as a false positive and
    # not considered further.
    variantFilteringPickling.VCFtoPickle(cases, controls, dbSNP_common, dbSNP_common_clinical, recurrence_level=2)



    # geneByPatient representations
    #cases is the list of vcf files which belong to cases
    #cases_list_name is the name with which the geneByPatient dict will be saved
    geneByPatient.geneByPatient(cases, cases_list_name)
    geneByPatient.geneByPatient(controls, controls_list_name)

    #here, cases_patient_list is the list of patient names.
    cases_gene_patient_csv = geneByPatient.geneByPatientBinaryCSV(cases_patient_list, cases_list_name)
    controls_gene_patient_csv = geneByPatient.geneByPatientBinaryCSV(controls_patient_list, controls_list_name)

    #composite_csv is the name of the csv that you specify, and the cases and
    # control data will be saved in this csv file
    case_control_gene_recurrence.case_control_gene_recurrence(composite_csv, cases_gene_patient_csv, controls_gene_patient_csv, cases_list_name, controls_list_name)






def main():
    # use sys args or call the function from another python script.
    VAFaR(cases, controls, dbSNP_common_all_vcf_file, dbSNP_common_clinical_vcf_filename, cases_list_name, controls_list_name, cases_patient_list, controls_patient_list, composite_csv)

if __name__ == '__main__':
    main()
