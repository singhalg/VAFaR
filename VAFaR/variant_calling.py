#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     15/05/2013
# Copyright:   (c) Gaurav Singhal 2013
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------


import sys
import pickle

def raw_bam_2_VCF():

    sampleList = open('SJ_renal_dysplasia_bam_files.csv', 'rU').readlines()
    fastq_file_list = open('fastq_files.csv', 'rU').readlines()


    fhoutName = 'SJ_Renal_variant_calling_unique.sh'





    fhout = open(fhoutName, 'w')

    wd = " /BlueArc-scratch/gsinghal/MWatson/SJ_renal/Dysplasia_1000Genes_Data/"
    options = "#!/bin/bash" + '\n'+ "#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=40gb \n" + "#PBS -q dque \n" + "#PBS -N   "+fhoutName[:-3] + '\n' + "#PBS -d   " + wd + '\n' + '#PBS -m abe '+'\n'
    fhout.write(options)
##    fastqFiles = {}

    for each in sampleList:

        sample_name = each.split(',')[0].strip()[:-4]


##    fastqFiles = open('listOfFiles', 'rU').readlines()
##
##    fhoutName = 'Samtools_sort_mpileup_vcf.sh'
##    fhout = open(fhoutName, 'w')
##    options1 = "#!/bin/bash" + '\n'+ "#PBS -l nodes=1:ppn=8,walltime=12:00:00,vmem=40gb \n" + "#PBS -N   "+fhoutName[:-3] + '\n' + "#PBS -d  /BlueArc-scratch/gsinghal/MWatson/Validation/hg19_Ref/ " + '\n' + '#PBS -m abe '+'\n'
##    fhout.write(options1)
##
##    for eachPair in fastqFiles:
##        files = eachPair.strip().split()
##        se = files[0].find('_L001')
##        sampleName = files[0][:se]


        samtools_sort = " /home/gsinghal/bin/samtools sort "+ sample_name +'_unique_reads.bam  '+ sample_name+'_sorted' + '\n'

        fhout.write(samtools_sort)
        fhout.write('sleep 2')
        fhout.write('\n')


        # picard remove duplicates
        inputfileName = sample_name+'_sorted.bam'
        outputfileName = inputfileName[:-4]+ '_rmdup.bam'

        picard_remove_duplicates = "java -jar /export/picard-tools-1.67/MarkDuplicates.jar INPUT=" + inputfileName + " OUTPUT=" +outputfileName + " METRICS_FILE=metrics_rmdup.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true "  '\n'

        fhout.write(picard_remove_duplicates)
        fhout.write('sleep 2')
        fhout.write('\n')


        samtools_index =  " /home/gsinghal/bin/samtools index "+ sample_name +'_sorted_rmdup.bam  '+ '\n'

        fhout.write(samtools_index)
        fhout.write('sleep 2')
        fhout.write('\n')




        mpileup = " /home/gsinghal/bin/samtools mpileup  -f  /BlueArc-scratch/gsinghal/MWatson/hg19/hg19.fa  "+ sample_name +'_sorted_rmdup.bam  | bcftools view -bvcg - > ' + sample_name+'_var.raw.bcf '+ '\n'

        fhout.write(mpileup)

        fhout.write('sleep 2')
        fhout.write('\n')

        vcf_calling = "home/gsinghal/bin/bcftools view " + sample_name + "_var.raw.bcf | vcfutils.pl varFilter -D10000 > "  + sample_name + "_unique_samtools.vcf"

##        VarScanCmd = 'java -jar  ./VarScan.v2.3.2.jar mpileup2snp ' + sampleName +  '_hg19_ref_mpileup.txt ' + '  --min-reads2  5  --min-var-freq  0.05  --p-value  0.05  --output-vcf  >  ' + sampleName + '_hg19_ref.vcf' + '\n'

        fhout.write(vcf_calling)
        fhout.write('\n')
        fhout.write('sleep 2')
        fhout.write('\n')

        inputfileName = sample_name + "_unique_samtools.vcf"
        outFile = sample_name + '_unique_samtools_snp_in_target.vcf'

        bedIntersect = 'home/gsinghal/myBin/BEDTools-Version-2.16.2/bin/bedtools intersect -a ' + inputfileName + ' -b  SJ_renal_1000Genes_miRNA_sorted_merged_nochr.bed  -wa >  ' + outFile
        fhout.write(bedIntersect)
        fhout.write('\n')
        fhout.write('sleep 2')
        fhout.write('\n')




def samtools_flagStat():
    fastqFiles = open('listOfFiles', 'rU').readlines()


    fhoutName = 'NovoAlign_amplicons_samtools_flagstat.sh'
    fhout = open(fhoutName, 'w')
    options1 = "#!/bin/bash" + '\n'+ "#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=40gb \n" + "#PBS -N   "+fhoutName[:-3] + '\n' + "#PBS -d  /BlueArc-scratch/gsinghal/MWatson/Validation/ " + '\n' + '#PBS -m abe '+'\n'
    fhout.write(options1)

    for eachPair in fastqFiles:
        files = eachPair.strip().split()
        se = files[0].find('_L001')
        sampleName = files[0][:se]

        samtools_flagstat_cmd = " /home/gsinghal/bin/samtools flagstat "+ sampleName +'_amplicon_ref.bam   >   '+ sampleName+'_amplicon_ref.stats' + '\n'

        print samtools_flagstat_cmd

        fhout.write(samtools_flagstat_cmd)

        fhout.write('sleep 2')
        fhout.write('\n')
    fhout.close()


def main():
    raw_bam_2_VCF()

if __name__ == '__main__':
    main()
