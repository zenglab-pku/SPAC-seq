######################################################
# Author: Augus Zhang
# Date: 2023-04-11
# Version: v6.0.0
# This script is used to run the SAW pipeline on a collection of paired-end reads.
######################################################

#!/bin/bash

if [[ $# -lt 6 ]];then
    echo "usage: sh $0 -speciesName -tissueType -splitCount -summaryFile -outDir -threads -doCellBin
    -speciesName : specie of the sample, currently support [homo] sapiens and [mus] musculus; also we have [homo_BGI] for BGI version of reference genome
    -tissueType : tissue type of the sample
    -splitCount : count of splited stereochip mask file, usually 16 for SE+Q4 fq data and 1 for PE+Q40 fq data
    -summaryFile :  a txt summarized the fastq file path, maskFile path ,ImageRecord[optional] and ImageCompressed[optional] file path, shall be ABSOLUTE path
    -outDir : output directory path, shall be ABSOLUTE path
    -doCellBin : [Y/N]
    -threads : the number of threads to be used in running the pipeline
    "
    exit
fi

# Get parameters
while [[ -n "$1" ]]
do
    case "$1" in
        -speciesName) species="$2"
            shift ;;
        -tissueType) tisType="$2"
            shift ;;
        -splitCount) splitCnt="$2"
            shift ;;
        -summaryFile) sumFile="$2"
            shift ;;
        -outDir) outDir="$2"
            shift ;;
        -threads) threads="$2"
            shift ;;
        -doCellBin) doCell="$2"
            shift ;;
    esac
        shift
done


# Set up software
path=$outDir
sifPath='/home/apps/SAW-6.0.0/SAW_06.0.0.sif'

# Reference genome: human: GRCH38, mouse: mm10; STAR: 100bp all from genecode
# For above reference genome detailed infor, please refer to the referance.jason file
# For homo_BGI, we use reference genome from https://nervous-puffin-fc6.notion.site/8d0a03113aa044afb90e4709a0e5eb01?v=bfb2de13631042a995871be9c899a718

if [[ $species = 'homo' ]];then
    refIndexDir=/home/apps/SAW-6.0.0/geneEx-GRCh38/star/
    annFile=/home/apps/SAW-6.0.0/geneEx-GRCh38/genes/genes.gtf
    Gsize=3211173485
elif [[ $species = 'homo_BGI' ]];then
    refIndexDir=/home/apps/SAW-6.0.0/geneEx-GRCh38-BGI/star/
    annFile=/home/apps/SAW-6.0.0/geneEx-GRCh38-BGI/genes/genes.gtf
    Gsize=3418782678
elif [[ $species = 'mus' ]];then
    refIndexDir=/home/apps/SAW-6.0.0/geneEx-mm10/star/
    annFile=/home/apps/SAW-6.0.0/geneEx-mm10/genes/genes.gtf
    Gsize=2796735707
else 
    echo "Currently SAW donot support your species, please check the spellings!"
    exit
fi

# Unlock max file limit
ulimit -c 10000000000

echo "=======>Your species is $species, reading the summary file..."

# Read the sample information, ingore the first line
cat $sumFile | tail -n +2 | while read line
do
    startTime=`date +%Y%m%d-%H:%M`
    sampleId=`echo $line | awk '{print $1}'`
    fqPath=`echo $line | awk '{print $2}'`
    fqFile1=`echo $line | awk '{print $3}'`
    fqFile2=`echo $line | awk '{print $4}'`
    maskFile=`echo $line | awk '{print $5}'`
    ImReFile=`echo $line | awk '{print $6}'`
    ImCoFile=`echo $line | awk '{print $7}'`
    
    echo "Processing sample $sampleId"
    echo "Running stereoPipleine for $sampleId"
    
    # Reorganize fq files with full path
    fq1Array=(`echo $fqFile1 | tr ',' ' '` ) 
    fq1New=""
    for var in ${fq1Array[@]}
    do
        fq=$fqPath$var
        fq1New=$fq1New$fq","
    done 
    fq1New=${fq1New%?}

    fq2Array=(`echo $fqFile2 | tr ',' ' '` ) 
    fq2New=""
    for var in ${fq2Array[@]}
    do
        fq=$fqPath$var
        fq2New=$fq2New$fq","
    done 
    fq2New=${fq2New%?}
    
    if [[ -n "$ImReFile" && -n "$ImCoFile" ]]; then
        echo "Image file provided[optional]! Will do rigister."
        echo "=======>All set up, running stereoPipeline..."
        sh 'stereoPipeline_v6.0.0.sh'\
                            -genomeSize $Gsize -splitCount $splitCnt \
                            -maskFile "${maskFile}"\
                            -fq1 "${fq1New}"\
                            -fq2 "${fq2New}"\
                            -imageRecordFile $ImReFile\
                            -imageCompressedFile $ImCoFile\
                            -refIndex $refIndexDir\
                            -speciesName $species -tissueType $tisType \
                            -annotationFile $annFile\
                            -outDir "${path}${sampleId}" \
                            -doCellBin $doCell -threads $threads -sif $sifPath
    else
        echo "Image file not provided[optional]!"
        echo "=======>All set up, running stereoPipeline..."
        sh 'stereoPipeline_v6.0.0.sh'\
                            -genomeSize $Gsize -splitCount $splitCnt \
                            -maskFile "${maskFile}"\
                            -fq1 "${fq1New}"\
                            -fq2 "${fq2New}"\
                            -refIndex $refIndexDir\
                            -speciesName $species -tissueType $tisType \
                            -annotationFile $annFile\
                            -outDir "${path}${sampleId}" \
                            -doCellBin $doCell -threads $threads -sif $sifPath
    fi
    endTime=`date +%Y%m%d-%H:%M`
    echo "$startTime ---> $endTime ""Finish processing $sampleId"
done

echo `date` " All done! "