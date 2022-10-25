# a module to integrate single-cell ATAC-seq data with genomic annotations to classify Peaks
import sys
import pandas as pd
import subprocess as sp
import csv

# The grabPeaks function uses the bigBed file containing Peaks from ATAC-seq
# to build a dataframe of Peaks for further analysis
# and provide some information on peak populations
# returns the Peaks dataframe
def grabPeaks(peaksfilepath):

    # create a df for storing Peaks
    peaks_df = pd.read_csv(peaksfilepath, delimiter='\t', usecols=[0, 1, 2])

    # create a new column with peak length
    peaks_df['peakLength'] = peaks_df.chromEnd - peaks_df.chromStart

    # create a new column with peak center
    peaks_df['peakCenter'] = ((peaks_df.chromEnd - peaks_df.chromStart)/2).round()

    return peaks_df

# The grabAnnotations function leverages the bigBedToBed tool from UCSC
# to access specific tracks from the ENCODE project at specific coordinates (#chrom  chromStart  chromEnd)
# returns the annotations linked to the input region
def grabAnnotations(trackName, chrom, chromStart, chromEnd):

    # initializing data frame to extract genomic annotations
    annotations_df = pd.DataFrame([])

    # leveraging the bigBedToBed tool to access the downloaded ENCODE tracks
    # save stout into .bed file
    annotations_bed_file = sp.getoutput('bigBedToBed ' + trackName + ' -chrom=' + str(chrom) + ' -start=' + str(chromStart) + ' -end=' + str(chromEnd) + ' stdout')

    # loading .bed file content into data frame
    annotations_csv_reader = csv.reader(annotations_bed_file.splitlines(), delimiter="\t")

    # populating the annotations data frame with rows from the standard output
    for row in annotations_csv_reader:

        series = pd.Series(row)
        annotations_df = annotations_df.append(series, ignore_index=True)
        print(annotations_df)

    # using encodeCcreCombined as genomic annotation track
    # either for humans(from http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb)
    # or mice (from http://hgdownload.soe.ucsc.edu/gbdb/mm10/encode3/ccre/encodeCcreCombined.bb)

    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'ccre', 'encodeLabel', 'zScore', 'ucscLabel', 'accessionLabel', 'description']
    annotations_df = annotations_df.rename(columns=pd.Series(header))

    return annotations_df

# The classifyPeaks function iterates over Peaks in the data frame
# collects corresponding tracks from the annotations dataframe
# and returns the Peaks labeled with annotation-based classification

# peaksDataFrame is the data frame containing Peaks to classify
# trackName is the name of the ENCODE track to consider for classification
# this script only handles cCRECombined tracks and in particular uscsLabels
# columnsNames is the list of columns from the specified tracks
# to be considered for peak classification.
# it saves the labeled Peaks to a .csv file
def classifyPeaks(peaksDataFrame, trackName, columnsNames):

    peakLabel = [''] * peaksDataFrame.shape[0]

    for peakIndex, row in peaksDataFrame.iterrows():

        annotations = grabAnnotations(trackName, row['#chrom'], row['chromStart'], row['chromEnd'])

        if not annotations.empty:

            peakLabel[peakIndex] = '\t'.join(annotations[columnsNames[0]].tolist())

    peaksDataFrame[str(trackName) + '_' + str(columnsNames[0])] = pd.Series(peakLabel[:])

    return peaksDataFrame

    
if __name__ == '__main__':

    # sys.argv[1] -> Peaks filename
    # sys.argv[2] -> annotations track filename

    peaks_path = sys.argv[1]
    annotation_track_path = sys.argv[2]
    dataset_name = os.path.basename(sys.argv[1])

    peaks = grabPeaks(peaks_path)

    labeledPeaks = classifyPeaks(peaks, annotation_track_path, columnsNames=['ucscLabel'])
    # saving the classifiedPeaks_df to the output file
    labeledPeaks.to_csv(path_or_buf= '../TMPResults/labeled_peaks/' + dataset_name + '/cCRE_labeled_peaks.tsv', sep='\t')

    print('The Peaks file at ', peaks_path, ' was labeled with genomic annotations from ', annotation_track_path)
