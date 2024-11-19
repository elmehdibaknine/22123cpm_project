# 22123cpm_project

## Get TCGA Data (rough guidelines)
Download breast cancer gene level identifiers:
https://xenabrowser.net/datapages/?host=https%3A%2F%2Fgdc.xenahubs.net&dataset=TCGA-BRCA.star_counts.tsv&allIdentifiers=true&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

Download all transcript cancer data:
https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_expected_count.gz

Run the follwing bash commands:
## Get headers for breast cancer only
awk 'NR==1 {for (i=1; i<=NF; i++) print $i}' ../TCGA-BRCA.star_counts.tsv > headers_small

## Remove the final letter of all rows (except first). For some reason, there is either an 'A' or 'B' at the end of each line
sed '1!s/.$//' headers_small > temp
mv temp headers_small

## Get headers of all transcripts
awk 'NR==1 {for (i=1; i<=NF; i++) print $i}' tcga_expected_count > headers_large

## Get the intersection of samples from ALL and BREAST (remember to add the "gene_id" col names)
awk 'NR==FNR {arr[$1]; next} $1 in arr' headers_large headers_small > intersection

## Subset the ALL data to only contain columns from BREAST (run from a .sh file):
```
# First, create a list of column names from headers.txt
columns=$(<intersection)

# Now use awk to extract the corresponding columns from GTEx.txt based on matching headers
awk -v cols="$columns" 'BEGIN {
    # Split the columns into an array
    split(cols, colArray, "\n")
    for (i in colArray) {
        columnName[colArray[i]] = 1
    }
}
NR == 1 {
    # For the header row, identify the columns that match the ones in headers.txt
    for (i = 1; i <= NF; i++) {
        if ($i in columnName) {
            selectedColumns[i] = $i
        }
    }
    # Print matching headers
    for (i = 1; i <= NF; i++) {
        if (i in selectedColumns) {
            printf "%s\t", $i
        }
    }
    print ""
}
NR > 1 {
    # For data rows, print the corresponding columns
    for (i = 1; i <= NF; i++) {
        if (i in selectedColumns) {
            printf "%s\t", $i
        }
    }
    print ""
}' tcga_expected_count > filtered_TCGA
```

## Get GTEx Data (rough guideline)
Download breast cancer gene level data:
https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/counts-by-tissue/gene_reads_v10_breast_mammary_tissue.gct.gz

Download all transcript cancer data:
https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_expected_count.txt.gz

Run the follwing bash commands:
## Get headers for breast cancer only
awk 'NR==1 {for (i=1; i<=NF; i++) print $i}' gene_reads_breast_mammary_tissue.gct > headers_small

## Get headers of all transcripts
awk 'NR==1 {for (i=1; i<=NF; i++) print $i}' GTEx_Analysis_2022-06-06_v10_RSEMv1.3.3_transcripts_expected_count.txt > headers_large

## Get the intersection of samples from ALL and BREAST (remember to add the "gene_id" col names)
awk 'NR==FNR {arr[$1]; next} $1 in arr' headers_large headers_small > intersection

## Subset the ALL data to only contain columns from BREAST (run from a .sh file):
```
# First, create a list of column names from headers.txt
columns=$(<intersection)

# Now use awk to extract the corresponding columns from GTEx.txt based on matching headers
awk -v cols="$columns" 'BEGIN {
    # Split the columns into an array
    split(cols, colArray, "\n")
    for (i in colArray) {
        columnName[colArray[i]] = 1
    }
}
NR == 1 {
    # For the header row, identify the columns that match the ones in headers.txt
    for (i = 1; i <= NF; i++) {
        if ($i in columnName) {
            selectedColumns[i] = $i
        }
    }
    # Print matching headers
    for (i = 1; i <= NF; i++) {
        if (i in selectedColumns) {
            printf "%s\t", $i
        }
    }
    print ""
}
NR > 1 {
    # For data rows, print the corresponding columns
    for (i = 1; i <= NF; i++) {
        if (i in selectedColumns) {
            printf "%s\t", $i
        }
    }
    print ""
}' GTEx_Analysis_2022-06-06_v10_RSEMv1.3.3_transcripts_expected_count.txt > filtered_GTEx
```