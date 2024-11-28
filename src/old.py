import logging
import xenaPython as xena

import pandas as pd


class DownloadManager:
    def __init__(self,
                 hub = "https://toil.xenahubs.net",
                 dataset = "tcga_rsem_isoform_tpm"):
        self.hub = hub
        self.dataset = dataset
        self.samples = None
    
    def get_samples(self, hub, dataset, amount=100):
        res = xena.dataset_samples(hub, dataset, amount)
        return res
    
    def get_identifiers(self, hub, dataset):
        """ These are transcripts ENST or gene ENSG """
        res = xena.dataset_field(hub, dataset)
        return(res)
    
    def get_field_codes(self, hub, dataset, feature):
        xena.field_codes(hub, dataset, [feature])

    def get_data(self, hub, dataset, samples):
        all_identifier = self.get_identifiers(hub, dataset)
        res = xena.gene_transcripts(hub, dataset, samples, all_identifier[:10])
        return res

def main():
    # Read in TCGA PANCAN metadata
    meta_path = "data/raw/TCGA_PANCAN_meta.txt"
    metadata = pd.read_csv(meta_path, sep='\t')

    BRCA_meta = metadata[metadata["cancer type abbreviation"] == "BRCA"]
    BRCA_samples = BRCA_meta["sample"].tolist()


    dl_manager = DownloadManager()
    dl_manager.samples = BRCA_samples
    #res_sample = dl_manager.get_samples(dl_manager.hub, dl_manager.dataset)

    #test = set(BRCA_samples) & set(res_sample)
    
    #res_transcripts = dl_manager.get_identifiers(dl_manager.hub, dl_manager.dataset)

    res_data = dl_manager.get_data(dl_manager.hub, dl_manager.dataset, BRCA_samples[:10])
    print(res_data)
    """
    def getFeatureCodes(hub, dataset, feature):
        # identify feature is categorical or continuous, if it is categorical, there will be codes associated with it, otherwise, codes will be none
        codes = xena.field_codes(hub, dataset, [feature])[0]['code']
        if codes:
            codes = codes.split('\t')
        return codes

    def getClinicalData(hub, cohort, target_feature):
        # find out all the datafiles (i.e. datasets) belong to a cohort
        datasets = xena.dataset_list(hub, [cohort])
        # filter to just clinicalMatrix type of data files
        datasets = list(filter(lambda x: x['type'] == 'clinicalMatrix', datasets))
        # collect all the phynotype features and their associated dataset from all the clinicalMatrix datasets
        features = []
        for dataset in datasets:
            for feature in xena.dataset_field(hub, dataset['name']):
                features.append([feature, dataset['name']])
        # find the target_feature among all the features, and the dataset it comes from 
        xenafield = list(filter(lambda f: f[0] == target_feature, features))
        if len(xenafield) == 0:
            print (target_feature, "not found")
            return [None, None]
        elif len(xenafield) == 1:
            dataset = xenafield[0][1]
            # query to get all the data from the target_feature
            # first, get all the samples in the cohort (a bit slower) 
            # samples = xena.cohort_samples(hub, cohort, None)
            # all the samples in the dataset (a bit faster), either will work, 		
            samples = xena.dataset_samples(hub, dataset, None)
            # second, get the data
            [position, [data]] = xena.dataset_probe_values (hub, dataset, samples, [target_feature])
            # thrid, identify feature is categorical or continuous, if it is categorical, there will be codes associated with it,
            codes = getFeatureCodes(hub, dataset, target_feature)
            if codes:
                data = [codes[int(x)] if x != 'NaN' else 'NaN' for x in data]
            return [samples, data]
        else:
            print ("there are more than one features match", target_feature)
            return [None, None]

    hub = 'https://tcga.xenahubs.net'
    cohort = 'TCGA Ovarian Cancer (OV)'
    target_feature = 'age_at_initial_pathologic_diagnosis'
    target_feature = 'sample_type'
    [samples, data] = getClinicalData(hub, cohort, target_feature)
    samples[:10], data[:10]
    """

    #xena.dataset_gene_probe_avg(hub, dataset, samples, genes)

if __name__ == "__main__":
    main()
    
    
from pathlib import Path
import logging
import sys

# Print all command-line arguments
print("Arguments passed to the script:", sys.argv)



class deepTMHMMRunner:
    def __init__(self,
                 fasta_path
        ):
        self.fasta_path = fasta_path
    
    def run_deeptmhmm(self):
        command = f"biolib run DTU/DeepTMHMM --fasta {self.fasta_path}"
        
    def remove_to_data_folder(self)
        command =

def main():
    # Access specific arguments
if len(sys.argv) > 1:
    print("First argument:", sys.argv[1])
    
    deeptmhmm_runner = deepTMHMMRunner(fasta_path = "data/processed/top100.fa")
    deeptmhmm_runner.run_deeptmhmm()

if __name__ == "__main__":
    main()
                
