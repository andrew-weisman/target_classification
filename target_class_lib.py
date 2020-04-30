# Load the data created by the extract_data() function in target_class_lib.sh
def load_tsv_files(datadir):

    # Import relevant modules
    import pandas as pd
    import json

    # Load the file data
    with open(datadir+'filedata.json') as f:
        filedata = json.load(f)

    # Load the datafiles into Pandas dataframes
    df = pd.concat( [ pd.read_csv(tsv_file, sep='\t').set_index('gene-pretty').drop(columns='gene-ugly') for tsv_file in filedata['tsv_files'] ] , axis='columns', join='outer', ignore_index=False, sort=False).transpose()
    filedata = pd.DataFrame({'filenames': filedata['filenames'], 'weblinks': filedata['weblinks'], 'normalization': df.index, 'idatafiles': filedata['idatafiles'], 'tsv_files': filedata['tsv_files']})
    df = df.reset_index(drop=True)

    # Load the metadata associated with the files
    with open(datadir+'metadata.json') as f:
        metadata = json.load(f)

    # Return all the loaded data in Python-friendly formats
    return(df, filedata, metadata)


# tsv_files = ['tmp1.txt', 'tmp2.txt']
# filenames = ['/data/HS9999.gene.quantification.txt', '/data/HS0825.gene.quantification.txt']
# weblinks = ['https://target-data.nci.nih.gov/HS9999.gene.quantification.txt', 'https://target-data.nci.nih.gov/HS0825.gene.quantification.txt']
# idatafiles = [0, 1]
#df, filedata = load_tsv_files(tsv_files, filenames, weblinks)

# Parameter
datadir = '/XXXX/'

# Load all the data created by the extract_data() function in target_class_lib.sh
df, filedata, metadata = load_tsv_files(datadir)
