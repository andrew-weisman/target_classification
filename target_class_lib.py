# Load the data created by the extract_data() function in target_class_lib.sh
def load_tsv_files(datadir):

    # Import relevant modules
    import pandas as pd
    import json

    # Load the file data
    with open(datadir+'filedata.json') as f:
        filedata = json.load(f)

    # Load the datafiles into Pandas dataframes
    df = pd.concat( [ pd.read_csv(tsv_file, sep='\t').set_index('gene-symbol').drop(columns='gene-id') for tsv_file in filedata['tsv_files'] ] , axis='columns', join='outer', ignore_index=False, sort=False).transpose()
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

# # Parameter
# datadir = '/XXXX/'

# # Load all the data created by the extract_data() function in target_class_lib.sh
# df, filedata, metadata = load_tsv_files(datadir)


# This function is modified from https://rest.ensembl.org/documentation/info/symbol_post; it basically performs a POST lookup of a list of symbols
def ensembl_symbol_lookup(symbols_list):

    # The maximum POST size (and therefore list size) is 1000
    if len(symbols_list) > 1000:
        print('ERROR: POST lookup list is too long ({})'.format(len(symbols_list)))
        exit

    # Import relevant modules
    import requests, sys
 
    # Set some constants
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    # Query the REST server
    #r = requests.post(server+ext, headers=headers, data='{ "symbols" : ["BRCA2", "BRAF" ] }')
    r = requests.post(server+ext, headers=headers, data=str({'symbols': symbols_list}).replace('\'','"'))
    
    # Ensure everything is okay
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    # Return the result
    #return(repr(r.json()))
    return(r.json())


# Return a list of symbols and corresponding IDs identified using the symbol lookup REST function of Ensembl given the symbols in a text file called to_lookup.txt
def get_lookup_list(project_dir, max_post_size=20, max_iter=-1):

    # Import relevant module
    import os

    # Initialize some overall variables
    iiter = 0
    lookup_list = []

    # Open the datafile containing names to look up using the basic symbol lookup in REST
    with open(os.path.join(project_dir, 'data', 'to_lookup.txt'), 'r') as f:

        # For each symbol in the file...
        iname = 0
        symbols_list = []
        result_processed = False
        for name in f:

            # Add the symbol to a list of symbols
            iname = iname + 1
            symbols_list.append(name.rstrip())

            # If the symbol list has become the maximum POST size...            
            if iname == max_post_size:

                # Run the query
                res = ensembl_symbol_lookup(symbols_list)

                # Process the result of the query, adding any found symbol-ID pairs to the main lookup list
                keys = res.keys()
                for key in keys:
                    lookup_list.append([key, res.get(key).get('id')])

                # Record that we've run a query
                result_processed = True
                iiter = iiter + 1

                # Print what we did
                print('[{}] Given {} symbols, we identified {} of them'.format(iiter, len(symbols_list), len(keys)))

                # If we've hit the maximum number of queries, get out of the loop
                if iiter == max_iter:
                    break

                # Reset the things we're incrementing and resume with the "for" loop
                iname = 0
                symbols_list = []
                result_processed = False

    # If the we left the loop because we ran out of names in the file but hadn't yet run a query, then we still have data to process
    if (not result_processed) and (iname != 0):
        res = ensembl_symbol_lookup(symbols_list)
        keys = res.keys()
        for key in keys:
            lookup_list.append([key, res.get(key).get('id')])
        result_processed = True
        iiter = iiter + 1
        print('[{}] Given {} symbols, we identified {} of them'.format(iiter, len(symbols_list), len(keys)))

    # Return the lookup list
    return(lookup_list)
