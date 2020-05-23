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


# This function performs a POST or GET from a list containing "names", which can be symbols, IDs, or other
# This is modified from https://rest.ensembl.org/documentation/info/symbol_post (and the other two endpoints' links)
def ensembl_name_request(endpoint, names_list, wait_time=1):

    # Import relevant modules
    import requests, sys, time

    # Wait wait_time seconds to not reach the API limit
    time.sleep(wait_time)
 
    # Set the constant
    server = "https://rest.ensembl.org"

    # Set four settings depending on the endpoint
    if endpoint == 'symbol':
        max_list_len = 1000 # maximum POST size
        ext = "/lookup/symbol/homo_sapiens"
        headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
        data = str({'symbols': names_list}).replace('\'','"')
    elif endpoint == 'id':
        max_list_len = 1000 # maximum POST size
        ext = "/lookup/id"
        headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
        data = str({'ids': names_list}).replace('\'','"')
    elif endpoint == 'xref':
        max_list_len = 1
        #ext = "/xrefs/symbol/homo_sapiens/BRCA2?"
        ext = "/xrefs/symbol/homo_sapiens/" + names_list[0]
        headers = { "Content-Type" : "application/json"}
        data = None

    # Check the names list length
    if len(names_list) > max_list_len:
        print('ERROR: Names list is too long ({}) for the {} endpoint (max length is {})'.format(len(names_list), endpoint, max_list_len))
        exit

    # Query the REST server
    if data is not None:
        r = requests.post(server+ext, headers=headers, data=data)
    else:
        r = requests.get(server+ext, headers=headers)
    
    # Ensure everything is okay
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    # Return the result
    return(r.json())


# Run and process the query, adding any found name-ID pairs to the main lookup list (or otherwise, add None)
def run_and_process_query(endpoint, names_list, lookup_list, iiter, nidentified_tot):

    # Run the query
    res = ensembl_name_request(endpoint, names_list)

    # If the endpoint is 'xref', convert the result to something that was already processable for the 'symbol' and 'id' endpoints
    if endpoint == 'xref':
        if res:
            res = {names_list[0]: res[0]}
        else:
            res = {}

    # Process the result of the query, adding any found name-ID pairs to the main lookup list
    keys = res.keys()
    bad_keys = 0
    for key in keys:
        key_dict = res.get(key)
        if key_dict is not None:
            lookup_list.append([key, key_dict.get('id')])
        else:
            lookup_list.append([key, None])
            bad_keys = bad_keys + 1

    # Additionally, add None to the names that were not found
    for key in (set(names_list)-set(keys)):
        lookup_list.append([key, None])

    # Record that we've run a query
    result_processed = True
    iiter = iiter + 1

    # Print what we did
    nidentified = len(keys) - bad_keys
    nidentified_tot = nidentified_tot + nidentified
    print('[{}] Given {} names, we identified {} of them'.format(iiter, len(names_list), nidentified))

    return(lookup_list, result_processed, iiter, nidentified_tot)


# Return a list of names and corresponding IDs (or None) identified using a REST function of Ensembl given the names in a supplied text file
def get_lookup_list(endpoint, names_file, max_list_len=None, max_iter=-1):

    # Import relevant module
    import os

    # Constant
    max_list_len_dict = {'symbol': 1000, 'id': 1000, 'xref': 1}

    # Set the maximum list length (basically maximum size of either POST or GET, which depends on the endpoint)
    if max_list_len is None:
        max_list_len = max_list_len_dict[endpoint]

    # Initialize some overall variables
    iiter = 0
    lookup_list = []
    nidentified = 0

    # Open the datafile containing names to query
    with open(names_file, 'r') as f:

        # Initialize the things we're incrementing
        iname = 0
        names_list = []
        result_processed = False

        # For each name in the file...
        for name in f:

            # Add the name to a list of names
            iname = iname + 1
            names_list.append(name.rstrip()) # remove the newline at the end of the string

            # If the name list has become the maximum list length...            
            if iname == max_list_len:

                # Run and process the query, adding any found name-ID pairs to the main lookup list (or otherwise, add None)
                lookup_list, result_processed, iiter, nidentified = run_and_process_query(endpoint, names_list, lookup_list, iiter, nidentified)

                # If we've hit the maximum number of queries, get out of the loop
                if iiter == max_iter:
                    break

                # Reset the things we're incrementing and resume with the "for" loop
                iname = 0
                names_list = []
                result_processed = False

    # If the we left the loop because we ran out of names in the file but hadn't yet run a query, then we still have data to process
    if (not result_processed) and (iname != 0):
        lookup_list, result_processed, iiter, nidentified = run_and_process_query(endpoint, names_list, lookup_list, iiter, nidentified)

    # Return the lookup list
    return(lookup_list, nidentified)


# Generate or read in the initial lookup lists that can be calculated by querying the Ensembl REST API
def get_initial_lookup_lists(project_dir):

    # Import relevant modules
    import os, sys
    gmb_dir = '/data/BIDS-HPC/private/projects/gmb/checkout'
    if gmb_dir not in sys.path:
        sys.path.append(gmb_dir)
    import time_cell_interaction_lib as tci # we need this to get the pickle functions

    # If the file containing the initial lookup lists does not already exist...
    if not os.path.exists(os.path.join(project_dir,'data','initial_lookup_lists.pkl')):

        # Query the Ensembl REST API to see whether we can identify the gene "names"
        id_lookup_list, id_nidentified = get_lookup_list('id', os.path.join(project_dir,'data','unique_ensembl_ids.txt'))
        symbol_lookup_list, symbol_nidentified = get_lookup_list('symbol', os.path.join(project_dir,'data','unique_other_names.txt'))

        # Run a check on the number of Nones in the lists; otherwise, save the lookup lists to disk
        if (sum([ 0 if x[1] is None else 1 for x in symbol_lookup_list ]) != symbol_nidentified) or (sum([ 0 if x[1] is None else 1 for x in id_lookup_list ]) != id_nidentified):
            print('ERROR: Inconsistent number of Nones in the lookup lists')
            exit
        else:
            tci.make_pickle([symbol_lookup_list, id_lookup_list], os.path.join(project_dir,'data'), 'initial_lookup_lists.pkl')

    # Otherwise, read it in
    else:
        [symbol_lookup_list, id_lookup_list] = tci.load_pickle(os.path.join(project_dir,'data'), 'initial_lookup_lists.pkl')

    return(symbol_lookup_list, id_lookup_list)


# Go through the missing lookups in the symbols lookup list and try to determine them using the xref endpoint of the Ensembl REST API
def get_missing_lookups(project_dir, symbol_lookup_list, max_num_names=-1):

    # Import relevant modules
    import os, sys
    gmb_dir = '/data/BIDS-HPC/private/projects/gmb/checkout'
    if gmb_dir not in sys.path:
        sys.path.append(gmb_dir)
    import time_cell_interaction_lib as tci # we need this to get the pickle functions

    # If the file containing the missing lookup list does not already exist...
    if not os.path.exists(os.path.join(project_dir,'data','missing_lookup_list.pkl')):

        # Initialize some overall variables
        missing_lookup_list = []
        iiter = 0
        missing_nidentified = 0

        # For each name in the symbol lookup list...
        imissing = 0
        for item in symbol_lookup_list:

            # If the name is missing an Ensembl ID, try to find it using xref 
            if item[1] is None:
                imissing = imissing + 1
                missing_lookup_list, _, iiter, missing_nidentified = run_and_process_query('xref', [item[0]], missing_lookup_list, iiter, missing_nidentified)
            
            # If we've tried a certain number of items in symbol_lookup_list, stop
            if imissing == max_num_names:
                break

        # Run a check on the number of Nones in the list; otherwise, save the lookup list to disk
        if (sum([ 0 if x[1] is None else 1 for x in missing_lookup_list ]) != missing_nidentified):
            print('ERROR: Inconsistent number of Nones in the lookup list')
            exit
        else:
            tci.make_pickle(missing_lookup_list, os.path.join(project_dir,'data'), 'missing_lookup_list.pkl')

    # Otherwise, read it in
    else:
        missing_lookup_list = tci.load_pickle(os.path.join(project_dir,'data'), 'missing_lookup_list.pkl')

    return(missing_lookup_list)
