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
    skip_query = False
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

        # Unfortunately, there can be slashes in synonyms (e.g., OK/SW-cl.4, see http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000181588;r=19:1554669-1568058), but this endpoint (and only this endpoint) will die if there is a slash. Thus, this code dies when searching for the two genes 'OK/SW-CL.36' and 'OK/SW-CL.58'. Thus, we have to skip these. However, they can be identified by placing them in the HTTP URL using the g option, e.g., http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=OK/SW-cl.4. Thus, we should just make a note of the presence of genes with slashes in them and say to manually identify them using this method, manually. Note upon doing this for those two genes, they are not in Ensembl!

        if names_list[0].find('/') != -1:
            print('WARNING: Gene name {} has a "/" in it and is being skipped because otherwise this endpoint ({}) dies. This gene must be queried manually, e.g., using the HTTP method (see the comment in the ensemble_name_request() function in target_class_lib.py). For the time being we\'re recording this gene as having no Ensembl ID.'.format(names_list[0], endpoint))
            skip_query = True

        max_list_len = 1
        #ext = "/xrefs/symbol/homo_sapiens/BRCA2?"
        ext = "/xrefs/symbol/homo_sapiens/" + names_list[0]
        headers = { "Content-Type" : "application/json"}
        data = None

    # Check the names list length
    if len(names_list) > max_list_len:
        print('ERROR: Names list is too long ({}) for the {} endpoint (max length is {})'.format(len(names_list), endpoint, max_list_len))
        exit

    # If we DO want to perform the REST query
    if not skip_query:

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

    # If we DON'T want to perform the REST query, return an empty list
    else:
        return([])


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


# Go through the names in the inputted list of names and try to determine their Ensembl IDs using the xref endpoint of the Ensembl REST API
# This function just saves the results in a pickle file if it doesn't already exist
def get_missing_lookups(pickle_dir, pickle_file, to_xref_single_list, max_num_names=-1):

    # Import relevant modules
    import os, sys
    gmb_dir = '/data/BIDS-HPC/private/projects/gmb/checkout'
    if gmb_dir not in sys.path:
        sys.path.append(gmb_dir)
    import time_cell_interaction_lib as tci # we need this to get the pickle functions

    # If the file containing the missing lookup list does not already exist...
    if not os.path.exists(os.path.join(pickle_dir, pickle_file)):

        print('Calculating and creating file {}'.format(pickle_file))

        # Initialize some overall variables
        missing_lookup_list = []
        iiter = 0
        missing_nidentified = 0

        # For each name in the symbol lookup list...
        imissing = 0
        for name in to_xref_single_list:

            # Try to find the Ensembl ID of the current name using xref 
            imissing = imissing + 1
            missing_lookup_list, _, iiter, missing_nidentified = run_and_process_query('xref', [name], missing_lookup_list, iiter, missing_nidentified)
            
            # If we've tried a certain number of names in to_xref_single_list, stop
            if imissing == max_num_names:
                break

        # Run a check on the number of Nones in the list; otherwise, save the lookup list to disk
        if (sum([ 0 if x[1] is None else 1 for x in missing_lookup_list ]) != missing_nidentified):
            print('ERROR: Inconsistent number of Nones in the lookup list')
            exit
        else:
            tci.make_pickle(missing_lookup_list, pickle_dir, pickle_file)

    else:
        print('Skipping file {} as it already exists'.format(pickle_file))

    return()


# Yield successive n-sized chunks from l
# Taken from https://www.geeksforgeeks.org/break-list-chunks-size-n-python
def divide_chunks(l, n): 
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 


# Go through symbol_lookup_list, split into chunks all the names that don't have Ensembl IDs, and try to determine them if the pickle file corresponding to each chunk doesn't yet exist
def save_missing_lookups_in_chunks(symbol_lookup_list, project_dir, chunk_size=1000):

    # Import relevant module
    import os

    # Get the names in the symbol lookup list that do not have Ensembl IDs
    to_xref = []
    for item in symbol_lookup_list:
        if item[1] is None:
            to_xref.append(item[0])

    # Split this list of names to try xref-ing into chunks of a certain size
    to_xref_lists = list(divide_chunks(to_xref, chunk_size))

    # Set the directory to save the missing lists in and create it if it doesn't already exist
    pickle_dir = os.path.join(project_dir,'data','missing_lookup_lists')
    if not os.path.exists(pickle_dir):
        os.mkdir(pickle_dir)

    # For each sub-list in to_xref, if the corresponding pickle file doesn't already exist, go through the names and try to determine their Ensembl IDs using the xref endpoint of the Ensembl REST API, and save the pickle file
    for ilist, to_xref_single_list in enumerate(to_xref_lists):
        pickle_file = 'missing_lookup_list_{:03d}.pkl'.format(ilist)
        get_missing_lookups(pickle_dir, pickle_file, to_xref_single_list, max_num_names=-1)

    return()
