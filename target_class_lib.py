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
# Note that the result appears to come out unsorted, so that when doing a self.keys() on the result, it's in a different order than names_list
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
# As long as the README.md steps are followed, everything should be uppercase
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
        #print('Skipping file {} as it already exists'.format(pickle_file))
        missing_lookup_list = tci.load_pickle(pickle_dir, pickle_file)

    return(missing_lookup_list)


# Yield successive n-sized chunks from l
# Taken from https://www.geeksforgeeks.org/break-list-chunks-size-n-python
def divide_chunks(l, n): 
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 


# Go through symbol_lookup_list, split into chunks all the names that don't have Ensembl IDs, and try to determine them if the pickle file corresponding to each chunk doesn't yet exist
def get_missing_lookups_in_chunks(symbol_lookup_list, project_dir, chunk_size=1000, pickle_dir_single='missing_lookup_lists'):

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
    pickle_dir = os.path.join(project_dir,'data',pickle_dir_single)
    if not os.path.exists(pickle_dir):
        os.mkdir(pickle_dir)

    # For each sub-list in to_xref, if the corresponding pickle file doesn't already exist, go through the names and try to determine their Ensembl IDs using the xref endpoint of the Ensembl REST API, and save the pickle file
    missing_lookups = []
    for ilist, to_xref_single_list in enumerate(to_xref_lists):
        pickle_file = 'missing_lookup_list_{:03d}.pkl'.format(ilist)
        missing_lookups.append(get_missing_lookups(pickle_dir, pickle_file, to_xref_single_list, max_num_names=-1))

    return(missing_lookups)


# From the two initial lookup lists and the set of chunks of "missing" lookups, create the full lookup table as a Pandas dataframe
# This contains a row corresponding to every unique name present in the data
# All the values are populated using the three endpoints of the REST API
# Everything should be uppercase
def get_data_lookup_table(id_lookup_list, symbol_lookup_list, missing_lookup_chunks):

    # Import relevant modules
    import pandas as pd
    import numpy as np

    # Combine the two "initial" lookup lists
    data = id_lookup_list + symbol_lookup_list

    # Create the initial Pandas lookup dataframe
    df = pd.DataFrame(data=[x[1] for x in data], index=[x[0] for x in data], columns=['id'])

    # Print the number of nulls/Nones in the initial lookup dataframe
    nnull_initial = df.isnull()['id'].sum()
    print('Out of {} entries in the initial lookup table, {} are null'.format(len(df), nnull_initial)) # this is a good check for the number of Nones in the two inital lists (id_lookup_list and symbol_lookup_list)

    # Create the indexes and values of lists to assign in the lookup table from the missing lookups
    indexes = []
    values = []
    for chunk in missing_lookup_chunks:
        for item in chunk:
            if item[1] is not None:
                #nnotnones = nnotnones + 1 # this was a good check for how many we identified (based on the .out files), so it gives us confidence that we're setting the indexes and values correctly
                indexes.append(item[0])
                values.append(item[1])

    # Add the "missing" entries to the "initial" lookup table
    df.loc[indexes,'id'] = values

    # Print the new number of nulls/Nones in the better-filled-in lookup dataframe after checking that the number added makes sense
    nnull_final = df.isnull()['id'].sum()
    nadded = len(values)
    assert(nnull_initial-nadded==nnull_final)
    print('{} entries have been added to the lookup table, so now {} are null'.format(nadded, nnull_final))

    # Add a label to the index
    df = df.rename_axis(index='name')

    # Return the final lookup table
    return(df)


# Create a lookup table using the output from the Biomart website
# This contains only rows from the HGNC table in which Ensembl IDs have been assigned; nothing pertains to the datafiles themselves
# Everything is uppercase
def get_hgnc_lookup_table(project_dir):

    # Import relevant modules
    import pandas as pd
    import os

    # Read the result.txt file created by the Biomart website
    df = pd.read_csv(os.path.join(project_dir,'data','gene_lookup_table.txt'), sep='\t', names=['hgnc','symbol','id'], header=0)

    # Delete rows without an Ensembl ID and make the contents of the symbol and id columns uppercase
    df = df.dropna(axis='index') # this should result in 38,956 non-header rows --> it does
    df['symbol'] = df['symbol'].apply(lambda x: x.upper())
    df['id'] = df['id'].apply(lambda x: x.upper())

    # Set the index to be the symbol column, delete that column and the HGNC ID column, and rename the index
    df = df.set_index('symbol')
    df = df.drop(columns='hgnc')
    df = df.rename_axis(index='name')

    # Duplicate the table but this time using the Ensembl IDs as the index labels
    df2 = df.set_index(df['id'])
    df2 = df2.rename_axis(index='name')

    # Return the two table stacked into one
    return(pd.concat([df2,df]))


# Test that string index and column values are all uppercase
def test_caps_equality(df):
    for icol, series in enumerate([df.index.to_series()] + [df.iloc[:,icol] for icol in range(df.shape[1])]):
        series2 = series[series.notnull()]
        series3 = series2.apply(lambda x: x.upper())
        print('Uppercase equality of column {}: {}'.format(icol, (series2==series3).sum(axis=0)==len(series2)))


# Incorporate some additional synonyms into the data lookup table using the HGNC database
def incorporate_hgnc_lookups(data_lookup, hgnc_lookup):

    # Get the names in the data lookup table that don't have Ensembl IDs
    set1 = set(data_lookup.index[data_lookup['id'].isnull()]) # the question is which of these are in HGNC lookup

    # Get all names in the HGNC lookup table
    set2 = set(hgnc_lookup.index)

    # Determine which of the former are in the latter
    intersection = set1 & set2

    print('There are {} names in the data lookup table lacking Ensembl IDs that DO have Ensembl IDs in the HGNC lookup table'.format(len(intersection)))
    print('However, not all of the corresponding Ensembl IDs in the HGNC lookup table are necessarily in the Ensembl database; that\'s what we\'re about to find out')
    print(intersection)

    # For each of the names just identified...
    pairs_to_test = []
    for name in intersection:

        # If the name isn't in Ensembl format (which we would have already identified), save the name-ID pair from the HGNC lookup table
        if name != hgnc_lookup.loc[name,'id']:
            pairs_to_test.append([hgnc_lookup.loc[name,'id'], name])

    # Save just the Ensembl IDs from the HGNC lookup table to look up in the Ensembl database
    id_list2 = [x[0] for x in pairs_to_test]

    # Look up in the Ensembl database using the REST API these Ensembl IDs to see if they're present in the current version of the database
    iiter = 0
    lookup_list = []
    nidentified = 0
    lookup_list, _, _, nidentified = run_and_process_query('id', id_list2, lookup_list, iiter, nidentified) # returned 11 as expected

    # Save the results as a dictionary
    tested_dict = dict(pairs_to_test)

    # If actual lookups exist, create the list of lookups to add to the data lookup table as a list
    lookups_to_add = []
    for item in lookup_list:
        if item[1] is not None:
            lookups_to_add.append([tested_dict[item[1]], item[1]])

    # Separate the list into the HGNC symbols and their corresponding Ensembl IDs
    symbols = [x[0] for x in lookups_to_add] # these don't show up in Ensembl
    ids = [x[1] for x in lookups_to_add] # these do

    # Print these entries in the data lookup table before adding them to the lookup table
    print(data_lookup.loc[symbols,'id'])

    # Add the entries to the data lookup table
    data_lookup.loc[symbols,'id'] = ids

    # Print the "afer" version
    print(data_lookup.loc[symbols,'id'])

    # Summarize what we did
    print('We added {} lookups to the data lookup table using the HGNC lookup table; now the number of null entries in the data lookup table is {}'.format(nidentified, data_lookup.isnull()['id'].sum()))

    # Separately, note that there are some duplicate Ensembl IDs in the HGNC database
    counts = hgnc_lookup.loc[:,'id'].value_counts()
    duplicates = list(counts.index[counts!=2])
    print('Note: There are {} (out of {}) duplicate Ensembl IDs in the HGNC database:'.format(len(duplicates), int(len(hgnc_lookup)/2)))
    print(data_lookup.loc[duplicates,'id'])

    return(data_lookup)
