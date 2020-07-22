# Prepare the annotation dataframe df_gencode_genes, particularly calculating the exon length of each gene (corresponding to its non-overlapping exons) and adding this as a column to the df_gencode_genes dataframe
# This takes about 10 minutes if the pickle file doesn't already exist
def calculate_exon_lengths(gencode_gtf_file):

    # Import relevant libraries
    import pandas as pd
    import numpy as np
    import os
    tci = get_tci_library()

    # Set the number of steps to output so we can evaluate progress
    nsteps = 100

    # Identify the data directory as the directory that the annotation file is in
    os_sep = os.sep
    data_dir = os_sep.join(gencode_gtf_file.split(sep=os_sep)[:-1])

    # If the file containing the annotation dataframe does not already exist...
    if not os.path.exists(os.path.join(data_dir,'annotation_dataframe.pkl')):

        # Read in the the GTF file from the Gencode website
        df_gencode = pd.read_csv(gencode_gtf_file, sep='\t', skiprows=5, header=None)
        df_gencode_genes = df_gencode[df_gencode[2]=='gene'].reset_index(drop=True)
        df_gencode_exons = df_gencode[df_gencode[2]=='exon'].reset_index(drop=True)

        # Format the df_gencode_genes dataframe for consistency
        df_gencode_genes['id'] = df_gencode_genes.apply(lambda x: x[8].split()[1].split('\"')[1], axis=1)
        df_gencode_genes['type'] = df_gencode_genes.apply(lambda x: x[8].split()[3].split('\"')[1], axis=1)
        df_gencode_genes['name'] = df_gencode_genes.apply(lambda x: x[8].split()[7].split('\"')[1], axis=1)
        df_gencode_genes = df_gencode_genes.rename({3: 'start', 4: 'end', 6: 'strand', 0: 'seqname'}, axis='columns')
        df_gencode_genes = df_gencode_genes.set_index('id')
        df_gencode_genes = df_gencode_genes.sort_index()

        # Format the df_gencode_exons dataframe for consistency
        # Takes about a minute
        df_gencode_exons['id'] = df_gencode_exons.apply(lambda x: x[8].split()[1].split('\"')[1], axis=1)
        df_gencode_exons['type'] = df_gencode_exons.apply(lambda x: x[8].split()[3].split('\"')[1], axis=1)
        df_gencode_exons['name'] = df_gencode_exons.apply(lambda x: x[8].split()[7].split('\"')[1], axis=1)
        df_gencode_exons = df_gencode_exons.rename({3: 'start', 4: 'end', 6: 'strand', 0: 'seqname'}, axis='columns')
        df_gencode_exons = df_gencode_exons.set_index('id')
        df_gencode_exons = df_gencode_exons.sort_index()

        # Set the step size in units of the size of the df_gencode_exons dataframe
        unit_len = int(len(df_gencode_exons) / nsteps)

        # Initialize some values
        istep = 0 # the step that we're on
        exon_lengths = [] # the array holding the final exon gene lengths (non-overlapping union of exon base pairs)
        prev_idx = '' # set the previous index to null

        # For every index in the ordered-by-index exons dataframe...
        for iidx, idx in enumerate(df_gencode_exons.index):

            # Get the current row of data in the dataframe
            curr_row = df_gencode_exons.iloc[iidx,:]
            
            # Output progress if the time is right
            if (iidx%unit_len) == 0:
                print('{}/{} complete...'.format(istep,nsteps))
                istep = istep + 1

            # If the current index is not equal to the previous index...
            if idx != prev_idx:

                # If the previous index is not null (i.e., if this isn't the very first loop iteration and therefore base_pairs has been initialized below), calculate and store the number of unique base pairs for the current unique idx
                if prev_idx != '':
                    exon_lengths.append(len(set(np.concatenate(base_pairs))))

                # Initialize the base_pairs holder (which will ultimately be a list of lists of base pairs)
                base_pairs = []

            # Always append the current set of base pairs corresponding to curr_row to the base_pairs list
            base_pairs.append(np.arange(curr_row['start'], curr_row['end']+1))

            # Set the previous index to the current index
            prev_idx = idx

        # Calculate and store the number of unique base pairs for the final unique idx
        exon_lengths.append(len(set(np.concatenate(base_pairs))))

        # Add a column of exon gene length to the genes dataframe
        df_gencode_genes['exon_length'] = exon_lengths

        tci.make_pickle(df_gencode_genes, data_dir, 'annotation_dataframe.pkl')

    # Otherwise, read it in
    else:
        df_gencode_genes = tci.load_pickle(data_dir, 'annotation_dataframe.pkl')

    # Return the main reference dataframe
    return(df_gencode_genes)


# Import the time cell interaction library
def get_tci_library():
    import sys
    gmb_dir = '/data/BIDS-HPC/private/projects/gmb/checkout'
    if gmb_dir not in sys.path:
        sys.path.append(gmb_dir)
    import time_cell_interaction_lib as tci # we need this to get the pickle functions e.g.
    return(tci)


# Get a list of series containing the actual counts of the samples
def get_counts_old(links_dir):
    # Sample call: srs_counts = tc.get_counts('/data/BIDS-HPC/private/projects/dmi2/data/all_gene_expression_files_in_target/links')

    # Import relevant libraries
    import pandas as pd
    import os

    # Define the sample HT-Seq datafiles
    file_counts = 'fffee315-9aa3-44d2-8c89-78a2c1d107e7.htseq_counts.txt'

    # Get a sample list of counts series
    srs_counts = [pd.read_csv(os.path.join(links_dir, file_counts), sep='\t', skipfooter=5, names=['id','intensity'])]

    # Format the sample series for consistency
    for isr, sr in enumerate(srs_counts):
        sr = sr.set_index('id')
        srs_counts[isr] = sr.sort_index().iloc[:,0]

    # Return a list of series of counts
    return(srs_counts)


# Given the counts for a sample, calculate the FPKM and FPKM-UQ values using the protein-coding genes for normalization
def calculate_fpkm(df_gencode_genes, sr_counts):

    # Import relevant library
    import numpy as np

    # Get the number of reads and gene lengths for the entire set of genes
    exon_lengths = df_gencode_genes['exon_length']

    # Get the number of reads and gene lengths for just the protein-coding genes
    pc_loc = df_gencode_genes['type'] == 'protein_coding'
    pc_counts = sr_counts[pc_loc]
    #pc_lengths = exon_lengths[pc_loc]

    # Calculate the normalizations for the FPKM and FPKM-UQ values
    pc_frag_count = pc_counts.sum()
    upper_quantile = np.percentile(pc_counts, 75) # equals pc_counts.sort_values()[int(pc_loc.sum()*.75)]

    # Calculate the normalized counts via https://github.com/NCI-GDC/htseq-tool/blob/master/htseq_tools/tools/fpkm.py
    tmp = sr_counts / exon_lengths * 1e9
    fpkm = tmp / pc_frag_count
    fpkm_uq = tmp / upper_quantile

    # Return the normalized counts
    return(fpkm, fpkm_uq)


# Run a few checks on some known data
def run_checks(df_gencode_genes, calculated_counts, fpkm, fpkm_uq):

    # Sample call: tc.run_checks(df_gencode_genes, srs_counts[0], fpkm, fpkm_uq)

    # Import relevant libraries
    import pandas as pd
    import os

    # Calculate the constants
    gdc_tsv_file = '/data/BIDS-HPC/private/projects/dmi2/data/gencode.gene.info.v22.tsv'
    file_fpkm = 'fffee315-9aa3-44d2-8c89-78a2c1d107e7.FPKM.txt'
    file_fpkm_uq = 'fffee315-9aa3-44d2-8c89-78a2c1d107e7.FPKM-UQ.txt'
    links_dir = '/data/BIDS-HPC/private/projects/dmi2/data/all_gene_expression_files_in_target/links'

    # Read in and process the final TSV file that GDC uses and contains exon lengths whose results we want to check against (our calculated values in df_gencode_genes)
    df_gdc = pd.read_csv(gdc_tsv_file, sep='\t')
    df_gdc = df_gdc.rename({'gene_id': 'id', 'gene_name': 'name', 'gene_type': 'type'}, axis='columns')
    df_gdc = df_gdc.set_index('id')
    df_gdc = df_gdc.sort_index()

    # Read in and process the files to check our calculated results against (known FPKM and FPKM-UQ values)
    df_fpkm = pd.read_csv(os.path.join(links_dir, file_fpkm), sep='\t', names=['id','intensity'])
    df_fpkm_uq = pd.read_csv(os.path.join(links_dir, file_fpkm_uq), sep='\t', names=['id','intensity'])
    srs_known = [df_fpkm, df_fpkm_uq]
    for idf, df in enumerate(srs_known):
        df = df.set_index('id')
        srs_known[idf] = df.sort_index().iloc[:,0]

    # Check for column equality between the two reference datafiles
    for colname in ['name', 'seqname', 'start', 'end', 'strand', 'type']:
        print('Columns equal between the 2 reference files?', df_gdc[colname].equals(df_gencode_genes[colname]))

    # Check that the ID columns of all five dataframes are exactly the same
    df_samples = [calculated_counts] + srs_known
    dfs = df_samples + [df_gdc, df_gencode_genes]
    ndfs = len(dfs)
    import numpy as np
    for idf1 in range(ndfs-1):
        for idf2 in np.array(range(ndfs-1-idf1)) + idf1+1:
            df1 = dfs[idf1]
            df2 = dfs[idf2]
            print('ID columns the same in all 5 dataframes?', idf1, idf2, df1.index.equals(df2.index))

    # Show that we've reproduced what GDC calls the "exon_length" and what I'm assuming is probably the "aggregate_length" as well
    print('Correct calculation of exon_/aggregate_length?', df_gencode_genes['exon_length'].equals(df_gdc['exon_length']))

    # Show that using these exon lengths we have achieved adjusted counts that are proportional to the FPKM values
    tmp = df_samples[0] / df_gencode_genes['exon_length'] / df_samples[1]
    tmp = tmp[tmp.notnull()]
    print('Adjusted counts using the calculated exon lengths proportional to the FPKM values?', tmp.std()/tmp.mean()*100, (tmp-tmp.mean()).abs().max()/tmp.mean()*100)

    # Print how well I reproduced the normalized values that I downloaded from the GDC data portal
    print('Maximum percent error in FPKM: {}'.format((fpkm-df_samples[1]).abs().max() / df_samples[1].mean() * 100))
    print('Maximum percent error in FPKM-UQ: {}'.format((fpkm_uq-df_samples[2]).abs().max() / df_samples[2].mean() * 100))


# Get a list of text files available for each sample in the links_dir
# This is no longer used
def get_files_per_sample(links_dir):
    # Sample call: files_per_sample = tc.get_files_per_sample('/data/BIDS-HPC/private/projects/dmi2/data/all_gene_expression_files_in_target/links')

    # Import relevant libraries
    import glob, os

    # Get a list of all files (with pathnames removed) in links_dir, except for the manifest file
    txt_files = set([ x.split('/')[-1] for x in glob.glob(os.path.join(links_dir,'*')) ]) - {'MANIFEST.txt'}

    # Get the corresponding sorted set of basenames (ostensibly, the unique sample names) from the file list
    basenames = sorted(set([ x.split('.')[0] for x in txt_files ]))

    # For each sample name, create a list of files having that sample name in the filename
    files_per_sample = []
    for basename in basenames:
        files_per_sample.append([ x for x in txt_files if basename in x ])

    # Return a list of text files available for each sample in the links_dir
    return(files_per_sample)


# Read in all the counts files and calculate the FPKM and FPKM-UQ values from them, checking the FPKM/FPKM-UQ values with known quantities if they're available
def get_intensities_old(files_per_sample, links_dir, df_gencode_genes, data_dir, nsamples=-1):

    # Import relevant libraries
    import pandas as pd
    import os
    tci = get_tci_library()

    # Constants (suffixes of the different file types in the links directory)
    # ['htseq.counts', 'htseq_counts.txt'] # these are the non-specified suffixes below
    star_counts_suffix = 'rna_seq.star_gene_counts.tsv'
    fpkm_suffix = 'FPKM.txt'
    fpkm_uq_suffix = 'FPKM-UQ.txt'

    # If the file containing the lists of series does not already exist...
    if not os.path.exists(os.path.join(data_dir,'series_lists.pkl')):

        # For the first namples_to_process samples in files_per_sample...
        srs_counts = []
        srs_fpkm = []
        srs_fpkm_uq = []
        nsamples = ( len(files_per_sample) if nsamples==-1 else nsamples )
        #for isample, files in enumerate(files_per_sample[:nsamples_to_process]):
        for isample, files in enumerate(files_per_sample[:nsamples]):

            # Initialize the descriptions of the sample (filenames namely) that we want to calculate
            counts_file = None
            fpkm_file = None
            fpkm_uq_file = None
            counts_type = None

            # For each file in the file list for the current sample...
            for ifile, x in enumerate([ curr_file.split('.')[1:] for curr_file in files ]):

                # Determine the suffix of the file
                suffix = '.'.join(x)

                # Run logic based on what the suffix of the current file is, calculating the scriptions of the sample (filenames namely) that we want to calculate
                if suffix == fpkm_suffix:
                    fpkm_file = files[ifile]
                elif suffix == fpkm_uq_suffix:
                    fpkm_uq_file = files[ifile]
                else:
                    if suffix == star_counts_suffix:
                        counts_type = 'STAR'
                    else:
                        counts_type = 'HTSeq'
                    counts_file = files[ifile]

            # Print the determined filenames and count filetype for the current sample
            # print('----')
            # print('Counts file ({}): {}'.format(counts_type, counts_file))
            # print('FPKM file: {}'.format(fpkm_file))
            # print('FPKM-UQ file: {}'.format(fpkm_uq_file))

            # Get counts dataframe for the current sample
            if counts_type == 'HTSeq':
                df_tmp = pd.read_csv(os.path.join(links_dir, counts_file), sep='\t', skipfooter=5, names=['id','intensity'])
            else:
                df_tmp = pd.read_csv(os.path.join(links_dir, counts_file), sep='\t', skiprows=5, usecols=[0,1], names=['id','intensity'])

            # Format the counts series and calculate FPKM and FPKM-UQ from it using the aggregate lengths in df_gencode_genes
            sr_counts = df_tmp.set_index('id').sort_index().iloc[:,0]
            sr_fpkm, sr_fpkm_uq = calculate_fpkm(df_gencode_genes, sr_counts)

            # Print how well I reproduced the FPKM values that I downloaded from the GDC data portal, if present
            if fpkm_file is not None:
                df_fpkm = pd.read_csv(os.path.join(links_dir, fpkm_file), sep='\t', names=['id','intensity'])
                sr_fpkm_known = df_fpkm.set_index('id').sort_index().iloc[:,0]
                perc_err = (sr_fpkm-sr_fpkm_known).abs().max() / sr_fpkm_known.mean() * 100
                #print('Maximum percent error in FPKM: {}'.format(perc_err))
                if perc_err > 1e-2:
                    print('ERROR: Maximum percent error ({}) in FPKM is too high!'.format(perc_err))
                    exit()

            # Print how well I reproduced the FPKM-UQ values that I downloaded from the GDC data portal, if present
            if fpkm_uq_file is not None:
                df_fpkm_uq = pd.read_csv(os.path.join(links_dir, fpkm_uq_file), sep='\t', names=['id','intensity'])
                sr_fpkm_uq_known = df_fpkm_uq.set_index('id').sort_index().iloc[:,0]
                perc_err = (sr_fpkm_uq-sr_fpkm_uq_known).abs().max() / sr_fpkm_uq_known.mean() * 100
                #print('Maximum percent error in FPKM-UQ: {}'.format(perc_err))
                if perc_err > 1e-5:
                    print('ERROR: Maximum percent error ({}) in FPKM-UQ is too high!'.format(perc_err))
                    exit()

            # Append the current calculated series to the lists of series
            srs_counts.append(sr_counts)
            srs_fpkm.append(sr_fpkm)
            srs_fpkm_uq.append(sr_fpkm_uq)

            print('\r', '{:3.1f}% complete...'.format((isample+1)/nsamples*100), end='')

        # Write a pickle file containing the data that take a while to calculate
        tci.make_pickle([srs_counts, srs_fpkm, srs_fpkm_uq], data_dir, 'series_lists.pkl')

    # Otherwise, read it in
    else:
        [srs_counts, srs_fpkm, srs_fpkm_uq] = tci.load_pickle(data_dir, 'series_lists.pkl')

    # Return the calculated lists of series
    return(srs_counts, srs_fpkm, srs_fpkm_uq)


# Obtain a Pandas dataframe from the fields of interest for all samples, essentially containing everything we'll ever need to know about the samples, including the labels themselves
def get_labels_dataframe(sample_sheet_file, metadata_file):

    # Sample call:
    #   sample_sheet_file = '/data/BIDS-HPC/private/projects/dmi2/data/gdc_sample_sheet.2020-07-02.tsv'
    #   metadata_file = '/data/BIDS-HPC/private/projects/dmi2/data/metadata.cart.2020-07-02.json'
    #   df_samples = tc.get_labels_dataframe(sample_sheet_file, metadata_file)

    # Import relevant libraries
    import pandas as pd
    import json, os

    # Constants
    htseq_suffixes = ['htseq.counts', 'htseq_counts.txt']
    labels_df_names = ['sample id', 'file list index', 'counts file name', 'average base quality', 'file id', 'project id', 'case id', 'sample type', 'contamination_error', 'proportion_reads_mapped', 'proportion_reads_duplicated', 'contamination', 'proportion_base_mismatch', 'state', 'platform', 'average_read_length', 'entity_submitter_id']
    desired_keys = ['average_base_quality', 'contamination_error', 'proportion_reads_mapped', 'proportion_reads_duplicated', 'contamination', 'proportion_base_mismatch', 'state', 'platform', 'average_read_length']

    # Read in the two datafiles
    df_samples = pd.read_csv(sample_sheet_file, sep='\t')
    with open(metadata_file) as f:
        metadata = json.load(f)

    # Get the corresponding filename mapping arrays
    filename_mapping_samples = df_samples['File Name']
    filename_mapping_metadata = []
    for curr_file in metadata:
        filename_mapping_metadata.append(curr_file['file_name'])

    # Get the full list of unique sample IDs
    samples = list(set(df_samples['Sample ID']))

    # For each unique sample ID...
    selected_fields_per_sample = []
    for sample in samples:

        # Store the dataframe for just the current sample
        df_sample = df_samples[df_samples['Sample ID']==sample]

        # Check that all relevant columns of the current sample are equal, as they should be since they all correspond to the same sample even though the rows correspond to different datafiles
        non_unique_values = (len(df_sample['Data Category'].unique())!=1) or (len(df_sample['Data Type'].unique())!=1) or (len(df_sample['Project ID'].unique())!=1) or (len(df_sample['Case ID'].unique())!=1) or (len(df_sample['Sample ID'].unique())!=1) or (len(df_sample['Sample Type'].unique())!=1)
        if non_unique_values:
            print('ERROR: All fields for the current sample are not equal over the files')
            exit()

        # Obtain the HTSeq counts files for the current sample (note: sometimes there are 2 instead of 1)
        htseq_files = [ (fn if ('.'.join(fn.split('.')[1:-1]) in htseq_suffixes) else None) for fn in df_sample['File Name'] ]
        counts_file_list = list(set(htseq_files) - {None})

        # If there is just a single HTSeq counts files for the current sample, then we have already identified the HTSeq counts file that we're looking for
        if len(counts_file_list) == 1:
            counts_file = counts_file_list[0]

        # If there are more than one HTSeq counts files for the current sample, it doesn't make sense to me to keep multiple analyses of the same sample, so just choose the analysis with the best base quality score (or the first if multiple files have the same best score)
        else:
           
            # Initialize variables in this small loop over possible HTSeq counts files
            best_score = -1
            best_counts_file = None

            # For each counts file...
            for counts_file in counts_file_list:

                # Obtain the corresponding index of the metadata list, using that to extract the average_base_quality field
                metadata_index = filename_mapping_metadata.index(counts_file)
                score = metadata[metadata_index]['analysis']['input_files'][0]['average_base_quality']

                # If the current score is better than the current best score, update the variables
                if score > best_score:
                    best_counts_file = counts_file
                    best_score = score

            # Rename the variable holding the HTSeq counts file that we were looking for
            counts_file = best_counts_file

        # For the best HTSeq counts file for the current sample, obtain the corresponding index of in the sample sheet and in the metadata, and check that they are the same
        samples_index = filename_mapping_samples[filename_mapping_samples==counts_file].index[0]
        metadata_index = filename_mapping_metadata.index(counts_file)
        if samples_index != metadata_index:
            print('ERROR: The File indexes for the sample sheet and the metadata file are different')
            exit()

        # Get shortcut variables for the sections of the sample sheet and the metadata that have some values that we want to add to our labels dataframe
        series = df_samples.loc[samples_index,:]
        md1 = metadata[metadata_index]['analysis']['input_files'][0]

        # Run a check to ensure at least a None value of all desired keys exists in the md1 dictionary
        current_keys = md1.keys()
        for desired_key in desired_keys:
            if desired_key not in current_keys:
                md1[desired_key] = None

        # Save the fields of interest to a running list
        selected_fields_per_sample.append([sample, samples_index, counts_file, md1['average_base_quality'], series['File ID'], series['Project ID'], series['Case ID'], series['Sample Type'], md1['contamination_error'], md1['proportion_reads_mapped'], md1['proportion_reads_duplicated'], md1['contamination'], md1['proportion_base_mismatch'], md1['state'], md1['platform'], md1['average_read_length'], metadata[metadata_index]['associated_entities'][0]['entity_submitter_id']])

    # Define the Pandas dataframe from the fields of interest for all samples, finishing by sorting by index
    df = pd.DataFrame(data=selected_fields_per_sample, columns=labels_df_names)
    df = df.set_index('sample id').sort_index()

    # Return this dataframe
    return(df)


# Plot histograms of the numerical columns of the samples/labels before and after cutoffs could theoretically be applied, and print out a summary of what we should probably do
def remove_bad_samples(df_samples, nstd=2):

    # Import relevant library
    import numpy as np

    # Generate the initial set of histograms on the numerical data in the samples dataframe
    ax_hist = df_samples.hist(figsize=(12,8))

    # "Constants" based on viewing the first set of histograms above so that we can figure out which ones to use to filter the data
    columns = ['average base quality', 'proportion_base_mismatch', 'proportion_reads_mapped']
    higher_is_better = [True, False, True]
    sp_locs = [(0,0), (1,1), (2,1)]

    # Start off filtering none of the data
    valid_ind = np.full((len(df_samples),), True)

    # For each plot containing data we'd like to use to filter our samples...
    for col, hib, sp_loc in zip(columns, higher_is_better, sp_locs):

        # Determine +1 or -1 depending on whether higher is better (if higher is better, use -1)
        sign = -2*int(hib) + 1

        # Get the data values of the current plot
        vals = df_samples[col]

        # Calculate the cutoff for the current plot using the inputted number of standard deviations from the mean as the cutoff
        cutoff = vals.mean() + sign*nstd*vals.std()

        # Determine the current axis in the overall histogram plot
        ax = ax_hist[sp_loc]

        # Determine the y limits of that plot
        ylim = ax.get_ylim()

        # Plot the calculated cutoff as a vertical red line
        ax.plot([cutoff,cutoff], ylim, 'r')

        # Determine a boolean array of where the values fall outside the cutoffs and print how many such bad values there are
        bad_vals = (sign*vals) > (sign*cutoff)
        print('There are {} bad values in the "{}" plot'.format(sum(bad_vals), col))

        # Update the boolean filtering array using the current set of bad values
        valid_ind[bad_vals]=False

    # Store some numbers for easier calculation below: the total number of samples and the aggregate number of bad samples
    ntot = len(valid_ind)
    nbad_tot = ntot - sum(valid_ind)

    # Print some output of the analysis and what we'd recommend we do in the future
    print('Most bad values are overlapping; taken together, there are {} bad values'.format(nbad_tot))
    print('We should likely use these cutoffs to remove the bad samples; this will only remove {:3.1f}% of the data, leaving {} good samples'.format(nbad_tot/ntot*100, ntot-nbad_tot))
    print('See for example the two generated images: the first is the original data with the cutoffs plotted in red, and the second is the filtered data with the cutoffs applied')
    # print('For the time being though, we are leaving the data untouched!')

    # Plot the same histograms using the filtered samples to show what would happen if we applied the calculated cutoffs
    _ = df_samples.iloc[valid_ind,:].hist(figsize=(12,8))

    return(df_samples.iloc[valid_ind,:], valid_ind)


# Return the samples dataframe with the samples removed that correspond to multiple cases (i.e., people)
def drop_multiperson_samples(df_samples):

    # Import relevant libraries
    import pandas as pd
    import numpy as np

    # Initialize the arrays of interest
    indexes_to_drop = [] # to store indexes of samples to drop
    samples_to_drop = [] # to store the samples themselves (Series) to drop
    indexes_to_keep = np.full((len(df_samples),), True) # to store the indexes to keep in the full samples array so that I can plug these logical indexes into other arrays

    # For every index in the sample dataframe...
    for isample, sample_index in enumerate(df_samples.index):

        # Save the current sample (as a Series)
        sample = df_samples.loc[sample_index]

        # If there are multiple cases (people; and I've confirmed that they can't be the same people, e.g., living white female and dead black male) corresponding to this one sample...
        if len(sample['case id'].split()) > 1:
            indexes_to_drop.append(sample_index) # save the sample index
            samples_to_drop.append(sample) # save the sample series
            indexes_to_keep[isample] = False

    # Create and print a Pandas dataframe of the samples to drop in order to visualize it nicely
    print('Dropping the following samples from the samples table:')
    df_samples_to_drop = pd.DataFrame(data=samples_to_drop).rename_axis(index='sample id')
    print(df_samples_to_drop)

    # Return the modified samples dataframe
    return(df_samples.drop(index=indexes_to_drop), indexes_to_keep, df_samples_to_drop)


# Perform exploratory data analysis on the sample labels
def eda_labels(df_samples):

    # Import relevant library
    import random

    # Add the index "column" as an actual column to the dataframe so we can analyze the index column in the same manner as the other columns
    df_samples[df_samples.index.name] = df_samples.index

    # Initialize the holder lists of the column types
    cols_unique = []
    cols_uniform = []
    cols_other = []

    # Get the total number of rows in the dataframe
    nsamples = len(df_samples)

    # Get a random index in the dataframe
    rand_index = random.randrange(nsamples)

    # Plot histograms of the numeric data
    _ = df_samples.hist(figsize=(12,8))

    # Determine the non-numeric columns
    non_numeric_cols = df_samples.select_dtypes(exclude='number').columns

    # Initialize the column name lengths
    max_col_len = -1

    # For every non-numeric column...
    for col in non_numeric_cols:

        # Determine the number of unique values in the column
        nunique = len(df_samples[col].unique())

        # Every row in the column is unique
        if nunique == nsamples:
            cols_unique.append([col, df_samples[col][rand_index]])

        # The column is completely uniform
        elif nunique == 1:
            cols_uniform.append([col, df_samples[col][0]])

        # The column is neither unique nor uniform
        else:
            cols_other.append([col, nunique])

        # Possibly update the maximum column name size (for pretty printing later)
        if len(col) > max_col_len:
            max_col_len = len(col)

    # Store the output format string depending on the maximum column name size
    output_col_str = ' . {:' + str(max_col_len+5) + '}{}'

    # Print the columns in which every row is unique
    print('Non-numeric columns with all unique values ({} of them), with sample values:\n'.format(nsamples))
    for col_data in cols_unique:
        print(output_col_str.format(col_data[0], col_data[1]))

    # Print the columns that are completely uniform
    print('\nNon-numeric columns with uniform values:\n')
    for col_data in cols_uniform:
        print(output_col_str.format(col_data[0], col_data[1]))

    # Print the columns (and supporting information) that are neither unique nor uniform
    print('\nNon-numeric columns with non-unique and non-uniform values:\n')
    for col_data in cols_other:
        print(output_col_str.format(col_data[0], col_data[1]), '\n')
        print(df_samples[col_data[0]].value_counts(), '\n')


# Read in the counts for all the samples in the samples dataframe df_samples
# the counts dataframe will be in the same order as df_samples
def get_counts(df_samples, links_dir):

    # Import relevant libraries
    import os
    import pandas as pd

    # Ensure that all values in the "counts file name" column of df_samples are unique as expected
    nsamples = len(df_samples)
    if not len(df_samples['counts file name'].unique()) == nsamples:
        print('ERROR: "counts file name" column of the samples dataframe does not contain all unique values')
        exit()

    # Strip the ".gz" off of the filenames in the "counts file name" column of the samples dataframe
    counts_filenames = [ x.split(sep='.gz')[0] for x in df_samples['counts file name'] ]

    # For every counts filename in the samples dataframe...
    srs_counts = []
    for isample, counts_fn in enumerate(counts_filenames):

        # Read in the counts data
        sr_counts = pd.read_csv(os.path.join(links_dir, counts_fn), sep='\t', skipfooter=5, names=['id','intensity']).set_index('id').sort_index().iloc[:,0] # assume this is of the HTSeq (as opposed to STAR) format
        
        # Append the read-in and calculated values to running lists
        srs_counts.append(sr_counts)

        print('\r', '{:3.1f}% complete...'.format((isample+1)/nsamples*100), end='')
    
    # Put the list of series into a Pandas dataframe
    df_counts = pd.DataFrame(srs_counts, index=df_samples.index)

    # Return the calculated lists of series
    return(df_counts)


# Convert the lists of Pandas series to Pandas dataframes
def make_intensities_dataframes(srs_list, index):
    import pandas as pd
    counts_list = []
    for srs in srs_list:
        counts_list.append(pd.DataFrame(srs, index=index))
    return(counts_list)


# Print some random data for us to spot-check in the files themselves to manually ensure we have a handle on the data arrays
def spot_check_data(df_samples, df_counts, df_fpkm, df_fpkm_uq, nsamples=4):

    # Import relevant library
    import random

    # Constants
    intensity_types = ['counts', 'FPKM', 'FPKM-UQ']

    # Variable
    intensities = [df_counts, df_fpkm, df_fpkm_uq]

    # Get some values from the intensity data
    nsamples_tot = intensities[0].shape[0]
    sample_names = intensities[0].index

    # For each intensity type...
    for iintensity, intensity_type in enumerate(intensity_types):

        # For each of nsamples random samples in the data...
        for sample_index in random.sample(range(nsamples_tot), k=nsamples):

            # Store the current sample name
            sample_name = sample_names[sample_index]

            # Get the non-zero intensities for the current sample
            srs = intensities[iintensity].iloc[sample_index,:]
            srs2 = srs[srs!=0]

            # Get a random index of the non-zero intensities and store the corresponding intensity and gene
            srs2_index = random.randrange(len(srs2))
            intensity = srs2[srs2_index]
            gene = srs2.index[srs2_index]

            # Get some important data from the samples dataframe
            project_id = df_samples.loc[sample_name, 'project id']
            sample_type = df_samples.iloc[sample_index, 6]

            # Print what we should see in the files
            print('Sample {} ({}, {}) should have a {} value of {} for gene {}'.format(sample_name, project_id, sample_type, intensity_type, intensity, gene))


# Load the data downloaded from the GDC Data Portal
def load_gdc_data(sample_sheet_file, metadata_file, links_dir):

    # Import the relevant libraries
    import os
    tci = get_tci_library()

    # Identify the data directory as the directory that the sample sheet file is in
    os_sep = os.sep
    data_dir = os_sep.join(sample_sheet_file.split(sep=os_sep)[:-1])

    # If the file containing the GDC data does not already exist...
    if not os.path.exists(os.path.join(data_dir,'gdc_data.pkl')):

        # Obtain a Pandas dataframe from the fields of interest for all samples, essentially containing everything we'll ever need to know about the samples, including the labels themselves
        df_samples = get_labels_dataframe(sample_sheet_file, metadata_file) # this will always be in alphabetical order of the sample IDs

        # Read in the counts for all the samples in the samples dataframe df_samples
        df_counts = get_counts(df_samples, links_dir) # the counts dataframe will be in the same order as df_samples

        # Write a pickle file containing the data that takes a while to calculate
        tci.make_pickle([df_samples, df_counts], data_dir, 'gdc_data.pkl')

    # Otherwise, read it in
    else:
        [df_samples, df_counts] = tci.load_pickle(data_dir, 'gdc_data.pkl')

    return(df_samples, df_counts)


# Calculate the FPKM and FPKM-UQ dataframes, and check them with known values if the needed datafiles are present
# the FPKM and FPKM-UQ dataframes will be in the same order as df_samples
# since df_counts is in the same order as df_samples, then the counts and FPKM/FPKM-UQ dataframes will also be aligned
# regardless, this shouldn't be a problem moving forward, since df_samples will always be in lexical order of the sample IDs!
# not to mention, whenever we're unsure, we should run the spot-checks!
def get_fpkm(df_counts, annotation_file, df_samples, links_dir):

    # Import relevant libraries
    import os
    import pandas as pd
    tci = get_tci_library()

    # Identify the data directory as the directory that the annotation file is in
    os_sep = os.sep
    data_dir = os_sep.join(annotation_file.split(sep=os_sep)[:-1])

    # If the file containing the FPKM/FPKM-UQ data does not already exist...
    if not os.path.exists(os.path.join(data_dir,'fpkm_data.pkl')):

        # Ensure that all values in the "counts file name" column of df_samples are unique as expected
        nsamples = len(df_samples)

        # Strip the ".gz" off of the filenames in the "counts file name" column of the samples dataframe
        counts_filenames = [ x.split(sep='.gz')[0] for x in df_samples['counts file name'] ]

        # Obtain a listing of all the files in the links directory
        files_in_links_dir = os.listdir(links_dir)

        # Prepare the annotation dataframe df_gencode_genes, particularly calculating the exon length of each gene (corresponding to its non-overlapping exons) and adding this as a column to the df_gencode_genes dataframe
        # This takes about 10 minutes if the pickle file doesn't already exist
        df_gencode_genes = calculate_exon_lengths(annotation_file)

        # For every counts filename in the samples dataframe...
        srs_fpkm = []
        srs_fpkm_uq = []
        for isample, counts_fn in enumerate(counts_filenames):

            # Read in the counts data
            sr_counts = df_counts.iloc[isample,:]

            # Use those counts data to calculate the FPKM and FPKM-UQ values
            sr_fpkm, sr_fpkm_uq = calculate_fpkm(df_gencode_genes, sr_counts)

            # Get the basename of the current counts file
            bn = counts_fn.split(sep='.')[0]

            # Determine the files in files_in_links_dir and their indexes matching the current basename
            bn_matches = []
            bn_matches_indexes = []
            for ifile, curr_file in enumerate(files_in_links_dir):
                if bn in curr_file:
                    bn_matches.append(curr_file)
                    bn_matches_indexes.append(ifile)

            # From the matching files, determine their suffixes in lowercase, finding where in them FPKM and FPKM-UQ strings match
            suffixes_lower = [ x.split(sep=bn+'.')[1].lower() for x in bn_matches ]
            fpkm_matches = [ 'fpkm.' in x for x in suffixes_lower ]
            fpkm_uq_matches = [ 'fpkm-uq.' in x for x in suffixes_lower ]

            # Ensure there aren't more than 1 match for either FPKM or FPKM-UQ for the current basename
            if sum(fpkm_matches)>1 or sum(fpkm_uq_matches)>1:
                print('ERROR: More than 1 FPKM or FPKM-UQ file matches the basename {}'.format(bn))
                exit()

            # If an FPKM file corresponding to the current basename is found...
            if sum(fpkm_matches) == 1:

                # Determine its filename
                fpkm_fn = files_in_links_dir[bn_matches_indexes[fpkm_matches.index(True)]]

                # Read in its data into a Pandas series
                sr_fpkm_known = pd.read_csv(os.path.join(links_dir, fpkm_fn), sep='\t', names=['id','intensity']).set_index('id').sort_index().iloc[:,0]

                # Determine how well our calculated values in sr_fpkm match those read in to sr_fpkm_known
                perc_err = (sr_fpkm-sr_fpkm_known).abs().max() / sr_fpkm_known.mean() * 100
                if perc_err > 1e-2:
                    print('ERROR: Maximum percent error ({}) in FPKM is too high!'.format(perc_err))
                    exit()

            # If an FPKM-UQ file corresponding to the current basename is found...
            if sum(fpkm_uq_matches) == 1:

                # Determine its filename
                fpkm_uq_fn = files_in_links_dir[bn_matches_indexes[fpkm_uq_matches.index(True)]]

                # Read in its data into a Pandas series
                sr_fpkm_uq_known = pd.read_csv(os.path.join(links_dir, fpkm_uq_fn), sep='\t', names=['id','intensity']).set_index('id').sort_index().iloc[:,0]

                # Determine how well our calculated values in sr_fpkm_uq match those read in to sr_fpkm_uq_known
                perc_err = (sr_fpkm_uq-sr_fpkm_uq_known).abs().max() / sr_fpkm_uq_known.mean() * 100
                if perc_err > 1e-5:
                    print('ERROR: Maximum percent error ({}) in FPKM-UQ is too high!'.format(perc_err))
                    exit()

            # Append the read-in and calculated values to running lists
            srs_fpkm.append(sr_fpkm)
            srs_fpkm_uq.append(sr_fpkm_uq)

            print('\r', '{:3.1f}% complete...'.format((isample+1)/nsamples*100), end='')

        # Put the lists of series into dataframes
        df_fpkm = pd.DataFrame(srs_fpkm, index=df_samples.index)
        df_fpkm_uq = pd.DataFrame(srs_fpkm_uq, index=df_samples.index)

        # Write a pickle file containing the data that takes a while to calculate
        tci.make_pickle([df_fpkm, df_fpkm_uq], data_dir, 'fpkm_data.pkl')

    # Otherwise, read it in
    else:
        [df_fpkm, df_fpkm_uq] = tci.load_pickle(data_dir, 'fpkm_data.pkl')

    return(df_fpkm, df_fpkm_uq)


# Calculate the TPM using the counts and gene lengths
def get_tpm(C_df, annotation_file):

    # Note: I've confirmed TPM calculation with get_tpm_from_fpkm() function below using both FPKM and FPKM-UQ via:
    # df_tpm = tc.get_tpm(df_counts, annotation_file)
    # df_tpm1 = tc.get_tpm_from_fpkm(df_fpkm)
    # df_tpm2 = tc.get_tpm_from_fpkm(df_fpkm_uq)
    # import numpy as np
    # print(np.amax(np.abs(df_tpm1-df_tpm).to_numpy(), axis=(0,1)))
    # print(np.amax(np.abs(df_tpm2-df_tpm).to_numpy(), axis=(0,1)))
    # print(np.amax(np.abs(df_tpm2-df_tpm1).to_numpy(), axis=(0,1)))
    # print(np.sqrt(np.mean(((df_tpm1-df_tpm)**2).to_numpy(), axis=(0,1))))
    # print(np.sqrt(np.mean(((df_tpm2-df_tpm)**2).to_numpy(), axis=(0,1))))
    # print(np.sqrt(np.mean(((df_tpm2-df_tpm1)**2).to_numpy(), axis=(0,1))))

    # Import relevant library
    import numpy as np

    # Calculate the aggregate exon lengths the way GDC does it
    # series of length ngenes
    L_srs = calculate_exon_lengths(annotation_file)['exon_length']

    # Ensure the gene order in the counts and lengths is consistent so that we can perform joint operations on them
    if not C_df.columns.equals(L_srs.index):
        print('ERROR: Order of genes in the counts dataframe is not the same as that in the lengths series')
        exit()

    # Extract the numbers of samples and genes; it seems like this may be unnecessary as seen in the comment after counts_norm, but doing things explicitly like this is significantly faster
    # C_df is a dataframe of shape (nsamples, ngenes)
    nsamples, ngenes = C_df.shape

    # Normalize the counts by their corresponding gene lengths
    # denominator: (nsamples,ngenes) --> repeats over axis=0 (i.e., depends only on gene, not sample) --> L_ij = L_j
    # numerator: (nsamples,ngenes) --> C_ij
    # Cn_ij = C_ij / L_j
    counts_norm = C_df / np.tile(np.expand_dims(L_srs, axis=0), (nsamples,1)) # this equals "C_df / L_srs" (which is simpler) but doing it this way is significantly faster

    # Calculate the normalization factor for each sample
    # D_ij = SUM(Cn_ij,j) --> repeats over axis=1 (i.e., depends only on sample, not gene) --> D_ij = D_i
    denom = np.tile(np.expand_dims(counts_norm.sum(axis=1), axis=1), (1,ngenes))

    # Calculate the TPM
    # Cn_ij / D_ij * 10^6
    # C_ij / L_j / SUM(C_ij/L_j,j) * 10^6
    # T_si = C_si / L_i / SUM(C_sk/L_k,k) * 10^6
    # This is perfectly consistent with boxed formula in sectino 2.3.1 of tpm_calculation.pdf
    tpm = counts_norm / denom * 1e6

    return(tpm)


# Calculate the TPM using FPKM or FPKM-UQ
def get_tpm_from_fpkm(F_df):

    # Import relevant library
    import numpy as np

    # Extract the number of genes
    # F_ij: (nsamples,ngenes)
    ngenes = F_df.shape[1]

    # Calculate the normalization factor for each sapmle
    # D_ij = SUM(F_ij,j) --> repeats over axis=1 (i.e., depends only on sample, not gene) --> D_ij = D_i
    denom = np.tile(np.expand_dims(F_df.sum(axis=1), axis=1), (1,ngenes))

    # Calculate the TPM
    # T_ij = F_ij / D_ij * 10^6
    # T_si = F_si / SUM(F_sk,k) * 10^6
    # This formula is perfectly consistent with the first lines in sections 2.3.1 and 2.3.2 of tpm_calculation.pdf
    tpm = F_df / denom * 1e6

    return(tpm)


# Write annotation and gene counts files (two files total) that are in the same format as the pasilla example so that we can follow the steps outlined at
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input
def write_sample_for_deseq2_input(srs_labels, df_counts, data_directory, reqd_string_in_label='primary tumor', nsamples_per_condition=[5,3,9]):

    # Sample call: write_sample_for_deseq2_input(df_samples['label 1'], df_counts, data_directory)

    # Import relevant libraries
    import numpy as np
    import os, random

    # Get a subset (using both a string in the condition names and a particular number of conditions) of the series of the value counts of the label of interest
    label_value_counts = srs_labels.value_counts()
    srs_subset = label_value_counts[[ reqd_string_in_label in x.lower() for x in label_value_counts.index ]][:len(nsamples_per_condition)]
    print('Using the following conditions (though not all of the samples for each label):')
    print(srs_subset)

    # Construct a list of indexes (actual numbers) to use as a sample of all our data
    all_indexes_to_use = []
    for label, nsamples in zip(srs_subset.index, nsamples_per_condition): # for each condition (label) and inputted number of samples to use for each condition...
        indexes = np.argwhere((srs_labels==label).to_numpy()).flatten() # get the numerical indexes of the current label
        #indexes_to_use = list(indexes[:nsamples]) # get just the number of numerical indexes that we want for the current condition - first nsamples of the list
        indexes_to_use = random.sample(list(indexes), nsamples) # get just the number of numerical indexes that we want for the current condition - random nsamples of the list
        print('\nHere are the {} indexes out of {} that correspond to the condition {}:'.format(len(indexes), len(srs_labels), label))
        print(indexes)
        #print('However, we\'re only using the first {}:'.format(nsamples))
        print('However, we\'re using just a random sample of {} items:'.format(nsamples))
        print(indexes_to_use)
        all_indexes_to_use = all_indexes_to_use + indexes_to_use
    print('\nHere is the final set of numerical indexes that we\'re using ({}={} of them):'.format(sum(nsamples_per_condition), len(all_indexes_to_use)))
    print(all_indexes_to_use)

    # Get just a sample of the labels/conditions and counts
    all_samples_to_use = srs_labels.index[all_indexes_to_use] # get the actual descriptive indexes from the numerical indexes
    labels_to_use = srs_labels[all_samples_to_use]
    counts_to_use = df_counts.loc[all_samples_to_use,:].transpose()

    # Delete rows of counts that are all zeros
    counts_to_use = counts_to_use[(counts_to_use!=0).any(axis=1)]

    # Do a quick check of the list of labels/conditions
    conditions_list = []
    for nsamples, label in zip(nsamples_per_condition, labels_to_use[np.cumsum(nsamples_per_condition)-1]):
        conditions_list = conditions_list + [label]*nsamples
    if conditions_list != labels_to_use.to_list():
        print('ERROR: The actual list of labels/conditions is not what\'s expected')
        exit()

    # Check that the indexes of the counts and labels that we're going to write out are the same    
    if not counts_to_use.columns.equals(labels_to_use.index):
        print('ERROR: Indexes/columns of the labels/counts are inconsistent')
        exit()

    # Write the annotation file in the same format as the pasilla example
    with open(file=os.path.join(data_directory, 'annotation.csv'), mode='w') as f:
        print('"file","condition"', file=f)
        for curr_file, condition in zip(labels_to_use.index, labels_to_use):
            print('"{}","{}"'.format(curr_file, condition), file=f)

    # Write the gene counts in the same format as the pasilla example
    with open(file=os.path.join(data_directory, 'gene_counts.tsv'), mode='w') as f:
        counts_to_use.to_csv(f, sep='\t', index_label='gene_id')
