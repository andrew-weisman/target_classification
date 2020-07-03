# Prepare the annotation dataframe df_gencode_genes, particularly calculating the exon length of each gene (corresponding to its non-overlapping exons) and adding this as a column to the df_gencode_genes dataframe
# This takes about 10 minutes
def calculate_exon_lengths(gencode_gtf_file, project_dir):

    # Import relevant libraries
    import pandas as pd
    import numpy as np
    import os
    tci = get_tci_library()

    # Set the number of steps to output so we can evaluate progress
    nsteps = 100

    # If the file containing the annotation dataframe does not already exist...
    if not os.path.exists(os.path.join(project_dir,'data','annotation_dataframe.pkl')):

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

        tci.make_pickle(df_gencode_genes, os.path.join(project_dir,'data'), 'annotation_dataframe.pkl')

    # Otherwise, read it in
    else:
        df_gencode_genes = tci.load_pickle(os.path.join(project_dir,'data'), 'annotation_dataframe.pkl')

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
def get_counts(links_dir):
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
def get_intensities(files_per_sample, links_dir, df_gencode_genes, project_dir, nsamples=-1):

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
    if not os.path.exists(os.path.join(project_dir,'data','series_lists.pkl')):

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
        tci.make_pickle([srs_counts, srs_fpkm, srs_fpkm_uq], os.path.join(project_dir,'data'), 'series_lists.pkl')

    # Otherwise, read it in
    else:
        [srs_counts, srs_fpkm, srs_fpkm_uq] = tci.load_pickle(os.path.join(project_dir,'data'), 'series_lists.pkl')

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
    import json

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

    # Define the Pandas dataframe from the fields of interest for all samples
    df = pd.DataFrame(data=selected_fields_per_sample, columns=labels_df_names)
    df = df.set_index('sample id')

    # Return this dataframe
    return(df)


# Plot histograms of the numerical columns of the samples/labels before and after cutoffs could theoretically be applied, and print out a summary of what we should probably do
def demonstrate_removal_of_bad_samples(df_samples, nstd=2):

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
    print('For the time being though, we are leaving the data untouched!')

    # Plot the same histograms using the filtered samples to show what would happen if we applied the calculated cutoffs
    _ = df_samples.iloc[valid_ind,:].hist(figsize=(12,8))
