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
def get_intensities(files_per_sample, links_dir, df_gencode_genes, project_dir, nsamples_to_process=-1):

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
        nsamples = ( len(files_per_sample) if nsamples_to_process==-1 else nsamples_to_process )
        for isample, files in enumerate(files_per_sample[:nsamples_to_process]):

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
