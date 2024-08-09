# IMPORT PACKAGES
#%pip install -q malariagen_data
import sys
import pandas as pd
import numpy as np
import dask
import dask.array as da
import xarray as xr
import malariagen_data
import subprocess
import os
import h5py

# don't need to typecast all args - happens in other scripts
WINDOW_WIDTH=sys.argv[1]
SITES_PER_WINDOW=sys.argv[2]
CHUNK=sys.argv[3]
WINDOWS_PER_CHUNK=sys.argv[4]
N=sys.argv[5]
RANDSAMP=sys.argv[6]
group_size=int(sys.argv[7])
sample_query_init=sys.argv[8]

# Define the grouping function
def group_population(df, group_size=group_size):
    grouped_data = []
    
    # Group by 'location' and 'year'
    grouped = df.groupby(['location', 'year'])
    
    for group_name, group_df in grouped:
        location, year = group_name
        # Only consider groups with 50 or more entries
        if len(group_df) >= group_size:
            num_subgroups = len(group_df) // group_size
            for i in range(num_subgroups):
                start_idx = i * group_size
                end_idx = (i + 1) * group_size
                subgroup = group_df.iloc[start_idx:end_idx].copy()
                subgroup['group_id'] = f"{location}-{year}-{i + 1}".replace(" ", "-").replace("(", "-").replace(")", "-")
                grouped_data.append(subgroup)
                
    return pd.concat(grouped_data)

def get_data(sample_query, POP_NAME):
    site_mask="gamb_colu"

    # LOAD METADATA
    metadata=ag3.sample_metadata(sample_query=sample_query)
    print("POP_NAME = " + POP_NAME)
    print("SAMPLE_QUERY = " + sample_query)
    print(metadata)

    # NO site_class TO ACCESS ALL SITES (0-FOLD, 4-FOLD, INTERGENIC, ETC)
    # JUST EXTRACT BIALLELIC SNPS (SO WE DON'T INCLUDE FIXED SITES)
    ds_snp_2R=ag3.biallelic_snp_calls(region="2R",
                    site_mask=site_mask, sample_query=sample_query)

    ds_snp_2L=ag3.biallelic_snp_calls(region="2L",
                    site_mask=site_mask, sample_query=sample_query)

    ds_snp_3R=ag3.biallelic_snp_calls(region="3R",
                    site_mask=site_mask, sample_query=sample_query)

    ds_snp_3L=ag3.biallelic_snp_calls(region="3L",
                    site_mask=site_mask, sample_query=sample_query)


    ds_snp_X=ag3.biallelic_snp_calls(region="X",
                    site_mask=site_mask, sample_query=sample_query)

    #Compute Genotype and position
    genotype_2R=ds_snp_2R['call_genotype'].compute()
    POS_2R=ds_snp_2R['variant_position'].compute()

    genotype_2L=ds_snp_2L['call_genotype'].compute()
    POS_2L=ds_snp_2L['variant_position'].compute()

    genotype_3R=ds_snp_3R['call_genotype'].compute()
    POS_3R=ds_snp_3R['variant_position'].compute()

    genotype_3L=ds_snp_3L['call_genotype'].compute()
    POS_3L=ds_snp_3L['variant_position'].compute()

    genotype_X=ds_snp_X['call_genotype'].compute()
    POS_X=ds_snp_X['variant_position'].compute()

    sample_id=metadata['sample_id']
    sample_id=list([bytes(id, "utf-8") for id in sample_id])
    sample_id

    # EXPORT THE GENOTYPES AND POSITIONS TO .hd5 SO THAT I CAN HANDLE IT WITH R
    outfile = '../data/' + POP_NAME + '.hdf5'
    f=h5py.File(outfile, 'a')
    f.create_dataset(name='/POS_2R', data=POS_2R)
    f.create_dataset(name='/genotype_2R',data=genotype_2R)

    f.create_dataset(name='/POS_2L', data=POS_2L)
    f.create_dataset(name='/genotype_2L',data=genotype_2L)

    f.create_dataset(name='/POS_3R', data=POS_3R)
    f.create_dataset(name='/genotype_3R',data=genotype_3R)

    f.create_dataset(name='/POS_3L', data=POS_3L)
    f.create_dataset(name='/genotype_3L',data=genotype_3L)

    f.create_dataset(name='/POS_X', data=POS_X)
    f.create_dataset(name='/genotype_X',data=genotype_X)\

    f.create_dataset(name='/sample_id', data=sample_id)
    f.close()

POP_GROUPS_FILE = "../data/1-METADATA_POPULATIONS.csv"
POP_GROUPS_FILE_SUMMARY = "../data/1-METADATA_POPULATIONS_SUMMARY.csv"
if os.path.exists(POP_GROUPS_FILE):
    result = pd.read_csv(POP_GROUPS_FILE)
    summary = pd.read_csv(POP_GROUPS_FILE_SUMMARY)
else:
    # CONNECT TO THE DATABASE
    ag3 = malariagen_data.Ag3()

    # Fetch sample metadata from initial sample query
    df = ag3.sample_metadata(sample_query=sample_query_init)

    # Group based on location & year, 
    result = group_population(df, group_size=group_size)

    # Save the grouped individuals and summary of groups to a new CSV file
    result.to_csv("../data/1-METADATA_POPULATIONS.csv", index=False)
    summary = result['group_id'].value_counts().reset_index()
    summary.columns = ['group_id', 'count']
    print("summary of populations selected: ")
    print(summary)
    summary.to_csv("../data/1-METADATA_POPULATIONS_SUMMARY.csv")

# download all hdf5 files from malariagen before moving onto further processing.
for POP_NAME in summary['group_id']:
    # Check if .hdf5 file has already been downloaded from malariagen, download with get_data() if not
    file_path = "../data/" + POP_NAME + ".hdf5"
    if os.path.isfile(file_path):
        print(f"{file_path} exists, using existing files for analysis...")
    else:
       sample_ids = result[result['group_id'] == POP_NAME]['sample_id'].to_list()
       sample_query = f"sample_id in ({', '.join([repr(sample_id) for sample_id in sample_ids])})"
       ag3 = malariagen_data.Ag3()
       print(sample_query)
       get_data(sample_query=sample_query, POP_NAME=POP_NAME)

# loop back over each population once hdf5s exist, run lidrunner.sh to convert, run, extract, and plot individually
for POP_NAME in summary['group_id']:
    results_path = "../results/" + "LDhat_out_" + POP_NAME + ".RData"
    if os.path.isfile(results_path):
        print(f"LDhat results file already exists, skipping.." + results_path)
    else:
        # run lidrunner.sh on pop (hdf5-convert, run-LDhat, plot-LDhat)
        command = ["./lidrunner.sh", POP_NAME, str(WINDOW_WIDTH), str(SITES_PER_WINDOW), str(CHUNK), str(WINDOWS_PER_CHUNK), str(N), str(RANDSAMP)]
        subprocess.run(command)
        

# Call 5-plot_all_LDhat.R (current version just takes whatever it finds in ../results/LDhat_out_*)
#pop_names = summary['group_id'].tolist()
#pop_names_arg = ','.join(pop_names)
command = ["Rscript", "5-plot_all.R"]
subprocess.run(command)
