# Lidrunner

**Date Created**: version 1.7 - 15 July 2024  
**Last Updated**: 9 August 2024

## Overview

`lidrunner` is a tool designed to automate the process of running LDhat on in sliding windows across genome data from the Malariagen API. It is well parallelized and should run comfortably on a personal laptop.

## Prerequisites

Before running `lidrunner`, ensure you have authenticated with the Malariagen API using `gcloud`. Follow these steps to authenticate your terminal:

```bash
gcloud auth application-default login
```

## Running the Program

The entire program is controlled by the script located at `code/runit.sh`. Before executing the script, you will need to modify certain parameters to suit your needs.

### Modifying Parameters

In the `runit.sh` script, adjust the following parameters on line 16:

```bash
WINDOW_WIDTH=sys.argv[1]       # Size of non-overlapping windows to run LDhat on
SITES_PER_WINDOW=sys.argv[2]   # Minimum threshold of biallelic SNPs per window
CHUNK=sys.argv[3]              # Size of non-overlapping chunk to sample windows from
WINDOWS_PER_CHUNK=sys.argv[4]  # Maximum number of windows to sample in each chunk
N=sys.argv[5]                  # Population size - *currently hard coded to 50 - leave as 50!*
RANDSAMP=sys.argv[6]           # Boolean - randomly sample from windows which meet the threshold?
group_size=int(sys.argv[7])    # Leave as 50
sample_query_init=sys.argv[8]  # Pandas sample query to create populations from
```

### Execution

Once you've set the desired parameters, run the script by executing the following command:

```bash
bash code/runit.sh
```

## Contributions

Contributions to this repository are welcome. Please submit issues or pull requests as needed.

## License
.
```

