###############################################################################
# 13/09/2011 - Luke Peck
#
# SERPent Input File
#


##############################################################################
# This file contains all the information the user must input before running SERPent.
# Check all the python variables as most will be relevant!

# To run SERPent, type in terminal: parseltongue SERPent.py


# Data Information:


'''
AIPS_user_number = 4000            # The AIPS user number the data is on.
Name = '2007+404'                # The uvdata name of the catalogue.
Klass = '15MIN2'                  # The uvdata klass of the catalogue.
Disk = 1                          # The uvdata disk number of the catalogue (integer).
Seq = 1                          # The uvdata sequence number of the catalogue (int).
'''


AIPS_user_number = 4000            # The AIPS user number the data is on.
Name = '2007+404'                # The uvdata name of the catalogue.
Klass = 'SPLAT'                  # The uvdata klass of the catalogue.
Disk = 1                          # The uvdata disk number of the catalogue (integer).
Seq = 1                          # The uvdata sequence number of the catalogue (int).


'''
AIPS_user_number = 4001            # The AIPS user number the data is on.
Name = '1407+284_2'                # The uvdata name of the catalogue.
Klass = 'SPLAT'                  # The uvdata klass of the catalogue.
Disk = 1                          # The uvdata disk number of the catalogue (integer).
Seq = 1                           # The uvdata sequence number of the catalogue (int).
'''

# Note: if you're running SERPent on a multi source file, please SPLIT the sources
# first to the same disk and seq number as your MULTTB file but input the Klass as the
# multi source file (e.g. MULTTB). SERPent will try to find the SPLIT files...


# Log Information:

write2log = 'no'                 # 'yes' only if using the multi.py script file to write timing
                                 # details to a log file (only used for performance testing)



# Parallelization Information:

NCPU = 16            # Define here the number of CPUs you want to use.
                    # Parallelization is currently implemented on a baseline basis
                    # thus the maximum number of CPUs ultilized will be the number of
                    # baselines in the data. i.e. for e-MERLIN a maximum of 21 CPUs
                    # will be used by SERPent (21 baselines).



# Directory Information:

#path2folder = '/import/cobrasarchive1/data/COBRAS/LEGACY/L-band-May12/serpent_test/Lovell/'

path2folder = '/import/cobrasarchive1/data/dmf/'

#path2folder = '/import/cobrasarchive1/data/COBRAS/serpent_test/'

                    # Directory path to where the output text files (FG text files
                    # and intermediate pickle files) will be written to.


# Phase Cal Information

phasecal = '2007+4029'     # If one of the sources (multi files) or source is the
                           # phase cal, please write the name of the source as this
                           # variable. Else put 'no'.
                           # This information is used for the Lovell Stationary Scans passage of
                           # SERPent if the source is the phasecal designated here, the telescope
                           # is e-MERLIN and the baseline contains the Lovell telescope.



# Execute the Zero Level code to flag telescope dropouts and other system failures where the visibilities drop to around 0 Jy in the same scans where good data is present.


zero_level = 'no'        # If you want to execute this code set this variable to 'yes'
                          # else set it to 'no'



# Baseline Information for Flagging

which_baselines = 'all'
                              # Variable to define whether to flag all baselines
                              # or a select few.
                              # Options are: 'all' or 'choose'.


# If you chose which_baselines = 'choose', add the baselines to the below list:
# Note you'll obviously have to know the antenna numbers for any given baseline...
# If 'all' is chosen then the baseline list below is ignored

baselines = ['2-7']#'5-7','5-8','5-9','6-7','6-8','6-9'] #['2-7']
                  # chosen baselines for flagging.
                  # Format: ['5-7', '7-8'], i.e. for baselines 5-7 and 7-8 (strings!)
                  # order of baselines does not matter



# Flagger Options:

# This option is needed irrespective of whether flagging_options = 'choose' or 'default' below!
coadd_polarization_flags = 'no'    # Combine flags for all polarizations = 'yes', 'no' flags all
                                    # polarizations separately

coadd_zero_flags = 'no'            # If you choose not to combine flags for all polarizations 
                                   # coadd_zero_flags = 'no' will keep the zero dropout flags 
                                   # separate for each stokes, whereas coadd_zero_flags = 'yes'
                                   # will combine them

# To increase the flagging options available, some of the variable inputs can be
# modified here. There is already a 'default' setting which flags well, so its
# not essential to change anything here, and you should consult the SERPent Cookbook
# before playing with the below to understand what each variable does...



flagging_options = 'choose'
                            # variable to define whether to use the flagging options
                            # below or the default options in the SERPent.py file
                            # 'default' ignores whatever variables are set in this file
                            # Options are: 'choose' or 'default'



aggressiveness_first_run = 25       # How aggressive the first run of the flagger is
                                    # Note: a lower number is more aggressive

max_subset_first_run = 16           # Maximum subset for first run of the flagger
                                    # Note: should be a binary number (1,2,4,8,16...)

aggressiveness_second_run = 25      # How aggressive the second run of the flagger is


max_subset_second_run = 128         # Maximum subset for second run of the flagger


rho = 1.5                           # Difference in coarseness between each threshold
                                    # level

kickout_sigma_level = 3.0           # The kickout clause is tested during before each subset run through
                                    # and takes the form:  median + kickout_sigma_level * MAD
                                    # lower coefficient = deeper flagging/ more aggressive.

flag_coinc_chans = 0                # This will include flag_coinc_chans number of channels between 
                                    # flagged channels when writing the flags to the table akin to REFLG
                                    # i.e. if flag_coinc_chans = 3 and chans 2-3 and 7-11 are flagged, 
                                    # channels 4,5 & 6 will be as well. Default is 0.


# REFLG CPARM variables:

# The values defined below are used if flagging_options == 'choose'.

cparm_1 = 0.90       # 1: interval

cparm_2 = 1       # 2: flags channels if between cparm[2] flagged channels

cparm_3 = 1.00    # 3: fraction of channels flagged at specific time - flag all channels

cparm_4 = 1.00    # 4: fraction of times flagged at specific channel - flag all times

cparm_5 = 1.00    # 5: fraction of baselines flagged at specific channel and time - flag all bline

cparm_6 = 1.00    # 6: fraction of baselines flagged at specific antenna - flag all bline with ant

cparm_7 = 0       # 7: cparm[7] = 1, flag cross hand pols, cparm[7] = 2, flag all pols


# Outputs from SERPent:

# SERPent will attach some FG tables to your data file, 1 from the original run of
# SERPent and version 2 after REFLG (if your AIPS version has it). It will also output
# a number of text files associated with the different reduction passages in SERPent.
# These can be manipulated at the users leisure but are automatically included in the final FG table.
