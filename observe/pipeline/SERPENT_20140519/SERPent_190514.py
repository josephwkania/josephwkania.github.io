#!/usr/bin/python

###############################################################################
# 13/09/2011 - Luke Peck
#
# Scripted E-merlin Rfi-mitigation PypelinE for iNTerferometry - SERPent
# Version 19/11/2013


##############################################################################
# AIPS and other Python Modules
##############################################################################


from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
# First three imports are starter modules but more might/ will be needed

from AIPSTV import AIPSTV
import Wizardry.AIPSData
from Wizardry.AIPSData import AIPSUVData
# Wizardry is a useful module for clever manipulation of AIPS data

import LocalProxy
from xmlrpclib import ServerProxy

import copy, optparse, os, sys
import re, string
import os.path
#import numarray.ma
#import numarray.ieeespecial
#import numarray.nd_image
import inspect
import warnings
warnings.defaultaction = "always"

import numpy
numpy.set_printoptions(threshold='nan')
import platform
import cPickle as pickle

import math, time, datetime
from time import gmtime, strftime, localtime
ti = time.time()    # To time the script


###############################################################################
# Global Functions
################################################################################


# Function to convert seconds to hh:mm:ss.ss format, returns a string
def time2hms(seconds):
    h=int(seconds/3600)
    m=int(seconds % 3600)/60
    s=seconds-(h*3600)-(m*60)
    h=`h`
    m=`m`
    s="%4.2f" % s
    hms=h.zfill(2)+":"+m.zfill(2)+":"+s.zfill(4)
    return hms


# Function to convert time to d/hh:mm:ss.s format which is what is required for FG tables
def time2dhms(seconds):
    conv = seconds * 86400
    d = int(conv/ 86400)
    h=int(conv % 86400)/ 3600
    m=int(conv % 3600)/60
    s=conv-(d*86400)-(h*3600)-(m*60)
    d=`d`
    h=`h`
    m=`m`
    s="%3.1f" % s
    dhms= d.zfill(1) + "/" + h.zfill(2) + ":" + m.zfill(2) + ":" + s.zfill(4)
    return dhms


# Definition to give REASON string to find the date and give REASON in the FG table in the format like IBLED does i.e. name of program + date + time...
def reason():
    reason_str = 'SERPent ' + strftime("%d-%b-%y") + ' ' + strftime("%H:%M", localtime())
    return reason_str.upper()


# Determines the array size
def array_size(data_mem):
    # Returns the size of the array with the appropriate suffix
    if len(str(data_mem)) <= 6:
        return "%.3f KB" % (data_mem / 10.0**3)
    elif len(str(data_mem)) <= 9:
        return "%.3f MB" % (data_mem / 10.0**6)
    elif len(str(data_mem)) <= 12:
        return "%.3f GB" % (data_mem / 10.0**9)
    elif len(str(data_mem)) <= 15:
        return "%.3f TB" % (data_mem / 10.0**12)


# Frequency prefix calculator
def freq_prefix(num):
    if num/(10**12) > 1:
        return "%.5f THz" % (num / 10.0**12)
    elif num/(10**9) > 1:
        return "%.5f GHz" % (num / 10.0**9)
    elif num/(10**6) > 1:
        return "%.5f MHz" % (num / 10.0**6)
    elif num/(10**3) > 1:
        return "%.5f KHz" % (num / 10.0**3)


# Frequency band calculator. Only works for freq between 1.2GHz - 26.5GHz but easy to add on more...
def freq_band(freq):
    if freq > 1.2E+9 and freq < 2.0E+9:
        return "L Band"
    elif freq > 2.0E+9 and freq < 4.0E+9:
        return "S Band"
    elif freq > 4.0E+9 and freq < 8.0E+9:
        return "C Band"
    elif freq > 8.0E+9 and freq < 12.0E+9:
        return "X Band"
    elif freq > 12.0E+9 and freq < 18.0E+9:
        return "Ku Band"
    elif freq > 18.0E+9 and freq < 26.5E+9:
        return "K Band"
    else:
        return "Out of range for calculator"


# Definition to find computer memory stats (works only for csh and bash at the moment)...
def computer_memory():
    memfile = open(path2folder + 'mem_stats.txt', 'wr')
    memfile.flush()
    file = path2folder + 'mem_stats.txt'
    #os.system('free -k >' + file)    # need to test this!
    os.system('free -k > mem_stats.txt')    # need to put in file path
    memfile = open(path2folder + 'mem_stats.txt', 'r')
    stats = []
    for line in memfile:
        string = ''
        for char in xrange(len(line)):
            if line[char].isdigit():
                string += line[char]
                if line[char+1] == ' ' or line[char+1] == '\n':
                    stats.append(int(string))
                    string = ''
    global mem_stats
    mem_stats = {}
    mem_stats['mem_total'] = stats[0]
    mem_stats['mem_used'] = stats[1]
    mem_stats['mem_free'] = stats[2]
    mem_stats['mem_shared'] = stats[3]
    mem_stats['mem_buffers'] = stats[4]
    mem_stats['mem_cached'] = stats[5]
    mem_stats['buf/cach_used'] = stats[6]
    mem_stats['buf/cach_free'] = stats[7]
    mem_stats['swap_total'] = stats[8]
    mem_stats['swap_used'] = stats[9]
    mem_stats['swap_free'] = stats[10]
    memfile.close()
    return mem_stats


# Task to find the source names and numbers if an SU table doesn't exist.
def indxr():
    indxr = AIPSTask('INDXR')
    indxr.inname = uvdata.name
    indxr.inclass = uvdata.klass
    indxr.inseq = uvdata.seq
    indxr.indisk = uvdata.disk
    #indxr.cparm = # sets length of scan, begin and end times of scan
    #indxr.bparm = # opacity and gain curve control
    indxr.inp()
    indxr.go()


###############################################################################
# SERPent Script
###############################################################################

version = "25/04/14"


print '\n Started running SERPent version: %s' % version, 'on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),'\n'


#print os.getcwd()
#path1=os.getcwd()
#if os.path.exists(path1+'/SERPent_input.py'):
    ## print 'hellow!'
    #execfile("SERPent_input.py")

#Execute the SERPent_input.py file with all the observation information and flagging parameters...
try:
    execfile("SERPent_input.py")
    #execfile("SERPent_input_multi.py")
    print AIPS_user_number
    
except:
    print "\n Could not find SERPent_input.py!"
    print " Make sure input file is in the same folder as this script (SERPent.py)"
    print " Aborting!"
    sys.exit(0)



# Test to see whether all variable have been set:
try:
    AIPS_user_number
    Name
    Klass
    Disk
    Seq
    NCPU
    path2folder
    if os.path.exists(path2folder) == False:
        print " Folder for outputs does not exist! Please check inputs and try again."
        print " Aborting SERPent!\n"
        sys.exit(0)
    phasecal
    zero_level
    which_baselines
    if which_baselines == 'choose':
        baselines
    coadd_polarization_flags
    if coadd_polarization_flags == 'no':
        coadd_zero_flags
        # coadd_lovell_flags
    flagging_options
    if flagging_options == 'choose':
        aggressiveness_first_run
        max_subset_first_run
        aggressiveness_second_run
        max_subset_second_run
        rho
        kickout_sigma_level
        flag_coinc_chans
except NameError:
    print " Please check you\'ve set all the variables in the input file."
    print " ABORTING SERPent!\n"
    sys.exit(0)



# Read parameters from input file to local variable
parameters = [aggressiveness_first_run, max_subset_first_run, aggressiveness_second_run, max_subset_second_run, rho, kickout_sigma_level]


AIPS.userno = AIPS_user_number
experiment = Name
experiment_klass = Klass
experiment_disk = Disk
experiment_seq = Seq


# print AIPS.userno
# print experiment, experiment_klass, experiment_disk, experiment_seq


cat = AIPSCat(Disk)


uvdata = AIPSUVData(Name, Klass, Disk, Seq)



# Destroy the fg table for trial runs (and mem_stats file)...
if os.path.exists(path2folder+experiment+'.fg') == True:
    os.remove(path2folder+experiment+'.fg')
if os.path.exists(path2folder+experiment+'_r.fg') == True:
    os.remove(path2folder+experiment+'_r.fg')
if os.path.exists(path2folder+'dummy.fg') == True:
    os.remove(path2folder+'dummy.fg')
if os.path.exists(path2folder+'mem_stats.txt') == True:
    os.remove(path2folder+'mem_stats.txt')
if os.path.exists(path2folder+experiment+'_lovell_dropouts.fg') == True:
    os.remove(path2folder+experiment+'_lovell_dropouts.fg')
if os.path.exists(path2folder+experiment+'_lovell_dummy.fg') == True:
    os.remove(path2folder+experiment+'_lovell_dummy.fg')
if os.path.exists(path2folder+experiment+'_combined.fg') == True:
    os.remove(path2folder+experiment+'_combined.fg')
if os.path.exists(path2folder+experiment+'_zeros_dropouts.fg') == True:
    os.remove(path2folder+experiment+'_zeros_dropouts.fg')
if os.path.exists(path2folder+experiment+'_zeros_dummy.fg') == True:
    os.remove(path2folder+experiment+'_zeros_dummy.fg')



# again destroy FG extensions for trial runs...
uvdata.zap_table('FG', -1)


source_list = {}

# First try and see whether there are any SU tables, and get source numbers and names from that.
try:
    nsource = len(uvdata.sources)
    print " Number of sources (nsource): %i" % nsource
    for s in xrange(len(uvdata.sources)):
        source_list[s+1] = str(uvdata.sources[s])
    print " Names of sources (with dictionary index):"
    for s in source_list:
        if not s == nsource:
            print " %i: %s," % (s, source_list[s]),
        else:
            print " %i: %s" % (s, source_list[s])
    #print " Sources:", source_list
except:
    print "\n No SU table found... "
    print " Assuming single source..."
    nsource = 1
    source_list[1] = Name.upper()
    print " %i: %s" % (1, source_list[1])


# Number of rows in the FG table. VERY IMPORTANT! DO NOT DELETE! AND MUST BE BEFORE THE SOURCE FOR LOOP.
no_rows = 0
final_num_vis = 0
global fg_row
fg_row = 1
pickle_list = {}
total_tscans = 1
nif = 1
nchan = 1
npol = 1
num_flagged_vis = 0

# Loop for all sources in the source list:
for source in source_list:
    srcname = source_list[source]
    if nsource > 1:
        try:
            uvdata = AIPSUVData(source_list[source], "SPLIT ", Disk, Seq)
        except:
            try:
                print "\n IGNORE ERROR MESSAGES. TRYING WITH DIFFERENT KLASS NAME..."
                uvdata = AIPSUVData(source_list[source], Klass, Disk, Seq)
            except:
                print "\n If you're running this on a multi source, SPLIT the sources first on the same disk and sequence number as the multi-source, with the klass either as \'SPLIT\' or the same klass as the multi-source file and then run SERPent on the multi-source file again..."
                print "\n Aborting!"
                sys.exit(0)
    else:
        #try:
            #uvdata = AIPSUVData(source_list[source], Klass, Disk, Seq)
        #except:
        print "\n IGNORE ERROR MESSAGES. TRYING WITH DIFFERENT UVDATA NAME..."
        uvdata = AIPSUVData(Name, Klass, Disk, Seq)

    # Information about the observations:
    name = uvdata.name
    print "\n RUNNING SCRIPT ON SOURCE: %s" % name



    # Name of Telescope:
    telescope = uvdata.header['telescop']
    print " Name of Telescope: %s" % telescope


    #uvdata = AIPSUVData('MULTI', 'UVDATA', 1, 1)    # sets up local variable
    nvis = len(uvdata)    # prints the number of visibilities
    print " Number of visibilities (nvis): %i" % nvis


    # Find the Frequency of the observation
    frequency = uvdata.header['crval'][uvdata.header['ctype'].index('FREQ')]
    print " Frequency of the Observation: %s" % freq_prefix(frequency)
    print " Frequency band of the Observation: %s" % freq_band(frequency)


    # Figure out number of Stokes and polarizations
    npol = len(uvdata.stokes)
    print " Number of Polarizations (npol): %i" % npol
    polnames = {}
    for x in xrange(npol):
        polnames[str(uvdata.stokes[x])] = x
    # Below is just to print out the stokes in a nice way :)
    string = ""
    for x in xrange(npol):
        if x != npol-1:
            string += (str(uvdata.stokes[x])+ ", ")
        else:
            string += str(uvdata.stokes[x])
    print " Polarizations: %s" % string





    # figure out number of IFs
    if 'IF' in uvdata.header['ctype']:
        nif = uvdata.header['naxis'][uvdata.header['ctype'].index('IF')]
        print " Number of IFs (nif): %i" % nif
    else:
        print " Keyword IF not found in header, assuming number of IFs is 1"
        nif = 1


    # figure out number of channels
    nchan = uvdata.header['naxis'][uvdata.header['ctype'].index('FREQ')]
    print " Number of Channels (nchan): %i" % nchan


    # Figure out number of baselines
    nant = len(uvdata.antennas)
    nbase = ((nant*(nant-1))/2)
    nacbase = nbase + nant
    print " Number of Possible Baselines (nbase): %i" % nacbase


    # Dictionary containing number of time scans for each baseline.
    if which_baselines == 'all':
        # print "\n ALL %s baselines have been selected for flagging..." % nbase

    # List of baselines. Note that this uses the actual visibility data to find the antenna numbers and doesn't assume they start at 1 and increase linearly by the same amount. This has been written because the VLA and early commissioning e-MERLIN data have a unique (i.e. stupid) numbering system.
    # Note that this is only in the 'all' baseline choice as it then needs to find all the baselines but doesn't if they've been stated before...
        baselines = []
        for vis in uvdata:
            bline = "%i-%i" % (vis.baseline[0], vis.baseline[1])
            if bline not in baselines:
                baselines.append(str(bline))
            elif len(baselines) == nacbase:
                break
        # print baselines
        print "\n ALL %s baselines have been selected for flagging..." % len(baselines)


        if len(baselines) > nbase:
            print 'Autocorrelations detected: these will be flagged'
#        print baselines

        ntimescans = {}
        ntime = 0
        for bline in baselines:
            ntimescans[bline] = 0
        for vis in uvdata:
            bline = "%i-%i" % (vis.baseline[0], vis.baseline[1])
            ntimescans[bline] += 1
            ntime += 1
        print "List of ALL baselines with corresponding number of time scans:", ntimescans

 

    # Testing for specified baselines.
    if which_baselines == 'choose':
        print "\n SPECIFIC baselines selected for flagging..."

        ntimescans = {}
        ntime = 0
        for bline in baselines:
            ntimescans[bline] = 0
        for vis in uvdata:
            bline = "%i-%i" % (vis.baseline[0], vis.baseline[1])
            if bline in baselines:
                ntimescans[bline] += 1
                ntime += 1
        print " List of CHOSEN baselines with corresponding number of time scans:", ntimescans


    # List for baselines to remove from SERPent job list and flag at end of script
    baseline_removal = {}

    # Removes Lovell-Mk2 baselines as this baselines is pretty much useless.
    antable = uvdata.table('AN', 0)
    lovellmk2 = 0
    # print baselines
    if telescope in ('e-MERLIN', 'E-MERLIN', 'eMERLIN', 'EMERLIN'):
        for ant in xrange(len(uvdata.antennas)):
            if uvdata.antennas[ant].upper() in ('LO', 'LOV', 'LOVE', 'LOVEL', 'LOVELL'):
                lo = antable[ant]['nosta']
            if uvdata.antennas[ant].upper() in ('MK2', 'M2', 'MARK2', 'MK'):
                mk = antable[ant]['nosta']
                #print 'MK is', mk
        try:
            lo
            mk
        except:
            lo = 0
            mk = 0
        for bline in baselines:
            if str(lo) in bline and str(mk) in bline:
                print "\n Removing Lovell-Mk2 baseline from baseline list."
                lovellmk2 = bline
                baseline_removal[bline] = 1
                baselines.remove(bline)
                break


    # Remove auto-correlations...
    antable = uvdata.table('AN', 0)
    for ant in xrange(len(uvdata.antennas)):
        ant_num = antable[ant]['nosta']
        auto_bline = str(ant_num) + '-' + str(ant_num)
        #print ant_num, auto_bline
        for bline in baselines:
            if auto_bline == bline:
                print " Removing Auto-correlation baseline %s from baseline list." % (auto_bline)
                baseline_removal[bline] = 1
                baselines.remove(bline)
                break


    # Checks for Linux-based OS then runs memory definition
    try:
        computer_memory()
    except:
        print "\n Sorry, computer_memory() definition does not work on your system!"


    if 'mem_total' in mem_stats:
        print "System Memory Information:"
        print "Total Memory  :    %s" % array_size(mem_stats['mem_total']*1000)
        print "Used Memory   :    %s" % array_size(mem_stats['mem_used']*1000)
        print "Free Memory   :    %s" % array_size(mem_stats['mem_free']*1000)
        print "Total Swap    :    %s" % array_size(mem_stats['swap_total']*1000)
        print "Used Swap     :    %s" % array_size(mem_stats['swap_used']*1000)
        print "Free Swap     :    %s" % array_size(mem_stats['swap_free']*1000)


    # The predicted array sizes, these predict the arrays when flagging all baselines simultaneously
    pred_array = nif*nchan*ntime*npol*3*8
    print "Total predicted memory usage: %s" % array_size(int(pred_array))

    decision = 1



###############################################################################
#   Flagging algorithm split into baselines and IFs (for each baseline work through the IFs individually then move on to the next baseline)


# Begin the parallelization block:

########################################################################################
# CREDIT: to John Conway, Ivan Marti-Vidal, Miguel Perez Torres and Eskil Varenius, whose code gave the foundation from which this parallelization is based on.
########################################################################################


    nbline = len(baselines)    # need to change
    total_jobs = nbline*nif


    total_tscans = 0
    #global pickle_list
    pickle_list = {}
    lovell_bline = []
    global lovell_row
    lovell_row = 1
    zeros_bline = []
    global zeros_row
    zeros_row = 1

    # Find the min number of CPUs needed either via NCPU or number of baselines
    nfork = min(total_jobs,NCPU)
    dscan = int(total_jobs/nfork)
    mybline = [range(dscan*i,dscan*(i+1)) for i in range(nfork)]
    # print "mybline", mybline

    # Find remaining jobs after initial designation
    remain = total_jobs-mybline[-1][-1]


    # Append remaining jobs onto cores equally. Note that the designation might not be neatly distributed but it is equally distributed.
    for core,i in enumerate(range(remain-1,0,-1)):
        mybline[core] += [total_jobs-i]
        # print mybline
        # print mybline[core]
    # print "mybline", mybline


    # Giving the job numbers from 0 to total number of IFs over all baselines to each baseline. Note that the job numbers do not refer to the IF number
    job = 0
    blineif = {}
    for b in baselines:
        blineif[b] = []
        for j in xrange(nif):
            blineif[b].append(job)
            job += 1

    # print "blineif",  blineif



    # Now finding the correct IF number for the jobs given for each core
    blist = {}
    for c in xrange(len(mybline)):
        for base in xrange(len(mybline[c])):
            blist[c] = {}
            for b in blineif:
                for job in blineif[b]:
                    if job in mybline[c]:
                        if job >= nif:
                            while job >= nif:
                                job -= (nif)
                        try:
                            blist[c][b].append(job)
                        except:
                            blist[c][b] = []
                            blist[c][b].append(job)
            break
        print "CPU", c, "will flag baselines and IFs: \n" , blist[c]

 
    pid = []
    for nf in range(nfork-1):
        if 0 not in pid:
            pid.append(os.fork())

    if 0 not in pid:
        mycore = 0
    else:
        mycore = len(pid)

    mynbline = len(mybline[mycore])



    # Code to append visibilities
    relsc = -1
    for bline in blist[mycore]:


        if decision == 1:


            # where does this list need to be? before or after the for bline in blist[mycore]?
            lovell_bline = []


            for j in blist[mycore][bline]:
                #lovell_bline = []
                relsc += 1

                lovell_row = 1


                print "Source:", name, ", Baseline:", bline, ", IF:", j+1, " on CPU", mycore, ", Job #:", relsc+1, " of", mynbline


                times_array = numpy.zeros([ntimescans[bline]], dtype='|S12')
                # print ntimescans[bline], bline
                # print len(times_array)
                appending_time = 0
                if 'RR' in polnames or 'XX' in polnames:
                    amp_time_RR = numpy.zeros([ntimescans[bline], nchan])
                    amp_freq_RR = numpy.zeros([ntimescans[bline], nchan])
                    flag_RR = numpy.zeros([ntimescans[bline], nchan])
                    appending_count_RR = 0
                if 'LL' in polnames or 'YY' in polnames:
                    amp_time_LL = numpy.zeros([ntimescans[bline], nchan])
                    amp_freq_LL = numpy.zeros([ntimescans[bline], nchan])
                    flag_LL = numpy.zeros([ntimescans[bline], nchan])
                    appending_count_LL = 0
                if 'RL' in polnames or 'XY' in polnames:
                    amp_time_RL = numpy.zeros([ntimescans[bline], nchan])
                    amp_freq_RL = numpy.zeros([ntimescans[bline], nchan])
                    flag_RL = numpy.zeros([ntimescans[bline], nchan])
                    appending_count_RL = 0
                if 'LR' in polnames or 'YX' in polnames:
                    amp_time_LR = numpy.zeros([ntimescans[bline], nchan])
                    amp_freq_LR = numpy.zeros([ntimescans[bline], nchan])
                    flag_LR = numpy.zeros([ntimescans[bline], nchan])
                    appending_count_LR = 0
                # for pol in polnames:
                    # print 'Arrays for %s polarization created!' % pol


                ap = time.time()
                memory_array = 0
                print " CPU", mycore, " Appending visibilities..."
                for vis in uvdata:
                    index = "%i-%i" % (vis.baseline[0], vis.baseline[1])
                    if index == bline:
                        count_time = appending_time
                        times_array[count_time] = float(vis.time)
                        appending_time += 1
                        if 'RR' in polnames or 'XX' in polnames:
                            if 'RR' in polnames:
                                pol = 'RR'
                            if 'XX' in polnames:
                                pol = 'XX'
                            amp_time_RR[appending_count_RR , 0:(nchan)] = numpy.where(vis.visibility[j][:,polnames[pol]][:,2] > 0.0, (numpy.sqrt(vis.visibility[j][:,polnames[pol]][:,0]**2 + vis.visibility[j][:,polnames[pol]][:,1]**2)), float('NaN'))
                            amp_freq_RR[appending_count_RR , 0:(nchan)] = numpy.where(vis.visibility[j][:,polnames[pol]][:,2] > 0.0, (numpy.sqrt(vis.visibility[j][:,polnames[pol]][:,0]**2 + vis.visibility[j][:,polnames[pol]][:,1]**2)), float('NaN'))
                            appending_count_RR += 1
                        if 'LL' in polnames or 'YY' in polnames:
                            if 'LL' in polnames:
                                pol = 'LL'
                            if 'YY' in polnames:
                                pol = 'YY'
                            amp_time_LL[appending_count_LL , 0:(nchan)] = numpy.where(vis.visibility[j][:,polnames[pol]][:,2] > 0.0, (numpy.sqrt(vis.visibility[j][:,polnames[pol]][:,0]**2 + vis.visibility[j][:,polnames[pol]][:,1]**2)), float('NaN'))
                            amp_freq_LL[appending_count_LL , 0:(nchan)] = numpy.where(vis.visibility[j][:,polnames[pol]][:,2] > 0.0, (numpy.sqrt(vis.visibility[j][:,polnames[pol]][:,0]**2 + vis.visibility[j][:,polnames[pol]][:,1]**2)), float('NaN'))
                            appending_count_LL += 1
                        if 'RL' in polnames or 'XY' in polnames:
                            if 'RL' in polnames:
                                pol = 'RL'
                            if 'XY' in polnames:
                                pol = 'XY'
                            amp_time_RL[appending_count_RL , 0:(nchan)] = numpy.where(vis.visibility[j][:,polnames[pol]][:,2] > 0.0, (numpy.sqrt(vis.visibility[j][:,polnames[pol]][:,0]**2 + vis.visibility[j][:,polnames[pol]][:,1]**2)), float('NaN'))
                            amp_freq_RL[appending_count_RL , 0:(nchan)] = numpy.where(vis.visibility[j][:,polnames[pol]][:,2] > 0.0, (numpy.sqrt(vis.visibility[j][:,polnames[pol]][:,0]**2 + vis.visibility[j][:,polnames[pol]][:,1]**2)), float('NaN'))
                            appending_count_RL += 1
                        if 'LR' in polnames or 'YX' in polnames:
                            if 'LR' in polnames:
                                pol = 'LR'
                            if 'YX' in polnames:
                                pol = 'YX'
                            amp_time_LR[appending_count_LR , 0:(nchan)] = numpy.where(vis.visibility[j][:,polnames[pol]][:,2] > 0.0, (numpy.sqrt(vis.visibility[j][:,polnames[pol]][:,0]**2 + vis.visibility[j][:,polnames[pol]][:,1]**2)), float('NaN'))
                            amp_freq_LR[appending_count_LR , 0:(nchan)] = numpy.where(vis.visibility[j][:,polnames[pol]][:,2] > 0.0, (numpy.sqrt(vis.visibility[j][:,polnames[pol]][:,0]**2 + vis.visibility[j][:,polnames[pol]][:,1]**2)), float('NaN'))
                            appending_count_LR += 1
                print "CPU", mycore, "Append time (hh:mm:ss):", time2hms(time.time()-ap)
                if 'RR' in polnames or 'XX' in polnames:
                    memory_array += 3*amp_time_RR.nbytes
                if 'LL' in polnames or 'YY' in polnames:
                    memory_array += 3*amp_time_LL.nbytes
                if 'RL' in polnames or 'XY' in polnames:
                    memory_array += 3*amp_time_RL.nbytes
                if 'LR' in polnames or 'YX' in polnames:
                    memory_array += 3*amp_time_LR.nbytes
                # print "Total memory size of arrays: %s" % array_size(memory_array)
                

                # Test to see if entire IF is Nan's from the correlator for each stoke
                print "\n CPU", mycore, 'Checking for NaN values in array...'
                percentage = 1.0    # 0 to 1.0
                if 'RR' in polnames or 'XX' in polnames:
                    nan_count = 0
                    nan_array_RR = 0
                    for col in xrange(amp_time_RR.shape[1]):
                        for row in xrange(amp_time_RR.shape[0]):
                            if str(amp_time_RR[row, col]) == "nan":
                                nan_count += 1
                    if nan_count == (percentage*amp_time_RR.shape[0]*amp_time_RR.shape[1]):
                        nan_array_RR = 1
                if 'LL' in polnames or 'YY' in polnames:
                    nan_count = 0
                    nan_array_LL = 0
                    for col in xrange(amp_time_LL.shape[1]):
                        for row in xrange(amp_time_LL.shape[0]):
                            if str(amp_time_LL[row, col]) == "nan":
                                nan_count += 1
                    if nan_count == (percentage*amp_time_LL.shape[0]*amp_time_LL.shape[1]):
                        nan_array_LL = 1
                if 'RL' in polnames or 'XY' in polnames:
                    nan_count = 0
                    nan_array_RL = 0
                    for col in xrange(amp_time_RL.shape[1]):
                        for row in xrange(amp_time_RL.shape[0]):
                            if str(amp_time_RL[row, col]) == "nan":
                                nan_count += 1
                    if nan_count == (percentage*amp_time_RL.shape[0]*amp_time_RL.shape[1]):
                        nan_array_RL = 1
                if 'LR' in polnames or 'YX' in polnames:
                    nan_count = 0
                    nan_array_LR = 0
                    for col in xrange(amp_time_LR.shape[1]):
                        for row in xrange(amp_time_LR.shape[0]):
                            if str(amp_time_LR[row, col]) == "nan":
                                nan_count += 1
                    if nan_count == (percentage*amp_time_LR.shape[0]*amp_time_LR.shape[1]):
                        nan_array_LR = 1

                #print "Baseline: %s, IF: %i, RR: %i, LL: %i, RL: %i, LR: %i." % (bline, j+1, nan_array_RR, nan_array_LL, nan_array_RL, nan_array_LR)


                ## removed mad_1d def from here ##



                # LOVELL STATIONARY SCAN CODE. ONLY FOR E-MERLIN.


                # Testing to see whether telescope is e-MERLIN and finding Lovell antenna number:
                lovell_num = 0
                covbase = [] 
                if telescope in ('e-MERLIN', 'E-MERLIN', 'eMERLIN', 'EMERLIN'):
                    print " Array is %s" % telescope
                    for ant in xrange(len(uvdata.antennas)):
                        if uvdata.antennas[ant].upper() in ('LO', 'LOV', 'LOVE', 'LOVEL', 'LOVELL'):
                            lovell_num = ant+1
                            print "Lovell number:", lovell_num
                            break
                        
                # print lovell_num, bline, lovell_bline, phasecal, srcname


                # If current baseline contains Lovell and hasn't been checked with a different IF then proceed.
                if str(lovell_num) in bline and bline not in lovell_bline and phasecal == srcname:
                    print "\n CPU ", mycore, "Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline

                    def find_source_index(source_list):
                        for s in source_list:
                            if source_list[s] == uvdata.name:
                                source = s
                                break
                        return source

                    ## removed mad_1d def from here ##

                    # Select a single channel to perform test on. Choosing channel in the middle of the IF for best quality amps...
                    test_chan = int(nchan/2)


                    con = numpy.zeros([len(times_array), 2], dtype='|S12')
                    test = len(times_array)
                    
                    for t in xrange(len(times_array)):
                        con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
                        
                    int_time = numpy.zeros([len(times_array)-1])
                    for t in xrange(1, len(times_array)):
                        int_time[t-1] = float(times_array[t]) - float(times_array[t-1])
                        
                    #print 'Lovell polarisation is', pol, j+1 
                    step = float(numpy.median(int_time))
                    print " Correlator/ Averaged integration time is:", time2dhms(step), 'seconds'

                    checkedpol = 0
                    dropouts = []


                    if 'LL' in polnames or 'YY' in polnames:
                        one = amp_time_LL[:,test_chan]*10**6
                        print 'Running Lovell check on LL or YY' 
                        checkedpol = 1
                    elif 'RR' in polnames or 'XX' in polnames:
                        one = amp_time_RR[:,test_chan]*10**6
                        print 'Running Lovell check on RR or XX' 
                        checkedpol = 2
                    elif 'RL' in polnames or 'XY' in polnames:
                        one = amp_time_RL[:,test_chan]*10**6
                        print 'Running Lovell check on RL or XY' 
                        checkedpol = 3
                    elif 'LR' in polnames or 'YX' in polnames:
                        one = amp_time_LR[:,test_chan]*10**6
                        print 'Running Lovell check on LR or YX' 
                        checkedpol = 4


                    continuenow = 1
                    if continuenow == 1:
                        # print 'continue'
                        # print 'checking pol', checkedpol

                        levels = {'high': [0,0], 'low': [10**10,0]}

                        # Firstly cycle through baseline and find mean, std of scans and set the lowest and highest levels.

                        count_drop = 0
                        begin_step = 1
                        end_step = 0
                        for i in xrange(1,len(one)):
                            # print i, one[i], con[i][0], con[i][1], con[i-1][0]
                            # Some test level of 5 * integration time...
                            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                                # print "Scan has changed level"
                                end_step = i
                                newo = one[begin_step-1:end_step]
                                avg = numpy.mean(newo[~numpy.isnan(newo)])
                                s = numpy.std(newo[~numpy.isnan(newo)])
                                if avg > levels['high'][0]:
                                    levels['high'][:] = [avg, s]
                                if avg < levels['low'][0]:
                                    levels['low'][:] = [avg, s]
                                begin_step = end_step
                                count_drop += 1
                            if i == len(one)-1:
                                # print "Last scan"
                                end_step = i
                                newo = one[begin_step-1:end_step+1]
                                avg = numpy.mean(newo[~numpy.isnan(newo)])
                                s = numpy.std(newo[~numpy.isnan(newo)])
                                if avg > levels['high'][0]:
                                    levels['high'][:] = [avg, s]
                                if avg < levels['low'][0]:
                                    levels['low'][:] = [avg, s]
                                begin_step = end_step


                        # now flag all scan levels within lowest avg + lowest s
                        count_drop = 0
                        begin_step = 1
                        end_step = 0
                        for i in xrange(1,len(one)):
                            # print i, one[i], con[i][0], con[i][1], one[0],one[1]
                            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                                # print float(con[i][0]) - float(con[i-1][0])
                                end_step = i
                                newo = one[begin_step-1:end_step]
                                avg = numpy.mean(newo[~numpy.isnan(newo)])
                                s = numpy.std(newo[~numpy.isnan(newo)])
                                # print "Average: %f, std: %f" % (avg, s)
                                if avg < levels['low'][0] + levels['low'][1] and avg > levels['low'][0] - levels['low'][1]:
                                    one[begin_step-1:end_step+1] = 0.0
                                    count_drop += 1
                                begin_step = end_step
                                # count_drop += 1
                            if i == len(one)-1:
                                end_step = i
                                newo = one[begin_step-1:end_step+1]
                                avg = numpy.mean(newo[~numpy.isnan(newo)])
                                s = numpy.std(newo[~numpy.isnan(newo)])
                                # print "Average: %f, std: %f" % (avg, s)
                                if avg < levels['low'][0] + levels['low'][1] and avg > levels['low'][0] - levels['low'][1]:
                                    one[begin_step-1:end_step+1] = 0.0
                                begin_step = end_step
                                # print end_step, one
                                
                                # Now collect indexes for amp == 0.0


                        if count_drop > 0:
                            dropouts = []
                            
                            # print 'len one is', len(one)
                            for i in xrange(len(one)):
                                #print i+1, one[i], con[i][0], con[i][1]
                                if one[i] == 0.0:
                                    dropouts.append(i)
                                # if i == len(one)-1:
                                  #   if one[i] == 0.0:
                                    #    dropouts.append(i+1)
                            #print dropouts



                        if count_drop > 0:
                            # for pol in polnames:
                            plname = str(name) + "__" + str(bline) + "_IF" + str(j+1) + "__lovell_dummy"
                            # infolist = numpy.array([[source], [con], [dropouts]])
                            pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                            pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                            pickle.dump(dropouts, open(path2folder + str(plname) + ".dropouts", "wb"))
                            covbase.append(j)
                            
                            if len(covbase) == len(blist[mycore][bline]):
                                # print 'hello!', blist[mycore]
                                lovell_bline.append(bline)
                        


                            times_array = numpy.delete(times_array, dropouts)
                            if ('RR' in polnames or 'XX' in polnames) and count_drop > 0:
                                amp_time_RR = numpy.delete(amp_time_RR, dropouts, axis=0)
                                amp_freq_RR = numpy.delete(amp_freq_RR, dropouts, axis=0)
                                flag_RR = numpy.delete(flag_RR, dropouts, axis=0)
                            if ('LL' in polnames or 'YY' in polnames) and count_drop > 0:
                                amp_time_LL = numpy.delete(amp_time_LL, dropouts, axis=0)
                                amp_freq_LL = numpy.delete(amp_freq_LL, dropouts, axis=0)
                                flag_LL = numpy.delete(flag_LL, dropouts, axis=0)
                            if ('RL' in polnames or 'XY' in polnames) and count_drop > 0:
                                amp_time_RL = numpy.delete(amp_time_RL, dropouts, axis=0)
                                amp_freq_RL = numpy.delete(amp_freq_RL, dropouts, axis=0)
                                flag_RL = numpy.delete(flag_RL, dropouts, axis=0)
                            if ('LR' in polnames or 'YX' in polnames) and count_drop > 0:
                                amp_time_LR = numpy.delete(amp_time_LR, dropouts, axis=0)
                                amp_freq_LR = numpy.delete(amp_freq_LR, dropouts, axis=0)
                                flag_LR = numpy.delete(flag_LR, dropouts, axis=0)

                                
                            

                # If same baseline as before, use same dropout indices for all pols.
                #elif str(lovell_num) in bline and bline in lovell_bline and phasecal == srcname:
                # elif bline in lovell_bline and phasecal == name:
                 #   print "\n CPU ", mycore, " Applying previously found Lovell dropouts on data for flagging."
                  #  if count_drop > 0:
                   #     times_array = numpy.delete(times_array, dropouts)
                    #    if ('RR' in polnames or 'XX' in polnames) and count_drop > 0:
                     #       amp_time_RR = numpy.delete(amp_time_RR, dropouts, axis=0)
                      #      amp_freq_RR = numpy.delete(amp_freq_RR, dropouts, axis=0)
                       #     flag_RR = numpy.delete(flag_RR, dropouts, axis=0)
                       # if ('LL' in polnames or 'YY' in polnames) and count_drop > 0:
                        #    amp_time_LL = numpy.delete(amp_time_LL, dropouts, axis=0)
                         #   amp_freq_LL = numpy.delete(amp_freq_LL, dropouts, axis=0)
                          #  flag_LL = numpy.delete(flag_LL, dropouts, axis=0)
                       # if ('RL' in polnames or 'XY' in polnames) and count_drop > 0:
                        #    amp_time_RL = numpy.delete(amp_time_RL, dropouts, axis=0)
                         #   amp_freq_RL = numpy.delete(amp_freq_RL, dropouts, axis=0)
                          #  flag_RL = numpy.delete(flag_RL, dropouts, axis=0)
                       # if ('LR' in polnames or 'YX' in polnames) and count_drop > 0:
                        #    amp_time_LR = numpy.delete(amp_time_LR, dropouts, axis=0)
                         #   amp_freq_LR = numpy.delete(amp_freq_LR, dropouts, axis=0)
                          #  flag_LR = numpy.delete(flag_LR, dropouts, axis=0)


                    print "Before Lovell:", test, "after Lovell:", len(times_array)







                # ZERO LEVEL DATA CODE. WHEN TELESCOPES DROPOUT FOR UNKNOWN REASONS.



                def mad_1d(array):
                    ab_dev = numpy.zeros_like(array)
                    coef = 1.4286
                    median_i = float(numpy.median(array[~numpy.isnan(array)]))
                    # median_i = float(numpy.median(array, axis=None))
                    # print median_i
                    # zero_avg = numpy.average(single[~numpy.isnan(single)])
                    ab_dev[:] = abs(array[:] - median_i)
                    # median_j = float(numpy.median(ab_dev, axis=None))
                    median_j = float(numpy.median(ab_dev[~numpy.isnan(ab_dev)]))
                    global final_mad
                    final_mad = coef*median_j
                    return final_mad



                if zero_level == 'yes':

                    print "\n CPU ", mycore, "Running Zero Level Dropout Passage on unprocessed baseline...", bline
                    #print "CPU",mycore, single

                    con = numpy.zeros([len(times_array), 2], dtype='|S12')
                    for t in xrange(len(times_array)):
                        con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]


                    zeros = []

                    for pol in polnames:
                        if pol == 'LL' or pol == 'YY':
                            mdat = numpy.ma.masked_array(amp_time_LL,numpy.isnan(amp_time_LL))
                            #print len(mdat),mdat.shape, mdat
                            single = numpy.ma.median(mdat, axis=1)*10**6
                            #break
                        elif pol == 'RR' or pol == 'XX':
                            mdat = numpy.ma.masked_array(amp_time_RR,numpy.isnan(amp_time_RR))
                            single = numpy.ma.median(mdat, axis=1)*10**6
                            #break
                        elif pol == 'RL' or pol == 'XY':
                            mdat = numpy.ma.masked_array(amp_time_RL,numpy.isnan(amp_time_RL))
                            single = numpy.ma.median(mdat, axis=1)*10**6
                            #break
                        elif pol == 'LR' or pol == 'YX':
                            mdat = numpy.ma.masked_array(amp_time_LR,numpy.isnan(amp_time_LR))
                            single = numpy.ma.median(mdat, axis=1)*10**6
                            #break


                        con = numpy.zeros([len(times_array), 2], dtype='|S12')
                        for t in xrange(len(times_array)):
                            con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]

                        #print single[0]
                        #print single[-1]
                        single = single

                        #print single[0]
                        #print single[-1]
                        #print amp

                        zero_avg = numpy.ma.average(single[~numpy.isnan(single)])
                        zero_std = numpy.ma.std(single[~numpy.isnan(single)])
                        zero_median = numpy.ma.median(single[~numpy.isnan(single)])
                        zero_mad = mad_1d(single)
                        # print 'zero_avg,zero_std,zero_median,zero_mad', zero_avg,zero_std,zero_median,zero_mad,pol
                        # print (0 + 3*zero_std), (zero_avg - 2*zero_std),pol
                        # print (zero_median - 2*zero_mad),pol
                        # print (zero_median - 2*zero_mad), (0 + 2*zero_mad),pol
                        # print (zero_median - zero_mad)

                        #print single
                        def alt_times(time):
                            time=float(time)
                            rhrs = (time*24.0)
                            hrs = math.floor(rhrs)
                            rmins = (rhrs-hrs)*60.0
                            mins = math.floor(rmins)
                            rsecs = (rmins-mins)*60.0
                            secs = math.floor(rsecs)
                            #print hrs,mins,rsecs
                            #print '%.2f:%.2f:%.2f' %(hrs,mins,rsecs)

                        for v in xrange(len(single)):
                            if single[v] < (0 + 3*zero_std) and single[v] < (zero_avg - 2*zero_std):
                                #print 'time1:', times_array[v],v; alt_times(times_array[v])
                                if v not in zeros:
                                    zeros.append(v)
                            if single[v] < (zero_median - 3*zero_mad):
                                #print 'time2:', times_array[v],v; alt_times(times_array[v])
                                if v not in zeros:
                                    zeros.append(v)
                            if single[v] < (zero_median - 3*zero_mad) and single[v] < (0 + 2*zero_mad):
                                #print 'time3:', times_array[v],v; alt_times(times_array[v])
                                if v not in zeros:
                                    zeros.append(v)

                            ## To get the Cambridge zero levels which affect 75% of observation?
                            #if single[v] < (0 + 5*zero_std):
                            #if single[v] < (zero_median - zero_mad):
                                #if v not in zeros:
                                    #zeros.append(v)

                        print "Dropouts flagged:", pol, 'IF'+str(j+1), len(zeros)
                        #print sorted(zeros), len(zeros)

                        
                        if pol == 'RR' or pol == 'XX' and len(zeros) > 0:
                            for i in xrange(len(single)):
                                if i in zeros:
                                    amp_time_RR[i,:] = zero_median
                                    amp_freq_RR[i,:] = zero_median
                            # print  bline,str(j+1),pol+':',sorted(zeros), len(zeros)
                        if pol == 'LL' or pol == 'YY' and len(zeros) > 0:
                            for i in xrange(len(single)):
                                if i in zeros:
                                    amp_time_LL[i,:] = zero_median
                                    amp_freq_LL[i,:] = zero_median
                            # print  bline,str(j+1),pol+':',sorted(zeros), len(zeros)
                        if pol == 'RL' or pol == 'XY' and len(zeros) > 0:
                            for i in xrange(len(single)):
                                if i in zeros:
                                    amp_time_RL[i,:] = zero_median
                                    amp_freq_RL[i,:] = zero_median
                            # print  bline,str(j+1),pol+':',sorted(zeros), len(zeros)
                        if pol == 'LR' or pol == 'YX' and len(zeros) > 0:
                            for i in xrange(len(single)):
                                if i in zeros:
                                    amp_time_LR[i,:] = zero_median
                                    amp_freq_LR[i,:] = zero_median
                            # print bline,str(j+1),pol+':',sorted(zeros), len(zeros)


                    
                    times_array = numpy.delete(times_array, zeros)
                    if ('RR' in polnames or 'XX' in polnames) and len(zeros) > 0:
                        amp_time_RR = numpy.delete(amp_time_RR, zeros, axis=0)
                        amp_freq_RR = numpy.delete(amp_freq_RR, zeros, axis=0)
                        flag_RR = numpy.delete(flag_RR, zeros, axis=0)
                    if ('LL' in polnames or 'YY' in polnames) and len(zeros) > 0:
                        amp_time_LL = numpy.delete(amp_time_LL, zeros, axis=0)
                        amp_freq_LL = numpy.delete(amp_freq_LL, zeros, axis=0)
                        flag_LL = numpy.delete(flag_LL, zeros, axis=0)
                    if ('RL' in polnames or 'XY' in polnames) and len(zeros) > 0:
                        amp_time_RL = numpy.delete(amp_time_RL, zeros, axis=0)
                        amp_freq_RL = numpy.delete(amp_freq_RL, zeros, axis=0)
                        flag_RL = numpy.delete(flag_RL, zeros, axis=0)
                    if ('LR' in polnames or 'YX' in polnames) and len(zeros) > 0:
                        amp_time_LR = numpy.delete(amp_time_LR, zeros, axis=0)
                        amp_freq_LR = numpy.delete(amp_freq_LR, zeros, axis=0)
                        flag_LR = numpy.delete(flag_LR, zeros, axis=0)



                    for pol in polnames:
                        if len(zeros) > 0:
                            #print pol, len(zeros)
                            plname = str(name) + "__" + str(bline) + "__IF" + str(j+1)+ "__" + pol + "__zeros_dummy"
                            # infolist = numpy.array([[source], [con], [dropouts]])
                            pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                            pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                            pickle.dump(zeros, open(path2folder + str(plname) + ".zeros", "wb"))
                    zeros_bline.append(bline)

               


                # Sum Threshold flagging sequence:


                def mad(array):
                    ab_dev = numpy.zeros_like(array)
                    coef = 1.4286
                    #median_i = float(numpy.median(array, axis=None))
                    median_i = float(numpy.median(array[~numpy.isnan(array)], axis=None))
                    ab_dev[:,:] = abs(array[:,:] - median_i)
                    #median_j = float(numpy.median(ab_dev, axis=None))
                    median_j = float(numpy.median(ab_dev[~numpy.isnan(array)], axis=None))
                    global final_mad
                    final_mad = coef*median_j
                    #print "The Median Absolute Deviation for this baseline is %f" % final_mad
                    return final_mad


                global_total = 0
                global_time = 0
                global_freq = 0




                def flagger(amp_time, amp_freq, flag):
                    if flagging_options == 'default':
                        aggressiveness_first_run = 25
                        max_subset_first_run = 32
                        aggressiveness_second_run = 25
                        max_subset_second_run = 256
                        rho = 1.5
                    else:
                        aggressiveness_first_run = parameters[0]
                        max_subset_first_run = parameters[1]
                        aggressiveness_second_run = parameters[2]
                        max_subset_second_run = parameters[3]
                        rho = parameters[4]
                        kickout_sigma_level = parameters[5]

                    # Initial run with subset size = 1 to remove v. strong RFI
                    # Note only one direction needed window = 1

                    # Need to consider any values the correlator has flagged. Written them in SERPent as 'nan'
                    for col in xrange(amp_time.shape[1]):
                        for row in xrange(amp_time.shape[0]):
                            if str(amp_time[row, col]) == "nan":
                                flag[row, col] = -1.0

                    print flag

                    median = numpy.median(amp_time[~numpy.isnan(amp_time)])
                    chi_1 = median + (20*mad(amp_time))

                    i = 1
                    if flagging_options == 'default':
                        if freq_band(frequency) in ('L Band' or 'S Band'):
                            kickout = median + 3.0*final_mad
                            #print "Kickout aggressive for low frequency..."
                        else:
                            kickout = median + 3.0*final_mad
                            #print "Kickout less aggressive for high frequency..."
                    # User defined kickout sigma level:
                    else:
                        kickout = median + kickout_sigma_level*final_mad

                    while i < amp_time.shape[0]+1 and i <= 1:
                        timing = time.time()
                        limit = chi_1/(rho**(math.log(i,2)))
                        # The below condition is to avoid flagging good data if the thresholding gets too close to the mean
                        if limit <= kickout:
                            break
                        for col in xrange(amp_time.shape[1]):# all the channels
                            for row in xrange(amp_time.shape[0]):
                                if flag[row,col] == 1.0:
                                    amp_time[row,col] = limit
                                    amp_freq[row, col] = limit
                            for row in xrange(1, amp_time.shape[0]-(i-2)):
                                if i == 1:
                                    if amp_time[row-1,col] > limit:
                                        if not flag[row-1, col] < 0.0:
                                            flag[row-1, col] = 1.0

                        #print "i =", i, "time (secs):", time.time()-timing
                        i *= 2
                    #
                    count_i = 0
                    # count_i = numpy.count_nonzero(flag)
                    count_i = (flag > 0.0).sum()
                    old_count = (flag != 0.0).sum()
                    print "Flagged with initial run:", count_i, 'old_count', old_count


                    # Run the full Threshold algorithm for the first time.

                    # Flagged vis will be at lowest threshold for the following statistics:

                    hm = time.time()
                    median = numpy.median(amp_time[~numpy.isnan(amp_time)])
                    chi_1 = median + (aggressiveness_first_run*mad(amp_time))

                    # Reset the amplitude arrays for the new statistics.
                    # Flagged visibilities get set to the first threshold level because the one within the while loop is for setting limits within that sumthreshold run

                    for col in xrange(flag.shape[1]):
                        for row in xrange(flag.shape[0]):
                            if flag[row,col] != 0.0:
                                amp_time[row,col] = chi_1
                                amp_freq[row,col] = chi_1


                    i = 1
                    #print "\n FLAGGING IN TIME..."
                    while i < amp_time.shape[0]+1 and i < max_subset_first_run+1:
                        timing = time.time()
                        limit = chi_1/(rho**(math.log(i,2)))
                        # The below condition is to avoid flagging good data if the thresholding gets too close to the mean
                        if limit <= kickout:
                            break
                        for col in xrange(amp_time.shape[1]):# all the channels
                            for row in xrange(amp_time.shape[0]):
                                if flag[row,col] == 2.0:
                                    amp_time[row,col] = limit
                                if flag[row,col] == -2.0:
                                    amp_time[row,col] = limit
                            for row in xrange(1, amp_time.shape[0]-(i-2)):
                                if i == 1:
                                    if amp_time[row-1,col] > limit:
                                        if not flag[row-1, col] < 0.0:
                                            flag[row-1, col] = 2.0
                                        elif flag[row-1, col] < 0.0:
                                            flag[row-1, col] = -2.0
                                            
                                else:
                                    sum = amp_time[row-1:row-1+i, col].sum()
                                    if (sum/i) > limit:
                                        flag[row-1:row-1+i,col]=numpy.where(flag[row-1:row-1+i,col] < 0.0, flag[row-1:row-1+i,col], 2.0)
                                        flag[row-1:row-1+i,col]=numpy.where(flag[row-1:row-1+i,col] >= 0.0, flag[row-1:row-1+i,col], -2.0)
                                        # flag[row-1:row-1+i, col] = 2.0
                        #print "i =", i, "time (secs):", time.time()-timing
                        i *= 2
                    #
                    count_t = 0
                    for col in xrange(flag.shape[1]):# all the channels
                        for row in xrange(flag.shape[0]):
                            if flag[row,col] == 2.0:
                                count_t += 1
                    #print "Flagged time:", count_t
                    i = 1
                    #print "\n FLAGGING IN FREQUENCY..."
                    while i < amp_freq.shape[1]+1 and i < max_subset_first_run+1:
                        timing = time.time()
                        limit = chi_1/(rho**(math.log(i,2)))
                        # The below condition is to avoid flagging good data if the thresholding gets too close to the mean
                        if limit <= kickout:
                            break
                        for row in xrange(amp_freq.shape[0]):
                            for col in xrange(amp_freq.shape[1]):
                                if flag[row,col] == 3.0:
                                    amp_freq[row,col] = limit
                                if flag[row,col] == -3.0:
                                    amp_freq[row,col] = limit
                            for col in xrange(1, amp_freq.shape[1]-(i-2)):
                                if i == 1:
                                    if amp_freq[row, col-1] > limit:
                                        if not flag[row, col-1] < 0.0:
                                            flag[row, col-1] = 3.0
                                        elif flag[row, col-1] < 0.0:
                                            flag[row, col-1] = -3.0
                                else:
                                    sum = amp_freq[row, col-1:col-1+i].sum()
                                    if (sum/i) > limit:
                                        flag[row,col-1:col-1+i]=numpy.where(flag[row,col-1:col-1+i] < 0.0, flag[row,col-1:col-1+i], 3.0)
                                        flag[row,col-1:col-1+i]=numpy.where(flag[row,col-1:col-1+i] >= 0.0, flag[row,col-1:col-1+i], -3.0)
                                        # flag[row, col-1:col-1+i] = 3.0
                        #print "i =", i, "time (secs):", time.time()-timing
                        i *= 2    # subset spacing
                    #
                    count_f = 0
                    for col in xrange(flag.shape[1]):# all the channels
                        for row in xrange(flag.shape[0]):
                            if flag[row,col] == 3.0:
                                count_f += 1
                    #print "Flagged freq:", count_f

                    # second run of flagger removing already found strong RFI flags and calculating thresholds with new stats
                    #if (really high max value a certain magnitude above average/median then run flagger again removing the strong RFIs from the arrays for the sumThreshold...)


                    median = numpy.median(amp_time)
                    maxvalue = numpy.amax(amp_time)
                    chi_1 = median + (aggressiveness_second_run*mad(amp_time))


                    for col in xrange(flag.shape[1]):
                        for row in xrange(flag.shape[0]):
                            if flag[row,col] != 0.0:
                                amp_time[row,col] = chi_1
                                amp_freq[row,col] = chi_1


                    i = 1
                    #print "\n FLAGGING IN TIME..."
                    if flagging_options == 'default':
                        if freq_band(frequency) in ('L Band' or 'S Band'):
                            kickout = median + 3.0*final_mad
                            #print "Kickout aggressive for low frequency..."
                        else:
                            kickout = median + 3.0*final_mad
                            #print "Kickout less aggressive for high frequency..."
                    # User defined kickout sigma level:
                    else:
                        kickout = median + kickout_sigma_level*final_mad
                    while i < amp_time.shape[0]+1 and i < max_subset_second_run+1:
                        if maxvalue > kickout:
                            if count_t > 0 or count_f > 0:
                                pass
                            else:
                                print " Nothing flagged. Skipping run..."
                                break
                        else:
                            print " Max value < kickout. Skipping run..."
                            break
                        timing = time.time()
                        limit = chi_1/(rho**(math.log(i,2)))
                        # The below condition is to avoid flagging good data if the thresholding gets too close to the mean
                        if limit <= kickout:
                            break
                        for col in xrange(amp_time.shape[1]):# all the channels
                            for row in xrange(amp_time.shape[0]):
                                if flag[row,col] == 4.0:
                                    amp_time[row,col] = limit
                                if flag[row,col] == -4.0:
                                    amp_time[row,col] = limit
                            for row in xrange(1, amp_time.shape[0]-(i-2)):
                                if i == 1:
                                    if amp_time[row-1,col] > limit:
                                        if not flag[row-1, col] < 0.0:
                                            flag[row-1, col] = 4.0
                                        elif flag[row-1, col] < 0.0:
                                            flag[row-1, col] = -4.0
                                else:
                                    sum = amp_time[row-1:row-1+i, col].sum()
                                    if (sum/i) > limit:
                                        flag[row-1:row-1+i,col]=numpy.where(flag[row-1:row-1+i,col] < 0.0, flag[row-1:row-1+i,col], 4.0)
                                        flag[row-1:row-1+i,col]=numpy.where(flag[row-1:row-1+i,col] >= 0.0, flag[row-1:row-1+i,col], -4.0)
                                        # flag[row-1:row-1+i, col] = 4.0
                        #print "i =", i, "time (secs):", time.time()-timing
                        i *= 2
                    #
                    for col in xrange(flag.shape[1]):# all the channels
                        for row in xrange(flag.shape[0]):
                            if flag[row,col] == 4.0:
                                count_t += 1
                    print "Flagged time:", count_t
                    i = 1
                    #print "\n FLAGGING IN FREQUENCY..."
                    if flagging_options == 'default':
                        if freq_band(frequency) in ('L Band' or 'S Band'):
                            kickout = median + 3.0*final_mad
                            #print "Kickout aggressive for low frequency..."
                        else:
                            kickout = median + 3.0*final_mad
                            #print "Kickout less aggressive for high frequency..."
                    # User defined kickout sigma level:
                    else:
                        kickout = median + kickout_sigma_level*final_mad
                    while i < amp_freq.shape[1]+1 and i < max_subset_second_run+1:
                        if maxvalue > kickout:
                            if count_t > 0 or count_f > 0:
                                pass
                            else:
                                print " Nothing flagged. Skipping run..."
                                break
                        else:
                            print " Max value < kickout. Skipping run..."
                            break
                        timing = time.time()
                        limit = chi_1/(rho**(math.log(i,2)))
                        # The below condition is to avoid flagging good data if the thresholding gets too close to the mean
                        if limit <= kickout:
                            break
                        for row in xrange(amp_freq.shape[0]):
                            for col in xrange(amp_freq.shape[1]):
                                if flag[row,col] == 5.0:
                                    amp_freq[row,col] = limit
                                if flag[row,col] == -5.0:
                                    amp_freq[row,col] = limit
                            for col in xrange(1, amp_freq.shape[1]-(i-2)):
                                if i == 1:
                                    if amp_freq[row, col-1] > limit:
                                        if not flag[row, col-1] < 0.0:
                                            flag[row, col-1] = 5.0
                                        elif flag[row, col-1] < 0.0:
                                            flag[row, col-1] = -5.0
                                else:
                                    sum = amp_freq[row, col-1:col-1+i].sum()
                                    if (sum/i) > limit:
                                        flag[row,col-1:col-1+i]=numpy.where(flag[row,col-1:col-1+i] < 0.0, flag[row,col-1:col-1+i], 5.0)
                                        flag[row,col-1:col-1+i]=numpy.where(flag[row,col-1:col-1+i] >= 0.0, flag[row,col-1:col-1+i], -5.0)
                                        # flag[row, col-1:col-1+i] = 5.0
                        #print "i =", i, "time (secs):", time.time()-timing
                        i *= 2    # subset spacing
                    #
                    for col in xrange(flag.shape[1]):# all the channels
                        for row in xrange(flag.shape[0]):
                            if flag[row,col] == 5.0:
                                count_f += 1
                    print "Flagged freq:", count_f
                    # flagged_total = numpy.count_nonzero(flag)
                    flagged_total = (flag > 0.0).sum()
                    print "Total number of flagged amplitudes: %i" % flagged_total
                    global global_total
                    global_total += flagged_total
                    global global_time
                    global_time += count_t
                    global global_freq
                    global_freq += count_f
                    print "CPU", mycore, "time (secs):", time.time() - hm

                    print flag
                    return (flag)

                
                # Execute the Flagger for all polarizations and baselines
                if 'RR' in polnames or 'XX' in polnames:
                    if 'RR' in polnames:
                        pol = 'RR'
                    if 'XX' in polnames:
                        pol = 'XX'
                    print " \n CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, uvdata.name)
                    if amp_time_RR.shape[0] > 0 and amp_freq_RR.shape[0] > 0:
                        print "\n CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, j+1)
                        flagger(amp_time_RR, amp_freq_RR, flag_RR)
                        # only need one direction test as the other will also be equal to 0.
                    del amp_time_RR
                    del amp_freq_RR

                if 'LL' in polnames or 'YY' in polnames:
                    if 'LL' in polnames:
                        pol = 'LL'
                    if 'YY' in polnames:
                        pol = 'YY'
                    print "\n CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, uvdata.name)
                    if amp_time_LL.shape[0] > 0 and amp_freq_LL.shape[0] > 0:
                        print "\n CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, j+1)
                        flagger(amp_time_LL, amp_freq_LL, flag_LL)
                    del amp_time_LL
                    del amp_freq_LL

                if 'RL' in polnames or 'XY' in polnames:
                    if 'RL' in polnames:
                        pol = 'RL'
                    if 'XY' in polnames:
                        pol = 'XY'
                    print "\n CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, uvdata.name)
                    if amp_time_RL.shape[0] > 0 and amp_freq_RL.shape[0] > 0:
                        print "\n CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, j+1)
                        flagger(amp_time_RL, amp_freq_RL, flag_RL)
                    del amp_time_RL
                    del amp_freq_RL

                if 'LR' in polnames or 'YX' in polnames:
                    if 'LR' in polnames:
                        pol = 'LR'
                    if 'YX' in polnames:
                        pol = 'YX'
                    print "\n CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, uvdata.name)
                    if amp_time_LR.shape[0] > 0 and amp_freq_LR.shape[0] > 0:
                        print "\n CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, j+1)
                        flagger(amp_time_LR, amp_freq_LR, flag_LR)
                    del amp_time_LR
                    del amp_freq_LR
                

                def find_source_index(source_list):
                    for s in source_list:
                        if source_list[s] == uvdata.name:
                            source = s
                            break
                    return source


                print "CPU %i has finished flagging Source: %s, IF: %i, Baseline: %s" % (mycore, name, j+1, bline)

                # Now combine all the stokes flag arrays together into one and when writing the flag table flag for all stokes. i.e. '1111'
                exists = 0

                if coadd_polarization_flags == 'yes':

                    # if nan_array is 0 do not add flags to master flag array as correlator already has all vis weights to 0... could still use other pols for continuum obs...
                    if ('RR' in polnames or 'XX' in polnames) and nan_array_RR == 0:
                        flag_ALL = flag_RR
                        exists = 1

                    if ('LL' in polnames or 'YY' in polnames) and nan_array_LL == 0:
                        if exists == 1:
                            flag_ALL = flag_ALL + flag_LL
                        else:
                            flag_ALL = flag_LL
                            exists = 1

                    if ('RL' in polnames or 'XY' in polnames) and nan_array_RL == 0:
                        if exists == 1:
                            flag_ALL = flag_ALL + flag_RL
                        else:
                            flag_ALL = flag_RL
                            exists = 1
                    if ('LR' in polnames or 'YX' in polnames) and nan_array_LR == 0:
                        if exists == 1:
                            flag_ALL = flag_ALL + flag_LR
                        else:
                            flag_ALL = flag_LR
                            exists = 1
                    try:
                        flag_ALL
                    except NameError:
                        flag_ALL = numpy.zeros([ntimescans[bline], nchan])

                    # num_flagged_vis = numpy.count_nonzero(flag_ALL)
                    #num_flagged_vis = (flag != 0).sum()
                    num_flagged_vis = (flag_ALL != 0).sum()
                    print "Total number of flags after co-adding all stokes: %i for IF %i" % (num_flagged_vis, j+1)
                    tsum = 0
                    tsum += ntimescans[bline]
                    total_tscans += tsum
                    #print "Total number of amplitudes: %i" % (tsum*nif*nchan*npol)
                    print "Percentage of amplitudes flagged for IF %i, baseline %s: %f%%" % (j+1, bline, (100.0*num_flagged_vis/ (tsum*nif*nchan*npol)))
                    final_num_vis += num_flagged_vis


                    # Pickle the flag_ALL array into a text file and save name and properties in a dictionary.
                    print "\n Source : %s, IF: %s, Baseline: %s" % (name, j+1, bline)
                    # need to count_nonzero of the flag_ALL array.
                    # then create text file with object name, bline, if and load array in
                    pname = str(name) + "__IF" + str(j+1) + "__" + str(bline) + "__all-stokes"
                    # print "pname", pname
                    #if numpy.count_nonzero(flag_ALL) > 0:
                    pickle.dump(flag_ALL, open(path2folder + str(pname)+".p", "wb"))
                    # save source num, ifnum, bline and time info with pickle name:
                    pickle_list[pname] = [source, j+1, bline, times_array]
                    pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))


                # do not co add polarizations
                else:
                    print 'Not co-adding stokes flags!'
                                        
                    if ('RR' in polnames or 'XX' in polnames) and nan_array_RR == 0:
                        if 'RR' in polnames:
                            pol = 'RR'
                        if 'XX' in polnames:
                            pol = 'XX'
                        print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, j+1, bline, pol)
                        pname = str(name) + "__IF" + str(j+1) + "__" + str(bline) + "__" + pol
                        pickle.dump(flag_RR, open(path2folder + str(pname)+".p", "wb"))
                        pickle_list[pname] = [source, j+1, bline, times_array]
                        pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                        #print len(flag_RR)

                    if ('LL' in polnames or 'YY' in polnames) and nan_array_LL == 0:
                        if 'LL' in polnames:
                            pol = 'LL'
                        if 'YY' in polnames:
                            pol = 'YY'
                        print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, j+1, bline, pol)
                        pname = str(name) + "__IF" + str(j+1) + "__" + str(bline) + "__" + pol
                        pickle.dump(flag_LL, open(path2folder + str(pname)+".p", "wb"))
                        pickle_list[pname] = [source, j+1, bline, times_array]
                        pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                        # print len(flag_LL)

                    if ('RL' in polnames or 'XY' in polnames) and nan_array_RL == 0:
                        if 'RL' in polnames:
                            pol = 'RL'
                        if 'XY' in polnames:
                            pol = 'XY'
                        print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, j+1, bline, pol)
                        pname = str(name) + "__IF" + str(j+1) + "__" + str(bline) + "__" + pol
                        pickle.dump(flag_RL, open(path2folder + str(pname)+".p", "wb"))
                        pickle_list[pname] = [source, j+1, bline, times_array]
                        pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                        # print len(flag_RL)

                    if ('LR' in polnames or 'YX' in polnames) and nan_array_LR == 0:
                        if 'LR' in polnames:
                            pol = 'LR'
                        if 'YX' in polnames:
                            pol = 'YX'
                        print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, j+1, bline, pol)
                        pname = str(name) + "__IF" + str(j+1) + "__" + str(bline) + "__" + pol
                        pickle.dump(flag_LR, open(path2folder + str(pname)+".p", "wb"))
                        pickle_list[pname] = [source, j+1, bline, times_array]
                        pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                        # print len(flag_LR)






###############################################################################
# DON'T INDENT BELOW THIS LINE!!!!!!!!!!!!!!!!!!!!!!!!!!
###############################################################################

    print "\n CPU", mycore, "has finished all jobs...\n"
# End Parallelization and wait for all processes to finish before reading everything into fg file.

    if 0 in pid:
# Childs die silently:
        os._exit(0)

    elif NCPU>1:
# Parent waits for childs and continues:
        for j,idi in enumerate(pid):
            os.waitpid(idi,0)
            print 'CPU',j+1,' died OK.'



# Header and FG text file read definitions


def fgheader(path2folder, no_rows):
    # Need number of rows which we're going to enter.
    global fgfile
    global fg_row
    fgfile = open(path2folder+experiment+'.fg', 'wr')#operation = 'wr'
    fgfile.flush()    # flush prints everything from the file onto screen
    # command below means print to this file
    print >> fgfile, "XTENSION= 'A3DTABLE'           / extension type\n",'BITPIX  =                    8 / printable ASCII codes\n','NAXIS   =                    2 / Table is a matrix\n','NAXIS1  =                  132 / Max. no. of characters/pass\n','NAXIS2  =             ','%7i' % no_rows,'/ Number of entries in table\n','PCOUNT  =                    0 / Random parameter count\n','GCOUNT  =                    1 / Group count\n','NOPASS  =                    1 / Number of passes thru table\n','TFIELDS =                    9 / Number of fields in each row\n',"EXTNAME = 'AIPS FG '           / AIPS table file\n",'EXTVER  =                    1 / Version Number of table\n','TBCOL1  =                    9 / Starting char. pos. of field\n',"TFORM1  = 'I11     '           / Fortran format of field  1\n",'TFDIM1  =                    1 / Dimension of field  1\n',"TTYPE1  = 'SOURCE  '           / type (heading) of field  1\n","TUNIT1  = '        '           / physical units of field  1\n",'TBCOL2  =                   20 / Starting char. pos. of field\n',"TFORM2  = 'I11     '           / Fortran format of field  2\n",'TFDIM2  =                    1 / Dimension of field  2\n',"TTYPE2  = 'SUBARRAY                '           / type (heading) of field  2\n","TUNIT2  = '        '           / physical units of field  2\n",'TBCOL3  =                   31 / Starting char. pos. of field\n',"TFORM3  = 'I11     '           / Fortran format of field  3\n",'TFDIM3  =                    1 / Dimension of field  3\n',"TTYPE3  = 'FREQ ID                 '           / type (heading) of field  3\n","TUNIT3  = '        '           / physical units of field  3\n",'TBCOL4  =                   42 / Starting char. pos. of field\n',"TFORM4  = 'I11     '           / Fortran format of field  4\n",'TFDIM4  =                    2 / Dimension of field  4\n',"TTYPE4  = 'ANTS                    '           / type (heading) of field  4\n","TUNIT4  = '        '           / physical units of field  4\n",'TBCOL5  =                   53 / Starting char. pos. of field\n',"TFORM5  = 'E18.8   '           / Fortran format of field  5\n",'TFDIM5  =                    2 / Dimension of field  5\n',"TTYPE5  = 'TIME RANGE              '           / type (heading) of field  5\n","TUNIT5  = 'DAYS    '           / physical units of field  5\n",'TBCOL6  =                   71 / Starting char. pos. of field\n',"TFORM6  = 'I9     '           / Fortran format of field  6\n",'TFDIM6  =                    2 / Dimension of field  6\n',"TTYPE6  = 'IFS                     '           / type (heading) of field  6\n","TUNIT6  = '        '           / physical units of field  6\n",'TBCOL7  =                   80 / Starting char. pos. of field\n',"TFORM7  = 'I9     '           / Fortran format of field  7\n",'TFDIM7  =                    2 / Dimension of field  7\n',"TTYPE7  = 'CHANS                   '           / type (heading) of field  7\n","TUNIT7  = '        '           / physical units of field  7\n",'TBCOL8  =                   89 / Starting char. pos. of field\n',"TFORM8  = 'X4      '           / Fortran format of field  8\n",'TFDIM8  =                    4 / Dimension of field  8\n',"TTYPE8  = 'PFLAGS                  '           / type (heading) of field  8\n","TUNIT8  = '        '           / physical units of field  8\n",'TBCOL9  =                  101 / Starting char. pos. of field\n',"TFORM9  = 'A24     '           / Fortran format of field  9\n",'TFDIM9  =                   26 / Dimension of field  9\n',"TTYPE9  = 'REASON                  '           / type (heading) of field  9\n","TUNIT9  = '        '           / physical units of field  9\n",'END\n','COL. NO.      1          2          3          4            5              6          7          8        9\n','     ROW   SOURCE     SUBARRAY   FREQ ID    ANTS         TIME RAN       IFS        CHANS      PFLAGS   REASON\n','  NUMBER                                                 DAYS\n','***BEGIN*PASS***'
#    dummy = open(path2folder+'dummy.fg', 'wr')
    dummy = open(path2folder+'dummy.fg', 'r')
    for line in dummy:
        line = line.rstrip('\n')
        print >> fgfile, line
    fgfile.close()




def lovell_header(path2folder, total_rows):
    # Need number of rows which we're going to enter.
    lovell_dropouts = open(path2folder+experiment+'_lovell_dropouts.fg', 'wr')#operation = 'wr'
    lovell_dropouts.flush()    # flush prints everything from the file onto screen
    # command below means print to this file
    print >> lovell_dropouts, "XTENSION= 'A3DTABLE'           / extension type\n",'BITPIX  =                    8 / printable ASCII codes\n','NAXIS   =                    2 / Table is a matrix\n','NAXIS1  =                  132 / Max. no. of characters/pass\n','NAXIS2  =             ','%7i' % total_rows,'/ Number of entries in table\n','PCOUNT  =                    0 / Random parameter count\n','GCOUNT  =                    1 / Group count\n','NOPASS  =                    1 / Number of passes thru table\n','TFIELDS =                    9 / Number of fields in each row\n',"EXTNAME = 'AIPS FG '           / AIPS table file\n",'EXTVER  =                    1 / Version Number of table\n','TBCOL1  =                    9 / Starting char. pos. of field\n',"TFORM1  = 'I11     '           / Fortran format of field  1\n",'TFDIM1  =                    1 / Dimension of field  1\n',"TTYPE1  = 'SOURCE  '           / type (heading) of field  1\n","TUNIT1  = '        '           / physical units of field  1\n",'TBCOL2  =                   16 / Starting char. pos. of field\n',"TFORM2  = 'I11     '           / Fortran format of field  2\n",'TFDIM2  =                    1 / Dimension of field  2\n',"TTYPE2  = 'SUBARRAY                '           / type (heading) of field  2\n","TUNIT2  = '        '           / physical units of field  2\n",'TBCOL3  =                   28 / Starting char. pos. of field\n',"TFORM3  = 'I11     '           / Fortran format of field  3\n",'TFDIM3  =                    1 / Dimension of field  3\n',"TTYPE3  = 'FREQ ID                 '           / type (heading) of field  3\n","TUNIT3  = '        '           / physical units of field  3\n",'TBCOL4  =                   40 / Starting char. pos. of field\n',"TFORM4  = 'I11     '           / Fortran format of field  4\n",'TFDIM4  =                    2 / Dimension of field  4\n',"TTYPE4  = 'ANTS                    '           / type (heading) of field  4\n","TUNIT4  = '        '           / physical units of field  4\n",'TBCOL5  =                   51 / Starting char. pos. of field\n',"TFORM5  = 'E18.8   '           / Fortran format of field  5\n",'TFDIM5  =                    2 / Dimension of field  5\n',"TTYPE5  = 'TIME RANGE              '           / type (heading) of field  5\n","TUNIT5  = 'DAYS    '           / physical units of field  5\n",'TBCOL6  =                   71 / Starting char. pos. of field\n',"TFORM6  = 'I9     '           / Fortran format of field  6\n",'TFDIM6  =                    2 / Dimension of field  6\n',"TTYPE6  = 'IFS                     '           / type (heading) of field  6\n","TUNIT6  = '        '           / physical units of field  6\n",'TBCOL7  =                   80 / Starting char. pos. of field\n',"TFORM7  = 'I9     '           / Fortran format of field  7\n",'TFDIM7  =                    2 / Dimension of field  7\n',"TTYPE7  = 'CHANS                   '           / type (heading) of field  7\n","TUNIT7  = '        '           / physical units of field  7\n",'TBCOL8  =                   89 / Starting char. pos. of field\n',"TFORM8  = 'X4      '           / Fortran format of field  8\n",'TFDIM8  =                    4 / Dimension of field  8\n',"TTYPE8  = 'PFLAGS                  '           / type (heading) of field  8\n","TUNIT8  = '        '           / physical units of field  8\n",'TBCOL9  =                  104 / Starting char. pos. of field\n',"TFORM9  = 'A24     '           / Fortran format of field  9\n",'TFDIM9  =                   26 / Dimension of field  9\n',"TTYPE9  = 'REASON                  '           / type (heading) of field  9\n","TUNIT9  = '        '           / physical units of field  9\n",'END\n','COL. NO.      1          2          3          4            5              6          7          8        9\n','     ROW   SOURCE     SUBARRAY   FREQ ID    ANTS         TIME RAN       IFS        CHANS      PFLAGS   REASON\n','  NUMBER                                                 DAYS\n','***BEGIN*PASS***'
#    dummy = open(path2folder+'dummy.fg', 'wr')
    try:
        lovell_dummy = open(path2folder+experiment+'_lovell_dummy.fg', 'r')
        for line in lovell_dummy:
            line = line.rstrip('\n')
            print >> lovell_dropouts, line
    except:
        print "No Lovell dropouts to flag..."
    lovell_dropouts.close()




def zeros_header(path2folder, total_rows):
    # Need number of rows which we're going to enter.
    zeros_dropouts = open(path2folder+experiment+'_zeros_dropouts.fg', 'wr')#operation = 'wr'
    zeros_dropouts.flush()    # flush prints everything from the file onto screen
    # command below means print to this file
    print >> zeros_dropouts, "XTENSION= 'A3DTABLE'           / extension type\n",'BITPIX  =                    8 / printable ASCII codes\n','NAXIS   =                    2 / Table is a matrix\n','NAXIS1  =                  132 / Max. no. of characters/pass\n','NAXIS2  =             ','%7i' % total_rows,'/ Number of entries in table\n','PCOUNT  =                    0 / Random parameter count\n','GCOUNT  =                    1 / Group count\n','NOPASS  =                    1 / Number of passes thru table\n','TFIELDS =                    9 / Number of fields in each row\n',"EXTNAME = 'AIPS FG '           / AIPS table file\n",'EXTVER  =                    1 / Version Number of table\n','TBCOL1  =                    9 / Starting char. pos. of field\n',"TFORM1  = 'I11     '           / Fortran format of field  1\n",'TFDIM1  =                    1 / Dimension of field  1\n',"TTYPE1  = 'SOURCE  '           / type (heading) of field  1\n","TUNIT1  = '        '           / physical units of field  1\n",'TBCOL2  =                   20 / Starting char. pos. of field\n',"TFORM2  = 'I11     '           / Fortran format of field  2\n",'TFDIM2  =                    1 / Dimension of field  2\n',"TTYPE2  = 'SUBARRAY                '           / type (heading) of field  2\n","TUNIT2  = '        '           / physical units of field  2\n",'TBCOL3  =                   31 / Starting char. pos. of field\n',"TFORM3  = 'I11     '           / Fortran format of field  3\n",'TFDIM3  =                    1 / Dimension of field  3\n',"TTYPE3  = 'FREQ ID                 '           / type (heading) of field  3\n","TUNIT3  = '        '           / physical units of field  3\n",'TBCOL4  =                   42 / Starting char. pos. of field\n',"TFORM4  = 'I11     '           / Fortran format of field  4\n",'TFDIM4  =                    2 / Dimension of field  4\n',"TTYPE4  = 'ANTS                    '           / type (heading) of field  4\n","TUNIT4  = '        '           / physical units of field  4\n",'TBCOL5  =                   53 / Starting char. pos. of field\n',"TFORM5  = 'E18.8   '           / Fortran format of field  5\n",'TFDIM5  =                    2 / Dimension of field  5\n',"TTYPE5  = 'TIME RANGE              '           / type (heading) of field  5\n","TUNIT5  = 'DAYS    '           / physical units of field  5\n",'TBCOL6  =                   71 / Starting char. pos. of field\n',"TFORM6  = 'I9     '           / Fortran format of field  6\n",'TFDIM6  =                    2 / Dimension of field  6\n',"TTYPE6  = 'IFS                     '           / type (heading) of field  6\n","TUNIT6  = '        '           / physical units of field  6\n",'TBCOL7  =                   80 / Starting char. pos. of field\n',"TFORM7  = 'I9     '           / Fortran format of field  7\n",'TFDIM7  =                    2 / Dimension of field  7\n',"TTYPE7  = 'CHANS                   '           / type (heading) of field  7\n","TUNIT7  = '        '           / physical units of field  7\n",'TBCOL8  =                   89 / Starting char. pos. of field\n',"TFORM8  = 'X4      '           / Fortran format of field  8\n",'TFDIM8  =                    4 / Dimension of field  8\n',"TTYPE8  = 'PFLAGS                  '           / type (heading) of field  8\n","TUNIT8  = '        '           / physical units of field  8\n",'TBCOL9  =                  101 / Starting char. pos. of field\n',"TFORM9  = 'A24     '           / Fortran format of field  9\n",'TFDIM9  =                   26 / Dimension of field  9\n',"TTYPE9  = 'REASON                  '           / type (heading) of field  9\n","TUNIT9  = '        '           / physical units of field  9\n",'END\n','COL. NO.      1          2          3          4            5              6          7          8        9\n','     ROW   SOURCE     SUBARRAY   FREQ ID    ANTS         TIME RAN       IFS        CHANS      PFLAGS   REASON\n','  NUMBER                                                 DAYS\n','***BEGIN*PASS***'
#    dummy = open(path2folder+'dummy.fg', 'wr')
    try:
        zeros_dummy = open(path2folder+experiment+'_zeros_dummy.fg', 'r')
        for line in zeros_dummy:
            line = line.rstrip('\n')
            print >> zeros_dropouts, line
    except:
        print "No Zero Amplitude dropouts to flag..."
    zeros_dropouts.close()



def combined_header(path2folder, combined_rows, no_rows, baseline_removal, nif, nchan): #lovellmk2):
    # Need number of rows which we're going to enter.
    #global fgfile
    #global fg_row
    combined_file = open(path2folder+experiment+'_combined.fg', 'wr')#operation = 'wr'
    combined_file.flush()    # flush prints everything from the file onto screen
    # command below means print to this file
    print >> combined_file, "XTENSION= 'A3DTABLE'           / extension type\n",'BITPIX  =                    8 / printable ASCII codes\n','NAXIS   =                    2 / Table is a matrix\n','NAXIS1  =                  132 / Max. no. of characters/pass\n','NAXIS2  =             ','%7i' % combined_rows,'/ Number of entries in table\n','PCOUNT  =                    0 / Random parameter count\n','GCOUNT  =                    1 / Group count\n','NOPASS  =                    1 / Number of passes thru table\n','TFIELDS =                    9 / Number of fields in each row\n',"EXTNAME = 'AIPS FG '           / AIPS table file\n",'EXTVER  =                    1 / Version Number of table\n','TBCOL1  =                    9 / Starting char. pos. of field\n',"TFORM1  = 'I11     '           / Fortran format of field  1\n",'TFDIM1  =                    1 / Dimension of field  1\n',"TTYPE1  = 'SOURCE  '           / type (heading) of field  1\n","TUNIT1  = '        '           / physical units of field  1\n",'TBCOL2  =                   20 / Starting char. pos. of field\n',"TFORM2  = 'I11     '           / Fortran format of field  2\n",'TFDIM2  =                    1 / Dimension of field  2\n',"TTYPE2  = 'SUBARRAY                '           / type (heading) of field  2\n","TUNIT2  = '        '           / physical units of field  2\n",'TBCOL3  =                   31 / Starting char. pos. of field\n',"TFORM3  = 'I11     '           / Fortran format of field  3\n",'TFDIM3  =                    1 / Dimension of field  3\n',"TTYPE3  = 'FREQ ID                 '           / type (heading) of field  3\n","TUNIT3  = '        '           / physical units of field  3\n",'TBCOL4  =                   42 / Starting char. pos. of field\n',"TFORM4  = 'I11     '           / Fortran format of field  4\n",'TFDIM4  =                    2 / Dimension of field  4\n',"TTYPE4  = 'ANTS                    '           / type (heading) of field  4\n","TUNIT4  = '        '           / physical units of field  4\n",'TBCOL5  =                   53 / Starting char. pos. of field\n',"TFORM5  = 'E18.8   '           / Fortran format of field  5\n",'TFDIM5  =                    2 / Dimension of field  5\n',"TTYPE5  = 'TIME RANGE              '           / type (heading) of field  5\n","TUNIT5  = 'DAYS    '           / physical units of field  5\n",'TBCOL6  =                   71 / Starting char. pos. of field\n',"TFORM6  = 'I9     '           / Fortran format of field  6\n",'TFDIM6  =                    2 / Dimension of field  6\n',"TTYPE6  = 'IFS                     '           / type (heading) of field  6\n","TUNIT6  = '        '           / physical units of field  6\n",'TBCOL7  =                   80 / Starting char. pos. of field\n',"TFORM7  = 'I9     '           / Fortran format of field  7\n",'TFDIM7  =                    2 / Dimension of field  7\n',"TTYPE7  = 'CHANS                   '           / type (heading) of field  7\n","TUNIT7  = '        '           / physical units of field  7\n",'TBCOL8  =                   89 / Starting char. pos. of field\n',"TFORM8  = 'X4      '           / Fortran format of field  8\n",'TFDIM8  =                    4 / Dimension of field  8\n',"TTYPE8  = 'PFLAGS                  '           / type (heading) of field  8\n","TUNIT8  = '        '           / physical units of field  8\n",'TBCOL9  =                  101 / Starting char. pos. of field\n',"TFORM9  = 'A24     '           / Fortran format of field  9\n",'TFDIM9  =                   26 / Dimension of field  9\n',"TTYPE9  = 'REASON                  '           / type (heading) of field  9\n","TUNIT9  = '        '           / physical units of field  9\n",'END\n','COL. NO.      1          2          3          4            5              6          7          8        9\n','     ROW   SOURCE     SUBARRAY   FREQ ID    ANTS         TIME RAN       IFS        CHANS      PFLAGS   REASON\n','  NUMBER                                                 DAYS\n','***BEGIN*PASS***'
#    dummy = open(path2folder+'dummy.fg', 'wr')
    dummy = open(path2folder+'dummy.fg', 'r')
    for line in dummy:
        line = line.rstrip('\n')
        print >> combined_file, line
    no_rows += 1
    try:
        lovell_dummy = open(path2folder+experiment+'_lovell_dummy.fg', 'r')
        print "Appending Lovell dropouts..."
        even = 0
        #no_rows += 1
        for line in lovell_dummy:
            line = line.rstrip('\n')
            if even == 0:
                num = "%8i" % no_rows
                line = num + line[8:]
                even = 1
                print >> combined_file, line
            else:
                num = "%8i" % no_rows
                line = num + line[8:]
                print >> combined_file, line
                even = 0
                no_rows += 1
    except:
        print "No Lovell dropouts to flag..."
    try:
        zeros_dummy = open(path2folder+experiment+'_zeros_dummy.fg', 'r')
        print "Appending Zero Amplitude dropouts..."
        even = 0
        #no_rows += 1
        for line in zeros_dummy:
            line = line.rstrip('\n')
            if even == 0:
                num = "%8i" % no_rows
                line = num + line[8:]
                even = 1
                print >> combined_file, line
            else:
                num = "%8i" % no_rows
                line = num + line[8:]
                print >> combined_file, line
                even = 0
                no_rows += 1
    except:
        print "No Zero Amplitude dropouts to flag..."
    #print lovellmk2
    #if lovellmk2 != 0:
    for bline in baseline_removal:
        print "Flagging Lovell-Mk2 baseline for all sources..."
        #lovellmk2 = lovellmk2.rsplit('-')
        auto_bline = bline.rsplit('-')
        print >> combined_file, '%8i' % no_rows, '%10i' % 0, '%10i' % 0, '%10i' % 0, '%10s' % auto_bline[0], '   %.8E' % 0, '%8i' % 1, '%8i' % 1, '     \'%s\'' % '1111', '  \'%s\'' % 'SERPENT BASELINE REMOVAL'
        print >> combined_file, '%8i' % no_rows, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % auto_bline[1], '   %.8E' % 9999.99, '%5i' % nif, '%8i' % nchan, '%11s' % '\'\'', '%25s' % '\'\''
        no_rows += 1
    #else:
        #print "Not flagging Lovell-Mk2  baseline..."
    combined_file.close()


# the code below contains a routine to group together flags with consecutive time scans to reduce number of rows in FG table. This could have been done with frequencies and stokes etc, but implementing this in conjunction with the time scans would be tough so as most of the flags are in time just grouping in time will reduce the number of rows in the FG table will be sufficient enough to flag everything without going over 1000000 rows...


try:
    count_limit = flag_coinc_chans
except NameError:
    count_limit = 0


def fg_create_dummy(path2folder, no_rows, tab_row, pflags, flag, times_array, source_num, ifnum, bline):
    # Need number of rows which we're going to enter.
    # Might need to do this more than once for each source. i.e. for each stokes... might need to put this definition at beginning of script
    global fg_row
    fgfile = open(path2folder+'dummy.fg', 'wr')#operation = 'wr'
    fgfile.flush()    # flush prints everything from the file onto screen
    
    source_index = source_num
    reason_def = reason()    # define locally so we only need to call upon the def once
    for col in xrange(flag.shape[1]):
        first = 0
        for row in xrange(flag.shape[0]):
            if flag[row,col] > 0.0 and row != flag.shape[0]-1:
                first += 1
                if flag[row+1,col] == 0.0 and first > count_limit:
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row+1-first])-0.00000116), '%8i' % ifnum, '%8i' % channum, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    flag[row+1-first:row+1, col] = 0.0
                    first = 0
            elif flag[row,col] > 0.0 and row == flag.shape[0]-1:
                if first > count_limit:
                    first += 1
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row+1-first])-0.00000116), '%8i' % ifnum, '%8i' % channum, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    flag[row+1-first:row+1, col] = 0.0
                    first = 0
                else:
                    #print first
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row-first])-0.00000116), '%8i' % ifnum, '%8i' % channum, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    flag[row-first:row+1, col] = 0.0
                    first = 0
    fg_row = tab_row
    fgfile.close()


def fg_append_dummy(path2folder, no_rows, tab_row, pflags, flag, times_array, source_num, ifnum, bline):
    # Need number of rows which we're going to enter.
    # Might need to do this more than once for each source. i.e. for each stokes... might need to put this definition at beginning of script
    global fg_row
    fgfile = open(path2folder+'dummy.fg', 'a')#operation = 'wr'
    fgfile.flush()    # flush prints everything from the file onto screen

    source_index = source_num
    reason_def = reason()    # define locally so we only need to call upon the def once
    for col in xrange(flag.shape[1]):
        first = 0
        for row in xrange(flag.shape[0]):
            if flag[row,col] > 0.0 and row != flag.shape[0]-1:
                first += 1
                if flag[row+1,col] == 0.0 and first > count_limit:
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row+1-first])-0.00000116), '%8i' % ifnum, '%8i' % channum, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    flag[row+1-first:row+1, col] = 0.0
                    first = 0
            elif flag[row,col] > 0.0 and row == flag.shape[0]-1:
                if first > count_limit:
                    first += 1
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row+1-first])-0.00000116), '%8i' % ifnum, '%8i' % channum, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    flag[row+1-first:row+1, col] = 0.0
                    first = 0
                else:
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row-first])-0.00000116), '%8i' % ifnum, '%8i' % channum, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    flag[row-first:row+1, col] = 0.0
                    first = 0
    fg_row = tab_row
    fgfile.close()
    

def fg_create_dummy_freq(path2folder, no_rows, tab_row, pflags, flag, times_array, source_num, ifnum, bline):
    # Need number of rows which we're going to enter.
    # Might need to do this more than once for each source. i.e. for each stokes... might need to put this definition at beginning of script
    global fg_row
    fgfile = open(path2folder+'dummy.fg', 'wr')#operation = 'wr'
    fgfile.flush()    # flush prints everything from the file onto screen

    source_index = source_num
    reason_def = reason()    # define locally so we only need to call upon the def once
    for row in xrange(flag.shape[0]):
        first = 0
        for col in xrange(flag.shape[1]):
            if flag[row,col] > 0.0 and col != flag.shape[1]-1:
                first += 1
                if flag[row,col+1] == 0.0:
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row])-0.00000116), '%8i' % ifnum, '%8i' % (channum-first), '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    first = 0
            elif flag[row,col] > 0.0 and col == flag.shape[1]-1:
                if first != 0:
                    first += 1
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row])-0.00000116), '%8i' % ifnum, '%8i' % (channum-first), '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    first = 0
                else:
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row])-0.00000116), '%8i' % ifnum, '%8i' % channum, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    first = 0
    fg_row = tab_row
    fgfile.close()


def fg_append_dummy_freq(path2folder, no_rows, tab_row, pflags, flag, times_array, source_num, ifnum, bline):
    # Need number of rows which we're going to enter.
    # Might need to do this more than once for each source. i.e. for each stokes... might need to put this definition at beginning of script
    global fg_row
    fgfile = open(path2folder+'dummy.fg', 'a')#operation = 'wr'
    fgfile.flush()    # flush prints everything from the file onto screen

    source_index = source_num
    reason_def = reason()    # define locally so we only need to call upon the def once

    for row in xrange(flag.shape[0]):
        first = 0
        for col in xrange(flag.shape[1]):
            if flag[row,col] > 0.0 and col != flag.shape[1]-1:
                first += 1
                if flag[row,col+1] == 0.0:
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row])-0.00000116), '%8i' % ifnum, '%8i' % (channum+1-first), '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    first = 0
            elif flag[row,col] != 0.0 and col == flag.shape[1]-1:
                if first != 0.0:
                    first += 1
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row])-0.00000116), '%8i' % ifnum, '%8i' % (channum+1-first), '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    first = 0
                else:
                    channum = col+1
                    print >> fgfile, '%8i' % tab_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(times_array[row])-0.00000116), '%8i' % ifnum, '%8i' % (channum), '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                    print >> fgfile, '%8i' % tab_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(times_array[row])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                    tab_row += 1
                    first = 0
    fg_row = tab_row
    fgfile.close()


def lovell_read_create(source_index, pflags, ifnum, channum, reason_def, con, dropouts, row, bline):

    global lovell_row
    lovell_dummy = open(path2folder +experiment+'_lovell_dummy.fg', 'wr')

    lovell_dummy.flush()

    start = dropouts[0]
    for i in xrange(1,len(dropouts)):
        if dropouts[i] - dropouts[i-1] != 1:
            print >> lovell_dummy, '%8i' % row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % 1, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
            print >> lovell_dummy, '%8i' % row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[dropouts[i-1]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
            start = dropouts[i]
            row += 1
            if i == len(dropouts)-1:
                print >> lovell_dummy, '%8i' % row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % 1, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                print >> lovell_dummy, '%8i' % row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[dropouts[i]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                row += 1
        elif i == len(dropouts)-1:
            print >> lovell_dummy, '%8i' % row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % 1, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
            print >> lovell_dummy, '%8i' % row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[dropouts[i]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
            row += 1
    lovell_row = row
    lovell_dummy.close()



def lovell_read_append(source_index, pflags, ifnum, channum, reason_def, con, dropouts, row, bline):

    global lovell_row
    lovell_dummy = open(path2folder +experiment+'_lovell_dummy.fg', 'a')

    lovell_dummy.flush()

    start = dropouts[0]
    for i in xrange(1,len(dropouts)):
        if dropouts[i] - dropouts[i-1] != 1:
            print >> lovell_dummy, '%8i' % row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % 1, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
            print >> lovell_dummy, '%8i' % row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[dropouts[i-1]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
            start = dropouts[i]
            row += 1
            if i == len(dropouts)-1:
                print >> lovell_dummy, '%8i' % row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % 1, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                print >> lovell_dummy, '%8i' % row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[dropouts[i]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                row += 1
        elif i == len(dropouts)-1:
            print >> lovell_dummy, '%8i' % row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % 1, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
            print >> lovell_dummy, '%8i' % row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[dropouts[i]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
            row += 1
    lovell_row = row
    lovell_dummy.close()





def zeros_read_create(source_index, pflags, ifnum, channum, reason_def, con, zeros, z_row, bline):

    global zeros_row
    zeros_dummy = open(path2folder +experiment+'_zeros_dummy.fg', 'wr')

    zeros_dummy.flush()
    start = zeros[0]
    #print con[start]
    #print start
    for i in xrange(1,len(zeros)):
        if zeros[i] - zeros[i-1] != 1:
            print >> zeros_dummy, '%8i' % z_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % ifnum, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
            print >> zeros_dummy, '%8i' % z_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[zeros[i-1]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
            start = zeros[i]
            z_row += 1
            if i == len(zeros)-1:
                print >> zeros_dummy, '%8i' % z_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % ifnum, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                print >> zeros_dummy, '%8i' % z_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[zeros[i]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                z_row += 1
        elif i == len(zeros)-1:
            print >> zeros_dummy, '%8i' % z_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % ifnum, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
            print >> zeros_dummy, '%8i' % z_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[zeros[i]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
            z_row += 1
    zeros_row = z_row
    zeros_dummy.close()



def zeros_read_append(source_index, pflags, ifnum, channum, reason_def, con, zeros, z_row, bline):

    global zeros_row
    zeros_dummy = open(path2folder +experiment+'_zeros_dummy.fg', 'a')

    zeros_dummy.flush()

    start = zeros[0]
    for i in xrange(1,len(zeros)):
        if zeros[i] - zeros[i-1] != 1:
            print >> zeros_dummy, '%8i' % z_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % ifnum, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
            print >> zeros_dummy, '%8i' % z_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[zeros[i-1]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
            start = zeros[i]
            z_row += 1
            if i == len(zeros)-1:
                print >> zeros_dummy, '%8i' % z_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % ifnum, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
                print >> zeros_dummy, '%8i' % z_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[zeros[i]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
                z_row += 1
        elif i == len(zeros)-1:
            print >> zeros_dummy, '%8i' % z_row, '%10i' % source_index, '%10i' % 0, '%10i' % 0, '%10s' % bline[0], '   %.8E' % (float(con[start][0])-0.00000116), '%8i' % ifnum, '%8i' % 1, '     \'%s\'' % pflags, '  \'%s\'' % reason_def
            print >> zeros_dummy, '%8i' % z_row, '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % '\'\'', '%10s' % bline[1], '   %.8E' % (float(con[zeros[i]][0])+0.00000116), '%8i' % ifnum, '%8i' % channum, '%11s' % '\'\'', '%27s' % '\'\''
            z_row += 1
    zeros_row = z_row
    zeros_dummy.close()



print " Creating first FG text file from pickle files..."


if coadd_polarization_flags == 'yes':
    print " Coadding polarizations for FG table..."

    lovell_row = 1
    for source in source_list:
        if nsource > 1:
            name = source_list[source]
        for bline in baselines:
            bline_name = bline
            bline = bline.rsplit('-')
            for j in xrange(nif):
                pname = str(name) + '__IF' + str(j+1) + '__' + str(bline_name) + '__all-stokes'
                print "Reading file:", pname
                if os.path.exists(path2folder+str(pname)+'.p') == True:
                    flag = pickle.load(open(path2folder + str(pname) +".p", "rb"))
                    times_array = pickle.load(open(path2folder + str(pname) +".info", "rb"))
                    source_num = source
                    ifnum = j+1
                    pflags = '1111'
                    if os.path.exists(path2folder+'dummy.fg') == True:
                        fg_append_dummy(path2folder, no_rows, fg_row, pflags, flag, times_array, source_num, ifnum, bline)
                        fg_append_dummy_freq(path2folder, no_rows, fg_row, pflags, flag, times_array, source_num, ifnum, bline)
                    else:
                        fg_create_dummy(path2folder, no_rows, fg_row, pflags, flag, times_array, source_num, ifnum, bline)
                        fg_append_dummy_freq(path2folder, no_rows, fg_row, pflags, flag, times_array, source_num, ifnum, bline)
                    os.remove(path2folder+str(pname) +".p")
                    os.remove(path2folder+str(pname) +".info")
                plname = str(name) + "__" + str(bline_name) + "_IF" + str(j+1) +"__lovell_dummy"
                if os.path.exists(path2folder + str(plname) + ".con") == True:
                    con = pickle.load(open(path2folder + str(plname) + ".con", "rb"))
                    dropouts = pickle.load(open(path2folder + str(plname) + ".dropouts", "rb"))
                    source_num = source
                    pflags = '1111'
                    ifnum = nif
                    channum = nchan
                    reason_def = 'LOVELL DROPOUT ' + strftime("%d-%b-%y")
                    if os.path.exists(path2folder+experiment+'_lovell_dummy.fg') == True:
                        lovell_read_append(source_num, pflags, ifnum, channum, reason_def, con, dropouts, lovell_row, bline)
                    else:
                        lovell_read_create(source_num, pflags, ifnum, channum, reason_def, con, dropouts, lovell_row, bline)
                    os.remove(path2folder + str(plname) + ".con")
                    os.remove(path2folder + str(plname) + ".info")
                    os.remove(path2folder + str(plname) + ".dropouts")
                zname = str(name) + "__" + str(bline_name) + "__IF" + str(j+1) + "__" + pol + "__zeros_dummy"
                #print zname
                if os.path.exists(path2folder + str(zname) + ".con") == True:
                    con = pickle.load(open(path2folder + str(zname) + ".con", "rb"))
                    zeros = pickle.load(open(path2folder + str(zname) + ".zeros", "rb"))
                    source_num = source
                    pflags = '1111'
                    ifnum = j+1
                    channum = nchan
                    reason_def = 'ZEROS DROPOUT ' + strftime("%d-%b-%y")
                    if os.path.exists(path2folder+experiment+'_zeros_dummy.fg') == True:
                        zeros_read_append(source_num, pflags, ifnum, channum, reason_def, con, zeros, zeros_row, bline)
                    else:
                        zeros_read_create(source_num, pflags, ifnum, channum, reason_def, con, zeros, zeros_row, bline)
                    os.remove(path2folder + str(zname) + ".con")
                    os.remove(path2folder + str(zname) + ".info")
                    os.remove(path2folder + str(zname) + ".zeros")


#if coadd_polarization_flags == 'no':

else:
    print " Creating separate polarizations flags for FG table..."

    lovell_row = 1
    for source in source_list:
        if nsource > 1:
            name = source_list[source]
        for bline in baselines:
            bline_name = bline
            bline = bline.rsplit('-')
            for j in xrange(nif):
                for pol in polnames:
                    pname = str(name) + '__IF' + str(j+1) + '__' + str(bline_name) + '__' + pol
                    print "Reading file:", pname
                    if os.path.exists(path2folder+str(pname)+'.p') == True:
                        flag = pickle.load(open(path2folder + str(pname) +".p", "rb"))
                        times_array = pickle.load(open(path2folder + str(pname) +".info", "rb"))
                        source_num = source
                        ifnum = j+1
                        if pol == 'RR' or pol == 'XX':
                            pflags = '1000'
                        if pol == 'LL' or pol == 'YY':
                            pflags = '0100'
                        if pol == 'RL' or pol == 'XY':
                            pflags = '0010'
                        if pol == 'LR' or pol == 'YX':
                            pflags = '0001'
                        #pflags = '1111'
                        if os.path.exists(path2folder+'dummy.fg') == True:
                            fg_append_dummy(path2folder, no_rows, fg_row, pflags, flag, times_array, source_num, ifnum, bline)
                            fg_append_dummy_freq(path2folder, no_rows, fg_row, pflags, flag, times_array, source_num, ifnum, bline)
                        else:
                            fg_create_dummy(path2folder, no_rows, fg_row, pflags, flag, times_array, source_num, ifnum, bline)
                            fg_append_dummy_freq(path2folder, no_rows, fg_row, pflags, flag, times_array, source_num, ifnum, bline)
                        os.remove(path2folder+str(pname) +".p")
                        os.remove(path2folder+str(pname) +".info")
                    # for pol in polnames:
                    zname = str(name) + "__" + str(bline_name) + "__IF" + str(j+1) + "__" + pol + "__zeros_dummy"
                    # print zname
                    if os.path.exists(path2folder + str(zname) + ".con") == True:
                        con = pickle.load(open(path2folder + str(zname) + ".con", "rb"))
                        zeros = pickle.load(open(path2folder + str(zname) + ".zeros", "rb"))
                        source_num = source
                        if coadd_zero_flags == 'no':
                            if pol == 'RR' or pol == 'XX':
                                pflags = '1000'
                            if pol == 'LL' or pol == 'YY':
                                pflags = '0100'
                            if pol == 'RL' or pol == 'XY':
                                pflags = '0010'
                            if pol == 'LR' or pol == 'YX':
                                pflags = '0001'
                        else:
                            pflags = '1111'
                        ifnum = j+1
                        channum = nchan
                        reason_def = 'ZEROS DROPOUT ' + strftime("%d-%b-%y")
                        if os.path.exists(path2folder+experiment+'_zeros_dummy.fg') == True:
                            zeros_read_append(source_num, pflags, ifnum, channum, reason_def, con, zeros, zeros_row, bline)
                        else:
                            zeros_read_create(source_num, pflags, ifnum, channum, reason_def, con, zeros, zeros_row, bline)
                        os.remove(path2folder + str(zname) + ".con")
                        os.remove(path2folder + str(zname) + ".info")
                        os.remove(path2folder + str(zname) + ".zeros")
                plname = str(name) + "__" + str(bline_name) + "_IF" + str(j+1) + "__lovell_dummy"
                if os.path.exists(path2folder + str(plname) + ".con") == True:
                    # print 'True'
                    con = pickle.load(open(path2folder + str(plname) + ".con", "rb"))
                    dropouts = pickle.load(open(path2folder + str(plname) + ".dropouts", "rb"))
                    source_num = source
                    pflags = '1111'
                    ifnum = nif
                    channum = nchan
                    reason_def = 'LOVELL DROPOUT ' + strftime("%d-%b-%y")
                    if os.path.exists(path2folder+experiment+'_lovell_dummy.fg') == True:
                        # print 'yes'
                        lovell_read_append(source_num, pflags, ifnum, channum, reason_def, con, dropouts, lovell_row, bline)
                    else:
                        # print 'no'
                        lovell_read_create(source_num, pflags, ifnum, channum, reason_def, con, dropouts, lovell_row, bline)
                    #os.remove(path2folder + str(plname) + ".con")
                    #os.remove(path2folder + str(plname) + ".info")
                    #os.remove(path2folder + str(plname) + ".dropouts")
                    


no_rows = int(fg_row-1)
print "Number of rows from SumThreshold flagging: %i" % no_rows


# Now print the FG header and append the table data below. Make sure this function is out of the for source loop when including every source in the flagging procedure i.e. making sure there's only one end pass (DON'T INDENT!).
fgheader(path2folder, no_rows)


# Writes the end pass line after appending all the flagged amps are written to the fgfile
def fgend():
    fgfile = open(path2folder+experiment+'.fg', 'a')
    fgfile.flush()
    print >> fgfile, '***END*PASS***'
    fgfile.close()


# Make sure this function is out of the for source loop when including every source in the flagging procedure i.e. making sure there's only one end pass 
fgend()


# Now print the FG header for the Lovell Dropout FG table
lovell_row -= 1    # because counting continues once after final FG row appendage
lovell_header(path2folder, lovell_row)


# Append the FG last line to the file
def lovell_end():
    lovell_dropouts = open(path2folder+experiment+'_lovell_dropouts.fg', 'a')
    lovell_dropouts.flush()
    print >> lovell_dropouts, '***END*PASS***'
    lovell_dropouts.close()


lovell_end()


# Now print the FG header for the Lovell Dropout FG table
zeros_row -= 1    # because counting continues once after final FG row appendage
zeros_header(path2folder, zeros_row)


# Append the FG last line to the file
def zeros_end():
    zeros_dropouts = open(path2folder+experiment+'_zeros_dropouts.fg', 'a')
    zeros_dropouts.flush()
    print >> zeros_dropouts, '***END*PASS***'
    zeros_dropouts.close()


zeros_end()


#print lovell_row, no_rows, zeros_row
#print lovell_row + no_rows + zeros_row + len(baseline_removal)

# Code to combine the flagger rows and Lovell dropout rows:
combined_rows = lovell_row + no_rows + zeros_row + len(baseline_removal)
#if lovellmk2 != 0:
    #combined_rows += 1
#combined_header(path2folder, combined_rows, no_rows, lovellmk2)

#combined_rows += len(baseline_removal)
combined_header(path2folder, combined_rows, no_rows, baseline_removal, nif, nchan)


def combined_end():
    header = open(path2folder + experiment + "_combined.fg", 'a')
    header.flush()
    print >> header, '***END*PASS***'
    header.close()

combined_end()


print "\n Number of rows for FG table: %i" % combined_rows


# Now remove the dummy file. After all... it's done it's job and was only a pawn. Pawn's can be promoted, but can never become players...
if os.path.exists(path2folder+'dummy.fg') == True:
    os.remove(path2folder+'dummy.fg')


if os.path.exists(path2folder+experiment+'_lovell_dummy.fg') == True:
    os.remove(path2folder+experiment+'_lovell_dummy.fg')

if os.path.exists(path2folder+experiment+'_zeros_dummy.fg') == True:
    os.remove(path2folder+experiment+'_zeros_dummy.fg')



# Definition to write text file to FG table.
def tbin():
    tbin = AIPSTask('TBIN')
    tbin.outname = experiment
    tbin.outclass = experiment_klass
    tbin.outseq = experiment_seq
    tbin.outdisk = experiment_disk
    tbin.intext = path2folder + experiment + '_combined.fg'     # text file name
    tbin.bcount = 0
    tbin.ecount = 0
    tbin.inp()
    tbin.go()


print " Reading flag txt file into AIPS via TBIN..."
tbin()


# On the new version of AIPS a verb called REFLG is available which concatenates rows in the FG table more efficiently. This will be very useful when flagging L band data as the RFI is no longer '1 dimensional' but 2D, but it also has a function where you can flag the entire time scan or channel if x% of the time scan or channel is flagged... Very useful as some RFI in a channel (e.g.) might be too weak to be picked up, but most of the rest of the channel is flagged... Makes the flagger more robust and efficient with its rows in the FG table...

def reflg(cparm_1, cparm_2, cparm_3, cparm_4, cparm_5, cparm_6, cparm_7):
    reflg = AIPSTask('REFLG')
    # changed from uvdata.name etc...
    reflg.inname = experiment
    reflg.inclass = experiment_klass
    reflg.inseq = experiment_seq
    reflg.indisk = experiment_disk
    reflg.flagver = 0
    #reflg.cparm[1:] = [0, 2, 0.70, 0.99, 0.70, 0.70, 2]
    if flagging_options == 'default':
        reflg.cparm[1:] = [0, 2, 0.70, 0.99, 0.70, 0.70, 2]
    else:
        reflg.cparm[1:] = [cparm_1, cparm_2, cparm_3, cparm_4, cparm_5, cparm_6, cparm_7]
    # These control the flagging conditions:
    # 1: interval
    # 2: flags channels if between cparm[2] flagged channels
    # 3: fraction of channels flagged at specific time - flag all channels
    # 4: fraction of times flagged at specific channel - flag all times
    # 5: fraction of baselines flagged at specific channel and time - flag all bline
    # 6: fraction of baselines flagged at specific antenna - flag all bline with ant
    # 7: cparm[7] = 1, flag cross hand pols, cparm[7] = 2, flag all pols
    reflg.inp()
    reflg.go()


# Write out the REFLG'd FG table to a text file...
def tbout():
    tbout = AIPSTask('TBOUT')
    tbout.inname = experiment
    tbout.inclass = experiment_klass
    tbout.inseq = experiment_seq
    tbout.indisk = experiment_disk
    tbout.inext = 'FG'
    tbout.invers = 0
    tbout.bcount = 0
    tbout.ecount = 0
    tbout.outtext = path2folder + experiment + "_r.fg"    # reflg'd file name
    tbout.inp()
    tbout.go()


# Try and use REFLG, if AIPS version does not contain verb then skip over it.

try:
    reflg(cparm_1, cparm_2, cparm_3, cparm_4, cparm_5, cparm_6, cparm_7)
    tbout()
except:
    print "\n Current AIPS Version does not have REFLG verb. Need a newer version."
    print " Passing over REFLG function..."




print "Time taken (hh:mm:ss):", time2hms(time.time()-ti)
print 'Finished on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),'\n'



if write2log == 'yes':
    try:
        log = open(path2folder + experiment + "_run_log", 'a')
        log.flush()
        t = time2hms(time.time()-ti)
        d = strftime("%d-%b-%y")
        dt = strftime("%H:%M:%S", localtime())
        print >> log, "NCPUs:", NCPU, ", Experiment Name:", experiment, ", Run #", runs, ", Time taken (hh:mm:ss):", time2hms(time.time()-ti), ", Date and Time: %s" % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()), "SERPent Version:", version, '\n'
        log.close()
    except:
        pass
