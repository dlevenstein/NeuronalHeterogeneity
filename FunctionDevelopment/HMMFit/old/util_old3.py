from scipy import io as scio
import numpy as np
import os
from scipy.interpolate import pchip_interpolate
from itertools import combinations
from scipy.optimize import curve_fit

# Utility functions

# Return sleep states as dictionary
def getStates(basepath):

    # Check for file existence
    basename = os.path.basename(basepath)
    filepath = os.path.join( basepath, basename+'.SleepState.states.mat' )
    if not os.path.isfile( filepath ):
        return -1

    states = scio.loadmat(filepath)
    # Extract state names out of the struct
    stateNames = []
    for kp in range( states['SleepState'][0][0][1][0][0][0][0].size ):
        if states['SleepState'][0][0][1][0][0][0][0][kp].size > 0:
            stateNames.append( states['SleepState'][0][0][1][0][0][0][0][kp][0] )

    # Extract state scoring as a dict
    ntup = len( states['SleepState'][0][0][0][0][0] )
    SleepState = {}
    for kp in range(ntup):
        SleepState[stateNames[kp]] = states['SleepState'][0][0][0][0][0][kp]

    return SleepState

# Return kilosort file, if exists
def getKilosortFolder(basepath):
    kilosort_file = ''
    for file in os.listdir(basepath):
        if file.startswith('Kilo'):
            kilosort_file = file
    return kilosort_file

# Return spikes
def getSpikes(basepath, loadKilo=False):
    basename = os.path.basename(basepath)
    # kilosort_file = getKilosortFolder(basepath)
    # if kilosort_file and loadKilo:
    #     spikespath = os.path.join(basepath, kilosort_file, basename+'.spikes.cellinfo.mat')
    # else:
    spikespath = os.path.join(basepath, basename+'.spikes.cellinfo.mat')
    out = scio.loadmat(spikespath)
    return out['spikes'][0,0]

# Return spikes
def getGammaFits(basepath):
    basename = os.path.basename(basepath)
    gammapath = os.path.join(basepath, basename+'.GammaFit.cellinfo.mat')
    out = scio.loadmat(gammapath)
    return out['GammaFit'][0,0]


def getStateDepISIs(states, spk):
    # Extract neuron's state dependent ISI sequences
    dictkey = list( states.keys() )

    n_isiSeq = {}
    spkSeq = {}

    # Loop over states
    for kp in range( len(dictkey) ):
        ints = states[dictkey[kp]]
        # Loop over epochs
        isiseq = []
        spkseq = []
        for jp in range( ints.shape[0] ):
            args = np.argwhere( (spk >= ints[jp,0]) & (spk <= ints[jp,1]) ).flatten()
            if np.any( args == 0):
                args = np.delete(args,np.argwhere( args == 0) )
            if args.size > 0:
                isiseq.append( np.diff( spk[ np.arange(args[0]-1, args[-1]+1) ] ) )
                spkseq.append( spk[ np.arange(args[0], args[-1]+1) ] )
        n_isiSeq[dictkey[kp]] = isiseq
        spkSeq[dictkey[kp]] = spkseq

    return n_isiSeq, spkSeq

# Convert logrates / cvs into  rate / shape (input to pomegranate)
def toAlphaBeta(lograte, cv):
    # NOTE: this anticipates arrays

    # Flatten rates, make sure it's numpy
    if isinstance(lograte, list):
        lograte = np.array(lograte).flatten()
    else:
        lograte = lograte.flatten()

    # Flatten cvs, make sure it's numpy
    if isinstance(cv, list):
        cv = np.array(cv).flatten()
    else:
        cv = cv.flatten()

    # Get the shape (alpha == k)
    k = np.divide( 1, cv )
    # Get the rate (beta == lambda)
    lambd = np.multiply( np.power( 10, lograte), k )

    return lambd, k

# Take rate / shape back to lograte / CV (interpretable as the rate / CV)
def toRateCV(lambd, k):
    # NOTE: this anticipates arrays

    # Flatten lambdas, make sure it's numpy
    if isinstance(lambd, list):
        lambd = np.array(lambd).flatten()
    else:
        lambd = lambd.flatten()

    # Flatten k, make sure it's numpy
    if isinstance(k, list):
        k = np.array(k).flatten()
    else:
        k = k.flatten()

    # Get CV
    cv = np.divide( 1, k )
    # Get rate
    lograte = np.log10( np.divide(lambd, k) )

    return lograte, cv

def get2folds( isi_seq ):

    # this can't really be changed.. 
    prop = 0.5

    nisis = np.array( [ isi_seq[x].size for x in range(len( isi_seq)) ] )
    nisis_tot = np.sum( nisis )

    sort_inds = np.flip( np.argsort(nisis, axis=None) )

    vtot = 0
    fold1 = []
    fold2 = []
    for kp in range(sort_inds.size):
        vtot += nisis[ sort_inds[kp] ]
        if vtot / nisis_tot > prop:
            vtot -= nisis[ sort_inds[kp] ]
            fold1.append( sort_inds[kp] )
        else:
            fold2.append( sort_inds[kp] )

    fold1_isis = [ isi_seq[x] for x in fold1 ]
    fold2_isis = [ isi_seq[x] for x in fold2 ]
    return (fold1_isis, fold2_isis)

def getOptNStates( x_orig, deviance ):

    # Describe likelihood curve with a model
    x = np.linspace(x_orig[0], x_orig[-1], num=1000)
    y = pchip_interpolate(x_orig, deviance, x)

    # Estimate derivatives by numerical differencing
    der1 = np.diff(y)
    der2 = np.diff(der1)
    der1 = der1[:-1]

    est_curvature = abs(der2)*(1 + der1**2)**-1.5
    max_curvature = np.argmax( est_curvature )

    # Number of states to return is that with greatest curvature
    return x_orig[ np.argmin( np.abs( x[ max_curvature ] - x_orig ) ) ]

# Describe likelihood curve with a model - exponential
def getOptNStates_v1( x_orig, ll ):
    
    x = np.linspace( x_orig[0], x_orig[-1], num=100)

    # Define exponential function
    def func(t, a, b, alpha):
        return a - b * np.exp(-alpha * t)

    # Initial parameters of exponential 
    a0 = ll[-1]
    b0 = ll[0]
    alpha0 = 1/x_orig[-1]

    # Coefficients and curve fit for curve
    popt4, pcov4 = curve_fit(func, x_orig, ll, p0=(a0, b0, alpha0))

    a, b, alpha = popt4
    y_hat = func(x, a, b, alpha)

    der1 = np.diff(y_hat)
    der2 = np.diff(der1)
    der1 = der1[:-1]

    est_curvature = abs(der2)*(1 + der1**2)**-1.5
    max_curvature = np.argmax( est_curvature )

    params = {}
    params['a'] = a
    params['b'] = b
    params['alpha'] = alpha
    # Number of states to return is that with greatest curvature
    optNstates = x_orig[ np.argmin( np.abs( x[ max_curvature ] - x_orig ) ) ]
    return optNstates, params

def getCombos(nums, rang):
    all_tup = []
    start = rang[0]
    stop = rang[1]+1
    for i in range(start,stop):
        all_tup.extend( list(combinations(nums, i)) )
    all_tup = np.array( all_tup )
    # Get all tuples with no 0 at the start
    to_delete = np.argwhere( [ all_tup[k][0] != 0 for k in range(len(all_tup)) ]  ).flatten()
    return np.delete(all_tup, to_delete)

# Split list into desired number of sublists
def get_sublists(original_list, number_of_sub_list_wanted):
    sublists = list()
    for sub_list_count in range(number_of_sub_list_wanted): 
        sublists.append(np.array( original_list[sub_list_count::number_of_sub_list_wanted] ) )
    return sublists

