from scipy import io as scio
import numpy as np
import os

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
def getSpikes(basepath):
    basename = os.path.basename(basepath)
    kilosort_file = getKilosortFolder(basepath)
    if not kilosort_file:
        spikespath = os.path.join(basepath, basename+'.spikes.cellinfo.mat')
    else:
        spikespath = os.path.join(basepath, kilosort_file, basename+'.spikes.cellinfo.mat')
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


        
