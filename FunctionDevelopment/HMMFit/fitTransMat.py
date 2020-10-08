#!/usr/bin/python

import sys
#import code.util as util
import util
import numpy as np
from scipy import io as scio
from pomegranate import *
import pickle
import os

def main():
    # Get command line arguments
    basepath_ind = int(sys.argv[1])-1
    UID = int(sys.argv[2])-1

    # Retrieve basepath
    base_dict = scio.loadmat('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/FunctionDevelopment/HMMFit/basepaths_test.mat')
    basepath = base_dict['basepaths'][basepath_ind,0][0]
    #basepath_original = base_dict['basepaths'][basepath_ind,0][0]
    #base_split = basepath_original.split(os.path.sep)
    #basepath = os.path.join('/gpfs/scratch/rh2618/rh_data', base_split[-1-2], base_split[-1-1], base_split[-1])

    # Load in data
    GammaFits = util.getGammaFits(basepath)
    spikes = util.getSpikes(basepath)
    SleepState = util.getStates(basepath)

    # Load spike train
    spk = spikes['times'][0,UID].flatten()

    brainstates = ['NREMstate', 'WAKEstate']
    brainstates1 = ['NREM', 'WAKE']

    out = {}
    outfile = os.path.join( basepath, 'GammaProcessed', str(UID+1)+'.mat' )

    # Loop over NREM / WAKE states
    for state in range(len(brainstates)):

        # Pull out UIDs
        UID_state = GammaFits[brainstates[state]][0,0]['cellstats'][0,0]['UID'].flatten() - 1

        # If this UID has been fit
        index = np.argwhere( UID_state == UID )
        if index.size:
            index = index[0][0]
            print('Fitting '+brainstates1[state])
            ## Get state specific ISIs and state specific GammaFit params specifi

            spk_stateISI, spk_stateT = util.getStateDepISIs(SleepState, spk)

            # Parameters for ground state (neuron specific)
            gscvs = GammaFits[brainstates[state]]['sharedfit'][0,0]['GSCVs'][0,0]
            gsrates = GammaFits[brainstates[state]]['sharedfit'][0,0]['GSlogrates'][0,0]
            lambda_gs, k_gs = util.toAlphaBeta(gsrates,gscvs)

            # Parameters for activated state (common across neurons)
            ascvs = GammaFits[brainstates[state]]['sharedfit'][0,0]['ASCVs'][0,0]
            asrates = GammaFits[brainstates[state]]['sharedfit'][0,0]['ASlogrates'][0,0]
            lambda_as, k_as = util.toAlphaBeta(asrates, ascvs)

            ## Construct model - freeze parameters of the emission distribution

            # Meaningful state names
            sep = " "
            state_names = [sep.join(("Activated state",str(x+1))) for x in range(lambda_as.size) ]
            state_names.insert(0, "Ground state")

            # Initialize gamma distributions
            dists = [GammaDistribution(k_as[x], lambda_as[x]) for x in range(lambda_as.size)]
            dists.insert(0, GammaDistribution(k_gs[index], lambda_gs[index]))

            # Random (and dense) transition matrix
            trans_mat = np.random.uniform(0,1, (len(dists),len(dists)))
            for kp in range(len(dists)):
                trans_mat[kp,:] = np.divide( trans_mat[kp,:], np.sum( trans_mat[kp,:] ) )

            starts = np.tile(1/len(dists), len(dists))
            ends = np.tile(0, len(dists))

            # Initialize and freeze emission distribution parameters
            model = HiddenMarkovModel.from_matrix(trans_mat, dists, starts, ends, state_names)
            model.freeze_distributions()
            # Fit the transition matrix
            model.fit( sequences=spk_stateISI[brainstates1[state]], algorithm='baum-welch')

            # Extract state parameters in order of storage
            k_all = [ model.states[x].distribution.parameters[0] for x in range(len( model.states )-2) ]
            lambda_all = [ model.states[x].distribution.parameters[1] for x in range(len( model.states )-2) ]
            state_names = [ model.states[x].name for x in range(len( model.states )-2) ]
            logrates, cvs = util.toRateCV(lambda_all, k_all)

            # Predict
            seq_len = len( spk_stateISI[brainstates1[state]] )
            decoded_mode = np.empty(seq_len,dtype=object)
            state_isi = np.empty(seq_len,dtype=object)
            state_spk = np.empty(seq_len,dtype=object)
            for seq_index in range( seq_len ):
                decoded_mode[seq_index] = np.array( model.predict(spk_stateISI[brainstates1[state]][seq_index], algorithm='viterbi')[1:] )+1
                state_isi[seq_index] = spk_stateISI[brainstates1[state]][seq_index]
                state_spk[seq_index] = spk_stateT[brainstates1[state]][seq_index]

            out[brainstates1[state]] = {'UID':UID+1,'basepath':basepath, 'decoded_mode':decoded_mode, 'state_isi':state_isi, 'state_label':state_names, 'state_spk':state_spk, 'logrates':logrates, 'cvs':cvs, 'trans_mat':model.dense_transition_matrix()[:logrates.size+1,:logrates.size]}

    scio.savemat(outfile, out)

if __name__ == "__main__":
    main()
