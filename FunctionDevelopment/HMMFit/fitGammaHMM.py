#!/usr/bin/python

import sys
#import code.util as util
import util
import numpy as np
from scipy import io as scio
from pomegranate import *
import pickle
import os
from sklearn.model_selection import KFold

def main():
    # Get command line arguments
    basepath_ind = int(sys.argv[1])-1
    UID = int(sys.argv[2])

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

    spk_index = util.UIDtoIndex(spikes['UID'].flatten(), UID)


    # Load spike train
    spk = spikes['times'][0,spk_index].flatten()


    brainstates = ['NREMstate', 'WAKEstate']
    brainstates1 = ['NREM', 'WAKE']

    # Number of restarts from a randomly initialized transition matrix
    nrestarts = 10
    nfolds = 6
    nstates = 8
    models = []
    lls = []
    dirname = 'GammaProcessed2'

    out = {}
    outfile = os.path.join( basepath, dirname, str(UID)+'.mat' )
    if not os.path.isdir( os.path.join( basepath, dirname)  ):
        os.mkdir( os.path.join( basepath, dirname) )

    # Loop over NREM / WAKE states
    for state in range(len(brainstates)):

        # Pull out UIDs of units that have been considered in mode decomposition up till this point
        try:
            UID_state = GammaFits[brainstates[state]][0,0]['cellstats'][0,0]['UID'].flatten()
        except:
            continue

        # If this UID has been fit
        gamma_index = util.UIDtoIndex(UID_state, UID)
        if gamma_index is not None:
            print('Fitting '+brainstates1[state])
            ## Get state specific ISIs and state specific GammaFit params specifi

            spk_stateISI, spk_stateT = util.getStateDepISIs(SleepState, spk)

            # Number of sequences in this brain state
            seq_len = len( spk_stateISI[brainstates1[state]] )


            # # Meaningful state names
            # sep = " "
            # state_names = [sep.join(("Activated state",str(x+1))) for x in range(lambda_as.size) ]
            # state_names.insert(0, "Ground state")

            # Cross validate to find the number of states
            isi_flat = np.hstack( spk_stateISI[brainstates1[state]] )
            kf = KFold(n_splits=nfolds)
            kf.get_n_splits(isi_flat)

            all_lls = []
            for nstate in range(1,nstates+1):
                print(nstate)
                test_lls = []
                fold = 1
                for train_index, test_index in kf.split(isi_flat):
                    isi_train, isi_test = isi_flat[train_index], isi_flat[test_index]

                    # Fit HMM with N states on train fold
                    counter = 1
                    while True:
                        try:
                            model = HiddenMarkovModel.from_samples(GammaDistribution, n_components=nstate, X=util.get_sublists(isi_train.tolist(), nfolds-1), algorithm='baum-welch')
                            print('Took', counter, 'attempts to fit fold',fold, 'of Nstate =', nstate)
                            break
                        except:
                            counter+=1
                    # Compute the log prob of sequences in test fold
                    test_lls.append( model.log_probability( isi_test )  )
                    fold+=1
                all_lls.append( test_lls )

            # Return optimal number of states - point of greatest curvature in the model deviance x Nstates curve
            all_lls = np.array([numpy.array(xi) for xi in all_lls])    # list of lists to numpy array
            ll = np.mean( all_lls, axis=1)

            optNstates, exp_params = util.getOptNStates_v1(range(1,nstates+1), ll)

            # Fit model multiple times
            for ip in range(nrestarts):

                # Fit from a random initialization
                counter = 1
                while True:
                    try:
                        model = HiddenMarkovModel.from_samples(GammaDistribution, n_components=optNstates, X=spk_stateISI[brainstates1[state]], algorithm='baum-welch')
                        print('Took', counter, 'attempts to fit',ip, '-th random initialization')
                        break
                    except:
                        counter+=1
                # Store the models and their associated log likelihoods
                models.append( model.copy() )
                lls.append( sum( [ model.log_probability(spk_stateISI[brainstates1[state]][k] ) for k in range( seq_len ) ] ) )

            model = models[ np.argmax(lls) ]

            # Extract state parameters in order of storage
            k_all = [ model.states[x].distribution.parameters[0] for x in range(len( model.states )-2) ]
            lambda_all = [ model.states[x].distribution.parameters[1] for x in range(len( model.states )-2) ]
            # state_names = [ model.states[x].name for x in range(len( model.states )-2) ]
            logrates, cvs = util.toRateCV(lambda_all, k_all)

            # Predict
            # seq_len = len( spk_stateISI[brainstates1[state]] )
            decoded_mode = np.empty(seq_len,dtype=object)
            prob_mode = np.empty(seq_len,dtype=object)
            state_isi = np.empty(seq_len,dtype=object)
            state_spk = np.empty(seq_len,dtype=object)
            for seq_index in range( seq_len ):
                decoded_mode[seq_index] = np.array( model.predict(spk_stateISI[brainstates1[state]][seq_index], algorithm='viterbi')[1:] )+1
                prob_mode[seq_index] = model.predict_proba( spk_stateISI[brainstates1[state]][seq_index] )
                state_isi[seq_index] = spk_stateISI[brainstates1[state]][seq_index]
                state_spk[seq_index] = spk_stateT[brainstates1[state]][seq_index]

            out[brainstates1[state]] = {'UID':UID,'basepath':basepath, 'decoded_mode':decoded_mode, 'prob_mode':prob_mode, 'state_isi':state_isi, \
            'state_spk':state_spk, 'logrates':logrates, 'cvs':cvs, 'trans_mat':model.dense_transition_matrix()[:logrates.size+1,:logrates.size], \
            'all_lls':all_lls, 'optNstates':optNstates, 'exp_params':exp_params}

    scio.savemat(outfile, out)

if __name__ == "__main__":
    main()
