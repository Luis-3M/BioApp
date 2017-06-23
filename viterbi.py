from hmm import HMM
import numpy as np
import time

#Viterbi algorithm
def viterbi(hmm, initial_dist, emissions):
    probs = hmm.emission_dist(emissions[0]) * initial_dist
    stack = []
    print
    print "\t\t-> Probability Matrix"
    print
    for emission in emissions[1:]:
        trans_probs = hmm.transition_probs * np.row_stack(probs)
        max_col_ixs = np.argmax(trans_probs, axis=0)
        probs = hmm.emission_dist(emission) * trans_probs[max_col_ixs, np.arange(hmm.num_states)]
        for item in probs:
            print str(item)+" | \t",
        print
        print
        stack.append(max_col_ixs)
    state_seq = [np.argmax(probs)]
    while stack:
        max_col_ixs = stack.pop()
        state_seq.append(max_col_ixs[state_seq[-1]])
    state_seq.reverse()
    return state_seq

def encode(seq):
    emissions = []
    for ch in seq.upper():
        if ch == "A":
            emissions.append(0)
        elif ch == "C":
            emissions.append(1)
        elif ch == "G":
            emissions.append(2)
        elif ch == "T":
            emissions.append(3)
    return emissions

def main():
    start = time.clock()
    transition_probs = np.array([
                                [0.0, 1.0, 0.0, 0.0],  #0 = Begin
                                [0.0, 0.9, 0.1, 0.0],  #1 = Exon
                                [0.0, 0.0, 0.0, 1.0],  #2 = Donor
                                [0.0, 0.0, 0.0, 1.0]   #3 = Intron
                                ])

    #                             A    C    G    T
    emission_probs = np.array([
                                [0.0, 0.0, 0.0, 0.0],       #0 = Begin
                                [0.25, 0.25, 0.25, 0.25],   #1 = Exon
                                [0.05, 0.0, 0.95, 0.0],     #2 = Donor
                                [0.4, 0.1, 0.1, 0.4]        #3 = Intron
                                ]) 

    initial_dist = np.array([[0.6, 0.4, 0.5, 0.3]])
    seq = raw_input("Please insert a sequence: ")
    emissions = encode(seq)
    hmm = HMM(transition_probs, emission_probs)
    stateSeq = viterbi(hmm, initial_dist, emissions)
    print
    print "Viterbi Best Path: ",
    print " -> ".join(["S"+str(x) for x in stateSeq])
    print
    end = time.clock()
    print "Program Running Time: "+str(abs(start-end))+" seconds"
    print
    return