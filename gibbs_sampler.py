#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random 
from collections import Counter
#%% Generates a random k-mer motif for each DNA string
def motif_selector(lines, k):
    motifs = []
    
    for line in lines:
        rand = random.randrange(0,500-k)
        motif = line[rand:rand+k]
        motifs.append(motif)
        
    return motifs
#%% Generates count matrix for each base counts in each column
def count_matrix_generator(motifs, k, del_motif):
    motif_matrix = motifs.copy()
    # remove the selected motif from the matrix
    del motif_matrix[del_motif]
    column_matrix = []
    # we assume that all the bases should have their counts at least 1 in each column
    count_matrix = {'A':[1]*k , 'C':[1]*k , 'G':[1]*k , 'T':[1]*k}
    
    # create a column list from the bases in the same location at the motifs 
    for j in range(k):
        column = []
        for motif in motif_matrix:
            column.append(motif[j])
        # count every base at that column and add the numbers to the corresponding column of count matrix
        count_column = Counter(column)
        column_matrix.append(count_column)
        count_matrix['A'][j] += column_matrix[j]['A']
        count_matrix['C'][j] += column_matrix[j]['C']
        count_matrix['G'][j] += column_matrix[j]['G']
        count_matrix['T'][j] += column_matrix[j]['T']
        
    return count_matrix
#%% Generates a profile matrix with the probabilities of the bases in each column
def profile_matrix_generator(motifs, k, del_motif):
    count_matrix = count_matrix_generator(motifs, k, del_motif)
    profile_matrix = {'A':[0]*k , 'C':[0]*k , 'G':[0]*k , 'T':[0]*k}
    col_sum = count_matrix['A'][0] + count_matrix['C'][0] + count_matrix['G'][0] + count_matrix['T'][0] 
    
    for j in range(k):
        profile_matrix['A'][j] = count_matrix['A'][j] / col_sum
        profile_matrix['C'][j] = count_matrix['C'][j] / col_sum
        profile_matrix['G'][j] = count_matrix['G'][j] / col_sum
        profile_matrix['T'][j] = count_matrix['T'][j] / col_sum
        
    return profile_matrix
#%% A weighted random generator
    # Selects a k-mer motif profile-randomly according to every line's bias values
def profile_random_generator(profile, dna_line, k):
    prob_list = [0]*(500-k)
    probs_sum = 0
    
    # compute every possible k-mer in selected line
    for x in range(500-k):
        probability = 1
        k_mer = dna_line[x:x+k]
        # compute the probability of that k-mer based on the values at profile matrix
        for y in range(k):
            probability *= profile[k_mer[y]][y]
        prob_list[x] = probability
        probs_sum += probability
    
    # select a profile-random k-mer based on the bias values
    rand = random.uniform(0, probs_sum)
    num = 0
    for x in range(500-k):
        bias = prob_list[x]
        num += bias
        if rand <= num:
            return dna_line[x:x+k]
#%% Finds the consensus string for a list of motifs
def find_consensus(motifs, k):
    consensus = []
    # create a column list from the bases in the same location at the motifs 
    for y in range(k):
        column = []
        for x in range(1, 10):
            column.append(motifs[x][y])
        # get the most common base in that column and append it to consensus string
        consensus.append(Counter(column).most_common(1)[0][0])
    return consensus
#%% Computes the score for a given list of motifs
def score(motifs, consensus, k):
    score = 0
    for y in range(k):
        column_score = 0
        for x in range(10):
            if not consensus[y] is motifs[x][y]:
                column_score += 1
        score += column_score
    return score
#%% Generates the list of k-mer motifs for each DNA string, 
    # finds the consensus string, and computes the score for the list of motifs
def gibbsSampler(dna_lines, k):
    # generate motifs
    motifs = motif_selector(dna_lines, k)
    # initialize the best_motifs as the first motifs list
    best_motifs = motifs.copy()
    final_consensus = ""
    
    # iterate the algorithm 1000 times
    for j in range(1000):
        # select a random motif for exception
        excepted_motif = random.randrange(0,10)
        line = lines[excepted_motif]
        # generate the profile matrix for given list of motifs except the selected motif
        profile = profile_matrix_generator(motifs, k, excepted_motif)
        # generate a new motif profile-randomly
        new_motif = profile_random_generator(profile, line, k)
        # replace the excepted motif with the new motif 
        motifs[excepted_motif] = new_motif
        # find the consensus for the new motifs list
        consensus = find_consensus(motifs, k)
        # find the consensus for the best motifs list
        best_consensus = find_consensus(best_motifs, k)
        
        # compare the scores of new motifs list and best motifs list
        if score(motifs, consensus, k) < score(best_motifs, best_consensus, k):
            # if the new motifs list gives a better score, make it the best motifs list
            best_motifs = motifs.copy()
            
    # return the consensus as a string
    for base in best_consensus:
        final_consensus = final_consensus + base
        
    return best_motifs, score(best_motifs, best_consensus, k), final_consensus
#%% 
file = open('input.txt', 'r')
lines = file.readlines()

##### 9-mer Motifs #####
nine_mer_motifs, nine_score, nine_consensus = gibbsSampler(lines, 9)
print("9-mer motifs:")
for motif in nine_mer_motifs:
    print(motif)
print("Consensus: {}\nScore: {}\n".format(nine_consensus, nine_score))

##### 10-mer Motifs #####
ten_mer_motifs, ten_score, ten_consensus = gibbsSampler(lines, 10)
print("10-mer motifs:")
for motif in ten_mer_motifs:
    print(motif)
print("Consensus: {}\nScore: {}\n".format(ten_consensus, ten_score))

##### 11-mer Motifs #####
eleven_mer_motifs, eleven_score, eleven_consensus= gibbsSampler(lines, 11)
print("11-mer motifs:")
for motif in eleven_mer_motifs:
    print(motif)
print("Consensus: {}\nScore: {}\n".format(eleven_consensus, eleven_score))
