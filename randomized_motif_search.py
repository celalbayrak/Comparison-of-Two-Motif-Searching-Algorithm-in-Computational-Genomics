# -*- coding: utf-8 -*-
"""
Created on Thu May 28 10:38:36 2020

@author: Celal
"""
import random
import collections
nucleotids=["A","T","G","C"]
#generates line of text, by using mutated k-mer. It inserts k-mer to random position of line.
def line_generator(mutated_k_mer):
    global nucleotids
    line=""
    for i in range(0,490):
        rand=random.randrange(4)
        line=line+nucleotids[rand]
    rand2=random.randrange(0,490)
    line=line[0:rand2]+mutated_k_mer+line[rand2:]
    return line

#generates input.txt file by using lines, returned from line_generator() function. 
def txt_generator(line_num,mutated_k_mers):
    lines=[]
    f = open("input.txt", "w")
    for i in range(0,line_num):
        line=line_generator(mutated_k_mers[i])
        lines.append(line)
    for l in lines:
        f.write("%s\n" %l)

#generates k-mer randomly.
def k_mer_generator(k):
    global nucleotids
    k_mer=""
    for i in range(0,k):
        rand=random.randrange(4)
        k_mer=k_mer+nucleotids[rand]
    return k_mer

#mutates k-mer's 4 nucleotides randomly. To make sure 4 mutations occur, it checks whether
#previous nucleotid equals to mutated nucleotid or not. If they are equal, it changes the mutated
#nucleotide.
def k_mer_mutator(k_mer):
    global nucleotids
    mutateds=[]
    for i in range(0,10):
        rands1=[]
        while len(rands1)!=4:
            rand1=random.randrange(10)
            if rand1 not in rands1:
                rands1.append(rand1)
    for i in range(0,10):
        temp = k_mer
        for rand in rands1:
            rand2=random.randrange(4)
            temp_arr=list(temp)
            if temp_arr[rand] != nucleotids[rand2]: #checking process which mentioned above.
                temp_arr[rand]=nucleotids[rand2]
            elif rand2 == 3:
                temp_arr[rand]=nucleotids[rand2 - 1] # changing the mutated nucleotide.
            else:
                temp_arr[rand]=nucleotids[rand2 +1]  # changing the mutated nucleotide.
            temp="".join(temp_arr)
        mutateds.append(temp)
    return mutateds

#randomized motif search algorithm. Inputs are lines of .txt file and k of k-mer.
def randomized_motif_search(lines,k):
    motifs=[]
    random_indexes=[]
    #randomly generating the initial 10 motifs. 1 motif from each line.
    for i in range(0,10):
        rand=random.randrange(491)
        random_indexes.append(rand)
        motif=lines[i][rand:rand+k]
        motifs.append(motif)
    prev_score=1000     #initial unreachable score.
    while True:         #infinite loop. if score does not improve in any iteration, loop ends.
        counters=[]
        score=0
        for i in range(0,k):    #these nested loops take the columns of motif matrix.
            arr=[]
            for motif in motifs:
                arr.append(motif[i]) # append the columns to a list.
            counter=collections.Counter(arr) # for each column counter object counts each element of column 
                                            #and keeps those variables in dictionary format.
            counters.append(counter)    # for each column there is 1 counter object.
            score=score+(10-counter.most_common(1)[0][1]) #score calculation. The score is summation of (element number of each column - occurence number of most common element of each column)
        if score < prev_score: # if score improves continue the loop.
            prev_score=score   # keeps the previous score.
            motifs=[]
            for line in lines:  #these nested loops for probability calculation of each k-mer in .txt file.
                max_prob=0      #initial maximum probability of each line
                best_word=""    #initial best k-mer of each line
                for index in range(0,500-k): 
                    word=line[index:index+k] #each k-mer of line
                    probability=1           #initial probability is 1. Because 1 is ineffective element for multiplication 
                    for j in range(0,k):
                        letter=word[j]      #each letter of k-mer
                        frequency=float(counters[j][letter])/float(k) #frequency of each letter of k-mer in profile matrix.
                        probability=frequency*probability        #probability is multiplication of each letters frequency
                    if probability > max_prob:        #for each line max_prob keeps the highest probability
                        max_prob=probability
                        best_word=word               #for each line best_word keeps the best k-mer
                motifs.append(best_word)            #motifs keeps the best k-mer of each line and loop returns back to top.
        else:
            return motifs,score      # if score does not improve function returns the best motifs and score.

#returns the consensus of motifs
def consensus(motifs,k):
    counters=[]
    consensus=""
    for i in range(0,k):
        arr=[]
        for motif in motifs:
            arr.append(motif[i])
        counter=collections.Counter(arr)
        counters.append(counter)
    for counter in counters:
        consensus=consensus+counter.most_common(1)[0][0]
    return consensus
#%%
ten_mer=k_mer_generator(10)
ten_mers=k_mer_mutator(ten_mer)
txt_generator(10,ten_mers)
#%%
file=open('input.txt', 'r')
lines=file.readlines()
#%%
motifs,score=randomized_motif_search(lines,9)
##### 9-mer Motifs #####
print("9-mer motifs:")
for motif in motifs:
    print(motif)
consensus_motif=consensus(motifs,9)
print("Consensus: "+consensus_motif)
print("Score: "+str(score))
print()
motifs_2,score_2=randomized_motif_search(lines,10)
##### 10-mer Motifs #####
print("10-mer motifs:")
for motif in motifs_2:
    print(motif)
consensus_motif=consensus(motifs_2,10)
print("Consensus: "+consensus_motif)
print("Score: "+str(score_2))
print()
motifs_3,score_3=randomized_motif_search(lines,11)
##### 11-mer Motifs #####
print("11-mer motifs:")
for motif in motifs_3:
    print(motif)
consensus_motif=consensus(motifs_3,11)
print("Consensus: "+consensus_motif)
print("Score: "+str(score_3))
