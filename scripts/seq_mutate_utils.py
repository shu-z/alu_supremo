#!/usr/bin/env python
# coding: utf-8



'''
Functions that accompany SuPreMo, to use after get_seq. Mutates GC content and motifs in sequences output from get_seq. 

'''

# # # # # # # # # # # # # # # # # # 
# # # # # Import packages # # # # #

import pandas as pd
import numpy as np

import os
import io


# # # # # # # # # # # # # # # # # # 


def mutate_gc(seq, variant_start, variant_end, posflank, endflank, mut_percent):

    if not (0 <= mut_percent <= 100):
        raise ValueError("Mutation percentage must be between 0 and 100")
        
    if (variant_start-posflank)<0 or (variant_end+endflank)>len(seq):
        raise ValueError("Flanking indices exceed available sequence length")
    
    mutated_sequence = list(seq)
    if posflank>0:
        mut_start=variant_start-posflank
        mut_end=variant_start
    if endflank>0:
        mut_start=variant_end
        mut_end=variant_end+endflank
    
    else:
        mut_start=variant_start
        mut_end=variant_end
    
    
    for i in range(mut_start, mut_end + 1):
        if random.randint(1, 100) <= mut_percent:
            current_nucleotide = mutated_sequence[i]
            if current_nucleotide in ['G', 'C']:
                mutated_sequence[i] = random.choice(['A', 'T'])
                


    return ''.join(mutated_sequence)

    
def gc_seq(seq, start, end):
    seq_substr=seq[start:(end+1)]
    gc_count = seq_substr.count('G') + seq_substr.count('C')
    return (gc_count / len(seq_substr)) 
    