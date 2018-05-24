import operator

def occurrence(string, length):
    kmers = {}
    
    for i in range(0,len(string) - length + 1):
        short = string[i:i+length]
        
        if short not in kmers:
            kmers[short] = 0
        
        kmers[short] += 1
        
    return kmers


def seek_relevant_occurrence(dna, min_length, max_length, threshold):
    results = {}
    
    for i in range(min_length, max_length):
        temp = occurrence(dna, i)
        sorted_temp = sorted(temp.items(), key=operator.itemgetter(1))
        
        i = 1
        while sorted_temp[-i][1] >= threshold:
            results[sorted_temp[-i][0]] = sorted_temp[-i][1]
            i += 1
        
    return results