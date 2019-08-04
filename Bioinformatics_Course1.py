# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 22:13:08 2019

@author: Kathryn
"""

folder = "data"

sample = open(folder+'/dataset_2_7.txt', 'r').read().splitlines()

sample

text = sample[0]
pattern = sample[1]

#count how often a particular sequence motif occurs in a range of text
#text is the sequence (string)
#pattern is the motif (string)
def PatternCount(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1): #want to shorten by the length of the pattern so you don't consider things shorter than the pattern
        if pattern == text[i:(i+len(pattern))]:
            count += 1
    return count

PatternCount(text, pattern)

#want to find the most frequently occuring motifs in a sequence
#text (string) is the sequence
#k (integer) is how long the motif is
def FrequentWords(text, k):
    freq_patt = {} #use a dictionary to store all the motifs and add counts
    count = 0
    for i in range(len(text)-k+1): 
        pattern = text[i:i+k] #identify the motif
        if pattern in freq_patt: #save computational power so you don't count for the same pattern multiple times
            break
        else:
            counti = PatternCount(text, pattern) #cycle through sequence to see how frequently it occurs
            freq_patt[pattern] = counti #set the count to the dictionary entry for that pattern
    pat = []
    m = max(freq_patt.values())
    for key in freq_patt:
        if freq_patt[key]==m:
            pat += [key]
    return pat

text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"

FrequentWords(text, 4)


sample = open(folder+'/dataset_2_10.txt', 'r').read().splitlines()
text = sample[0]
k = sample[1]

FrequentWords(text, 12)

#Make the reverse complement of a DNA sequence
#text1 is the sequence (string)
def ReverseComplement(text1):
    output = ""
    text = text1[::-1] #reverse string
    for i in range(len(text)): #replace the nucleotide with its base pair
        if text[i]=="A":
            output += "T"
        if text[i]=="T":
            output += "A"
        if text[i] == "C":
            output+="G"
        if text[i]=="G":
            output+="C"
    return output

sample = open(folder+'/dataset_3_2.txt', 'r').read().splitlines()
text = sample[0]

ReverseComplement(text)

#print at what positions in a string a pattern occurs
#pattern (string) is the motif of interest
#text (string) is the sequence
def PatternMatching(pattern, text):
    pos = []
    for i in range(len(text)-len(pattern)+1): #go through the sequence
        if text[i:i+len(pattern)]==pattern: #see if this portion of the text matches the pattern
            pos.append((i)) #if so, document position
    print(*pos, sep=" ") #print out the positions separated by a space

pattern = "ATAT"
text = "GATATATGCATATACTT"

PatternMatching(pattern, text)

sample = open(folder+'/dataset_3_5.txt', 'r').read().splitlines()
pattern = sample[0]
text = sample[1]
PatternMatching(pattern, text)


#Find k-mers that are repeated at least t times in a l length clump
#k (integer) is the lenght of the motif we're looking for
#genome is the sequence (string) that contain potential motifs
#t (integer) is how frequently a motif must appear to be considered
#l (integer) is the length of sequence that we're looking at for the repetition
def Clump(genome, k, l, t):
    dic = {}
    sequences = {}
    for i in range(l-k+1): #for the first area with the length of interest, record all the morifs and how frequently they occur
        if genome[i:(i+k)] in dic: 
            dic[genome[i:(i+k)]] += 1
        else:
            dic[genome[i:(i+k)]] = 1
        #print(dic)
    for i in range(l, len(genome)): #now want to look at the remainder of the genome after the first l basepairs
        for key in dic:
            if dic[key] >= t: #if a key has at least t occurrences in the first length l, then add it to the sequences dictionary
                sequences[key] = ""
                #print(sequences)
                #print(dic)
        dic[genome[i-l:(i-l+k)]] -= 1 #now going to go after l basepairs and need to remove a sequence count that's now out of that range.
        if genome[i-k+1:(i+1)] in dic: #for the new sequence, see if it's already in the dictionary. If so, add a count.
            dic[genome[i-k+1:(i+1)]] += 1        
        else: #if it's not in the dictionary, add it to the dictionary
            dic[genome[i-k+1:(i+1)]] = 1
    print(*list(sequences.keys()), sep=" ") #print out all the sequences that are "clumped"
            

sample = open(folder+'/dataset_4_5.txt', 'r').read().splitlines()
text = sample[0]
sample[1]

Clump(text, 8, 28, 3)

ecoli = open(folder+'/E_coli.txt', 'r').read().splitlines()
ecoli1=ecoli[0]
ecoli1[0:30]
Clump(ecoli1, 9, 500, 3)

#Assign a pattern to a number to supposedly save memory
#pattern is a sequence (string)
def PatternToNumber(pattern):
    pattern = pattern[::-1] #reverse the sequence
    q = len(pattern) - 1 #set a number to count down for the number of iterations
    number = 0
    while q >= 0 : #different sequence values equal different numbers
        if pattern[q] == 'A':
            num = 0
        if pattern[q] == 'C':
            num = 1
        if pattern[q] == 'G':
            num = 2
        if pattern[q] == 'T':
            num = 3
        number += (num * (4**q)) #keep adding a value for each sequence and position
        q -= 1
    return(number)

pattern = "ATGCAA"
print(PatternToNumber(pattern))


#Return a sequence for each number
#index is an integer that stands for the sequence
#k is how long the sequence is (integer)
def NumberToPattern(index, k):
    num = index
    rem = 0 #the remainder
    pattern = ""
    q = k - 1 #set how many iterations to cycle through
    while q >= 0:
        rem = num%4 #get the remainder of the input number
        num = num/4 #reset the next input number to be the previous value divded by 4
        if rem < 1: #depending on what the remainder is, get a different letter output
            pattern += "A"
        if rem <2 and rem >=1:
            pattern += "C"
        if rem <3 and rem >= 2:
            pattern += "G"
        if rem < 4 and rem >=3:
            pattern += "T"
        q -= 1
    pattern2 = pattern[::-1] #reverse the pattern
    return pattern2

NumberToPattern(5437, 7)
NumberToPattern(5437, 8)

#An alternate way to count how many occurences of a pattern there are
#text is the sequence of interest (string)
#k is how long of a motif you're interested in
def ComputingFrequencies(text, k):
    freq = {}
    pattern = ""
#    list = []
    for i in range(4**k):
        freq[i] = 0 #set up dictionary with all possible values
    for i in range(len(text)-k+1): 
        pattern = text[i:i+k]
        j = PatternToNumber(pattern) #turn all patterns into numbers
        freq[j]+=1 #add a count for each detected pattern that you generated earlier
#    for i in range(4**k):
#        list.append(freq[i]) 
    return freq

ComputingFrequencies("ACGCGGCTCTGAAA", 2)

sample = open(folder+"/dataset_2994_5.txt", 'r').read().splitlines()
text = sample[0]
k = sample[1]

x = ComputingFrequencies(text, 5)
f = open(folder+"/answer.txt", 'w')
f.write()


#Identify frequent motifs in a sequence faster by converting them into numbers
#text is the sequence (string)
#k is the length of the motif (integer)
def FasterFrequentWords(text, k):
    freq_patt = []
    freq_array = ComputingFrequencies(text, k)
    max_count = max(freq_array.values())
    for i in range(4**k):
        if freq_array[i] == max_count: #if one of the integers stands for a motif that has a high count, change it back into a pattern
            pattern = NumberToPattern(i, k)
            freq_patt.append(pattern)
    return freq_patt

FrequentWords(text, 4)
FasterFrequentWords(text, 4)

#Find the last letter in a pattern (string)
def LastSymbol(pattern):
    pattern2=pattern[len(pattern)-1]
    return pattern2

#Find everything but the last letter in a pattern (string)
def Prefix(pattern):
    pattern2 = pattern[0:(len(pattern)-1)]
    return pattern2

#for a certain letter, return a number
def SymbolToNumber(symbol):
    if symbol == "A":
        return 0
    if symbol == "C":
        return 1
    if symbol == "G":
        return 2
    if symbol == "T":
        return 3


#Use sorting to find frequent motifs. Return list.
#Text is sequence input (string)
#k is length of motif (integer)
def FindingFrequentWordsBySorting(text, k):
    freq_patt = []
    count = {}
    for i in range(0, len(text)-k+1):
        pattern = text[i:i+k] #collect patterns in sequence of length k
        Index = PatternToNumber(pattern) #convert patterns to integer
        if Index in list(count.keys()): #if this pattern already is a key in the dictionary, add 1 to it
            count[Index] += 1
        else: #if pattern not already a key in the dictionary, set it as a key with value 1
            count[Index]=1
    SortedIndex = sorted(count) #sort dictionary keys ascending order
    max_count = max(count.values()) #find highest dictionary values
    for key in count: 
        if count[key]==max_count: #if key is one of the most frequently occurring
            pattern = NumberToPattern(key, k)  #convert number back to pattern
            freq_patt.append(pattern) #add to list of frequent patterns
    return freq_patt

FindingFrequentWordsBySorting(text, 4)    

#Find k-mers that are repeated at least t times in a l length clump. Supposedly a bit faster because using numbers instead of letters.
#k (integer) is the lenght of the motif we're looking for
#genome is the sequence (string) that contain potential motifs
#t (integer) is how frequently a motif must appear to be considered
#l (integer) is the length of sequence that we're looking at for the repetition
def BetterClumpFinding(genome, k, l, t):
    freq_patt = []
    clump = {}
    for i in range(0, 4**k-1):
        clump[i]=0
    text = genome[0:l]
    freq_array=ComputingFrequencies(text,k)
    for key in clump: #count motif repetitions
        if freq_array[key] >= t:
            clump[key] += 1
    for i in range(1, len(genome)-l):
        first_patt = genome[i-1:i-1+k] 
        index = PatternToNumber(first_patt)
        freq_array[index] -= 1 #remove first pattern encountered from counts
        last_patt = genome[i+l-k:i+l]
        index = PatternToNumber(last_patt)
        freq_array[index] += 1 #add new pattern to counts
        if freq_array[index] >= t: #for things that appear above the threshold, set value to 1
            clump[index] = 1
    for key in clump: 
        if clump[key] ==1:
            pattern = NumberToPattern(key, k) #translate numbers back to patterns
            freq_patt.append(pattern) 
    return freq_patt
            

#Determine if the sequence skews more G or C
#genome is a sequence (string)
def Skew(genome):
    skew = 0
    G = 0
    C = 0
    result = [0]
    for i in range(len(genome)): #keep a running count of G or C occurrences
        if genome[i] == "G":
            G +=1
        if genome[i] == "C":
            C +=1
        skew = G - C 
        result.append(skew) #adjust skew count for each base pair
    return result

result = Skew("GAGCCACCGCGATA")
print(*result, sep=" ")


#return the position (integer) where the skew is the lowest
#input is the sequence of interest (string)
def MinSkew(genome):
    skew = Skew(genome)
    low = min(skew)
    pos = []    
    for i in range(len(skew)):
        if skew[i]==low:
            pos.append(i)
    return pos

MinSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")

sample = open(folder+"/dataset_7_6.txt", 'r').read().splitlines()
text = sample[0]

MinSkew(text)

#Calculate how different two sequences are from each other. The bigger the number, the more different
#p and q are both sequence strings that should be the same length
#Return an integer
def HammingDistance(p, q):
    dist = 0
    if len(p)==len(q): #check to make sure sequences are the same length
        for i in range(len(p)):
            if p[i]==q[i]: #if the sequence has the same letter at the same position, don't add to score
                dist+=0
            else:
                dist+=1
    return dist



sample=open(folder+"/dataset_9_3.txt", 'r').read().splitlines()
p = sample[0]
q=sample[1]

HammingDistance(p,q)


#Identify the position in a sequence where a pattern has less than d mismatches
#pattern is the motif of interest (string)
#text is the sequence you're probing (string)
#d is the number of allowable mismatches (integer)
#Output is an integer that indicates the position of the pattern match
def ApproxPatternMatch(pattern, text, d):
    result = []
    n = len(pattern)
    for i in range(len(text)-n+1):
        x = HammingDistance(pattern, text[i:i+n])
        #print(x)
        if x <= d:
            result.append(i)
    return result


ApproxPatternMatch("ATTCTGGA","CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3)

ApproxPatternMatch("CCC", "CATGCCATTCGCATTGTCCCAGTGA", 2)

sample=open(folder+"/dataset_9_4.txt", 'r').read().splitlines()
pattern = sample[0]
text=sample[1]
sample[2]
ApproxPatternMatch(pattern, text, 6)



#Need to make all possible iterations with the 4-nucleotides
import itertools

bases=['A','T','G','C']
x = [''.join(p) for p in itertools.product(bases, repeat=4)] #Make all possible 4-mers

#Find the most frequent k-mers that have less than d mismatches
#text is the sequence of interest (string)
#k is how long the motif (or k-mer) is (integer)
#is the number of allow mismatches
#returns a list of k-mers
def FrequentMismatch(text, k, d):
    bases=['A','T','G','C']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)] #generate all possible k-mers of a certain length
    dic = {}        
    for j in range(len(kmers)): #cycle through every kmer
        for i in range(len(text)-k+1): #compare that k-mer to every position in the sequence
            x = HammingDistance(kmers[j], text[i:i+k]) #see how different that k-mer is from the sequence
            if x <= d: #if there's less than d mismatches
                if kmers[j] in dic: #if it's already in the dictionary, add to it's score
                    dic[kmers[j]]+=1
                else: #if not, add it to the dictionary and give it a score of 1
                    dic[kmers[j]]=1
    keys = []
    for key in dic: #now find the keys with the highest scores
        if dic[key]==max(dic.values()):
            keys.append(key)
    return keys

text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
   
FrequentMismatch(text, 4, 1)

sample = open(folder+"/dataset_9_7.txt", 'r').read().splitlines()
text=sample[0]
sample[1]

x = FrequentMismatch(text, 6, 3)
print(*x, sep=" ")

#Similar to Frequent Mismatch, but this also considers the reverse complement
#(DNA strings are double stranded, with one sequence going one direction and the reverse complement going the other)
#If a motif matches the reverse complement, it's still relevant
#Returns a list of the most frequent motifs of length k with less than d mismatches
def FrequentMismatchRC(text, k, d):
    bases=['A','T','G','C']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    dic = {}        
    RCText = ReverseComplement(text)
    for j in range(len(kmers)):
        for i in range(len(text)-k+1):
            x = HammingDistance(kmers[j], text[i:i+k])
            y = HammingDistance(kmers[j], text[i:i+k])
            if x <= d:
                if kmers[j] in list(dic.keys()):
                    dic[kmers[j]]+=1
                else:
                    dic[kmers[j]]=1
        for i in range(len(RCText)-k+1):
            y = HammingDistance(kmers[j], RCText[i:i+k])            
            if y <= d:
                if kmers[j] in list(dic.keys()):
                    dic[kmers[j]]+=1
                else:
                    dic[kmers[j]]=1                
    keys = []
    for key in dic:
        if dic[key]==max(dic.values()):
            keys.append(key)
    return keys

test = "ACGTTGCATGTCGCATGATGCATGAGAGCT"

FrequentMismatchRC(test, 4, 1)

FrequentMismatch(test, 4, 1)

sample=open(folder+"/dataset_9_8.txt", 'r').read().splitlines()
text=sample[0]
sample[1]
FrequentMismatchRC(text, 5, 2)


#So, the previous FrequentMismatches are grossly inefficient if you have a k-mer that's not super short
#The general idea of this function is that instead of generating all possible k-mers,
#you just generates ones that are d-number of mismatches from the initial pattern
#To be truly honest, I copied this code and do not understand the recursion element of it.
#returns a list of strings with d-number of mismatches from the inital pattern (string)
def Neighbors(pattern, d):
    if d == 0: #if you're allowed no mismatches, just return the pattern
        return pattern
    elif len(pattern) == 1: #if the pattern is only 1nt long, return all the possible nucleotides
        return {'A', 'C', 'G', 'T'}
    neighborhood = []
    suffix_neighbors = Neighbors(pattern[1:], d) #recursion. Logically goes off the rails for me here.
    for text in suffix_neighbors:
        nucleotides = {'A', 'C', 'G', 'T'}
        if HammingDistance(pattern[1:], text) < d:
            for x in nucleotides:
                neighborhood.append(x + text)
        else:
            neighborhood.append(pattern[0] + text)
    return neighborhood


#This is all practice for me trying to figure out how to do MotifEnumeration
prac = ["dog", "liger", "kitten", "dog"]
prac2 = {"dog":1, "liger":0, "kitten":3}
for key in prac2:
    if prac2[key]==0:
        prac.remove(key)
        print(prac)

list(set(prac))
prac.remove("dog")

prac["dog"]=[]

prac["dog"].append(1)
prac["zebra"]=1

prac.defaultdict("dog", []).append(1)
prac.remove("liger")
##end of practice


#Find motifs of k-length that occur in each dna string
#dna is a list of strings. Can also just be a string. 
#k is an integer to indicate length
#d is an integer for the number of mismatches allowable
#Returns motifs. motifs is a list of strings (each k-length)
def MotifEnumeration(dna, k, d):
    motifs = []
    for i in range(len(dna)): #cycle through the dna list
        if i == 0: #for the first string
            for j in range(len(dna[i])-k+1): #cycle through the first string
                pattern = (dna[i][j:j+k]) #identify a sequence of k-length in the first string
                pattern_p=Neighbors(pattern, d) #use that pattern to identify all possibilities with d mismatches
                if isinstance(pattern_p, str)==True:  #If there's only one value, it's a string not a list.
                    motifs.append(pattern_p) #add to potential motif list
                else:
                    motifs.extend(pattern_p) #add to potential motif list
    motifs = list(set(motifs)) #get rid of duplicates
    temp = {}
    removemotif = []
    for j in range(1, len(dna)): #now look at the other dna strings
        for motif in motifs: #cycle through the motifs and see if they're in the other dna strings
            #print(motif) #checking that things are working
            for i in range(len(dna[j])-k+1):
                if HammingDistance(motif, dna[j][i:i+k]) <= d: #identify if motif occurs in the string
                    if motif in list(temp.keys()):
                        temp[motif]+=1
                    else:
                        temp[motif]=1
                else:
                    if motif not in list(temp.keys()):
                        temp[motif]=0
            if temp[motif] == 0: #identify motifs that should be removed because they don't occur in the string
                removemotif.append(motif)
        #print(motifs) #checking that things are working
        temp = {} #reset temp before cycling to the next dna string
    for motif in list(set(removemotif)): #remove all the motifs not in every string
        motifs.remove(motif)
    return motifs

test = ["ACGT", "ACGT", "ACGT"]

Neighbors("ACG", 0)

MotifEnumeration(test, 3, 0)

text = ["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]
MotifEnumeration(text, 3, 1)

text = open(folder+"/dataset_156_8.txt", 'r').read().splitlines()
thing = text[1:len(text)-1]
text[0]

MotifEnumeration(thing, 5, 1)


#Count how frequently a symbol occcurs at a certain position when given a list of motifs of the same length
#motifs is a list of strings of the same length
#count output is a dictionary with A C G T as the keys
def Count(motifs):
    count = {} # initializing the count dictionary
    k = len(motifs[0]) #get the length of each motif row (x axis)
    for symbol in "ACGT": #Make A C G T a key
        count[symbol] = [] #make an empty list for each key
        for j in range(k): #for the length of the motif...
            count[symbol].append(0) #...for the key, put a zero for each position in the list for now
    t = len(motifs) #get the number of keys in the dictionary (y axis)
    for i in range(t): #for every key...
        for j in range(k): #for list position in each key...
            symbol = motifs[i][j] #for this specific position for a key...
            count[symbol][j] += 1 #Add one every time this position occurs in the matrix at this position
    return count

#similar to Count, but get a ratio instead. All columns should add to 1.
#motifs input is a list of strings of the same length
#profile output is a dictionary with A C G T as keys
def Profile(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = Count(motifs)
    sub= Count(motifs)
    for x in range(t): #for key in position x...
        for y in range(k): #for list position in that key...
            symbol=motifs[x][y]
            profile[symbol][y]=(sub[symbol][y])/t
    return profile

motifs = [
"TCGGGGGTTTTT",
"CCGGTGACTTAC",
"ACGGGGATTTTC",
"TTGGGGACTTTT",
"AAGGGGACTTCC",
"TTGGGGACTTCC",
"TCGGGGATTCAT",
"TCGGGGATTCCT",
"TAGGGGAACTAC",
"TCGGGTATAACC"
]

Profile(motifs)


#Need this for calculating log2 values
import math

#Calculate the randomness of a motif
#uses a list of dna strings as the input
#output is a number standing for the entropy
#This wasn't a very useful function.
def entropy(motifs):
    prof = Profile(motifs)
    #result = []
    result2 = 0
    for i in range(len(prof["A"])):
        temp = 0
        for key in prof:
            if prof[key][i] != 0: #don't want to take the log of zero values
                temp += abs(prof[key][i]*math.log2(prof[key][i]))
            else:
                temp += 0
        #result.append(temp) #this is really just here to see that things were working
        result2 += temp
    return result2

0.5*math.log2(0.5)*2

0.25*math.log2(0.25)*4

1*math.log2(1)

0.25*math.log2(0.25)*2+0.5*math.log2(0.5)

entropy(motifs)




#Identify the motif with the lowest score
#pattern is a string
#dna is a list of strings
#return a list with the best match for each dna string 
def Motifs(pattern, dna):
    motifs = []
    k= len(pattern)
    if isinstance(dna, str)==True: #verify if the input is a list or a string
        new_dna = []
        new_dna.append(dna)
        dna = new_dna #If a string, make it a list so it works with the rest of the function
    for seq in dna:
        score = k
        temp_motif = {}
        for i in range(len(seq)-k+1): #find the motifs with the lowest score
            new_score = HammingDistance(pattern, seq[i:i+k])
            if new_score < score:
                temp_motif[seq[i:i+k]]=new_score
        if len(temp_motif)==0: #if you get nothing lower, just return the first k-mer
            motifs.append(seq[0:k])
        else:
            low_score = min(temp_motif.values()) #Find the lowest score
            for key in temp_motif: #identify the motif with the lowest score
                if temp_motif[key] == low_score:
                    motifs.append(key)
                    break #just want one motif per dna string, so end search once you find it
    return motifs

dna = ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG"]    

Motifs("AAA", dna)

Motifs("GATTCTCA", "GCAAAGACGCTGACCAA")


#Return a score for how different a pattern is from the dna input
#pattern is a string
#dna is a list of strings
#output score is an integer
def d(pattern, dna):
    motifs = Motifs(pattern, dna)
    if len(motifs)==0:
        return len(pattern)
    score = 0
    for motif in motifs:
        x = HammingDistance(pattern, motif)
        score += x
    return score

d("AAA", dna)

d("GAC", dna)

#Search for the k-mer that is in at least one string and has the lowest score for all strings
#dna can be either a string or a list of strings.
#k is an integer for how long the output is supposed to be
#output is median, which is a string
def MedianString(dna, k):
    kmers = []
    median = ""
    if isinstance(dna, str)==True: #If dna is a string, make it a list
        new_dna = []
        new_dna.append(dna)
        dna = new_dna 
    distance = len(dna[0])*len(dna)
    for seq in dna: #for each string in the dna list
        for i in range(len(seq)-k+1): #make a list of all kmers
            kmers.append(seq[i:i+k])
    kmers = list(set(kmers)) #remove duplicates
    for kmer in kmers: 
        if distance > d(kmer, dna): #if the differences are less than the total length
            distance = d(kmer, dna) #update the d output as the new min
            median = kmer #return the kmer with the lowest distance
        elif distance == d(kmer, dna):
            median = kmer
    return median

dna = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"]

MedianString(dna, 3)

#Similar to previous, but search kmers not in the main text (uses Neighbors)
#output is now a list of strings
def MedianString2(dna, k):
    kmers = []
    median = []
    if isinstance(dna, str)==True: #If dna is a string, make it a list
        new_dna = []
        new_dna.append(dna)
        dna = new_dna 
    distance = len(dna[0])*len(dna) #technically, an entry could poorly match everything
    for seq in dna:
        for i in range(len(seq)-k+1):
            kmers.append(seq[i:i+k]) #make a list of all kmers in the first dna string
    kmers = list(set(kmers)) #remove duplicates
    kmers2 = []
    for i in range(len(kmers)):#add more kmers that are not just in the main text
        kmers2.append(kmers[i])
        newkmers = Neighbors(kmers[i], k%2+1) #arbitrarily decided on this since no number of mismatches is stated. Essentially, adding things that can over half mismatches.
        kmers2.extend(newkmers) #add new kmers to the list
    kmers2 = list(set(kmers2)) #remove duplicates
    for kmer in kmers2:
        if distance > d(kmer, dna): #if the new kmer has a small distance...
            distance = d(kmer, dna) #reset the distance to the new min
            median = [] #reset the median to get rid of anything that performed worse
            median.append(kmer)
        elif distance == d(kmer, dna): #if distance is the same, add it to the list so we know there are multiple motifs that perform equally well
            median.append(kmer)
    return median

test = ["ATA","ACA","AGA","AAT","AAC"]

MedianString2(test, 3)

MedianString2(dna, 3)

sample = open(folder+"/dataset_158_9.txt", 'r').read().splitlines()
text = sample[0]
text2 = sample[1]
dna = ['TGTTGGAAATATCACGTGGACCCTCTCGTGCAGTGGATATAA',
 'TGTCCGCGGCAACACCCTAGTGCCGGAGCCGGATCTGTGCAA',
 'TGGGCGTGGACGTACCCTGACCCCCAATTTCGTGAGTATCAA',
 'CTTCATAATTATAACCCTATTTTCAGTCACGAGCACCCTAAG',
 'TACCCTCCAATACAGGACGTACCTAGGTGCGCGCGCCAGACG',
 'TGTATTAACCAAAGATAGTCCCGGTACCCTGCGCGCTCGGCA',
 'AACCCTCCACTCCGGAATAATGTTTACATACCCCCTTTAATT',
 'CTGCTAGCATAGAGCATTGCTTGTTTGGCTGCCCTACACCCT',
 'GTGCAGGGCGTAAGATAGAACCCTGCTACGAGCCGTGGAAAC',
 'CAGAACACCCGCGCAAGAGACCCTTCGGTAGGAATCATGCCT']

MedianString(dna, 6)
MedianString2(dna, 6)

#Find the likelihood of a certain sequence occuring within a profile of dna strings (probability of each nucleotide at each position)
#text is a dna string
#profile is a dictionary with keys A, C, G, T and the same number of entries as text
#output is an integer
def Pr(text, profile):
    p=1
    for i in range(len(text)):
        p=p*profile[text[i]][i]
    return p

profile = {

    'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
}

Pr("TCGTGGATTTCC", profile)

#Identify the most probably kmer from the profile thats in the current text
#text is a dna string
#k is the length of the motif (integer)
#profile is a diciontary with keys A C G T
#output is a string of length k
def ProfileMostProbableKmer(text, k, profile):
    score=-1 #start negative because you want to go as close to 1 as possible, but result could still be 0
    pattern=text[0:k] #set first pattern
    for i in range(len(text)-k+1): #cycle through the dna string
        score2=Pr(text[i:i+k], profile) #get probability of pattern occuring in the text
        if score2>score: #if probability is higher than previous, update the threshold
            score=score2
            pattern=text[i:i+k] #choose the next pattern
        if score2==score: #if score is the same, update the pattern
            pattern=pattern
    return pattern

profile = {
    'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    'T': [0.1, 0.2, 0.1, 0.1, 0.2]
}

ProfileMostProbableKmer("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT",5, profile)

sample = open(folder+"/dataset_159_3.txt", 'r').read().splitlines()
text=sample[0]

profile = {
         "A": [0.303, 0.237, 0.25, 0.197, 0.171, 0.25, 0.342, 0.237, 0.092, 0.224, 0.289, 0.303, 0.289],
 'C' : [0.184, 0.316, 0.289, 0.276, 0.329, 0.237, 0.079, 0.303, 0.263, 0.276, 0.25, 0.237, 0.303],
 'G' : [0.289, 0.237, 0.224, 0.25, 0.237, 0.25, 0.355, 0.211, 0.25, 0.276, 0.237, 0.211, 0.263],
 'T' : [0.224, 0.211, 0.237, 0.276, 0.263, 0.263, 0.224, 0.25, 0.395, 0.224, 0.224, 0.25, 0.145]}

ProfileMostProbableKmer(text, 13, profile)

#Previous counts had a problem with 0s. Made calculating probability 0 for everything.
#This will never leave something as 0.
#takes an input of the list of strings motifs
#output is an dictionary with counts of occurrences at each position, but instead of lowest value being 0, it's 1
def CountWithPseudocounts(motifs):
    count = {} # initializing the count dictionary
    k = len(motifs[0]) #get the length of each motif row (x axis)
    for symbol in "ACGT": #Make A C G T a key
        count[symbol] = [] #make an empty list for each key
        for j in range(k): #for the length of the motif...
            count[symbol].append(1) #...for the key, put a one for each position in the list for now
    t = len(motifs) #get the number of keys in the dictionary (y axis)
    for i in range(t): #for every key...
        for j in range(k): #for list position in each key...
            symbol = motifs[i][j] #for this specific position for a key...
            count[symbol][j] += 1 #Add one every time this position occurs in the matrix at this position
    return count

#Similar to the previous profile function, but using pseudocounts instead of normal counts
#Once again, eliminates 0s
#takes an input of the list of strings motifs
#Output is a dictionary with keys A C G T with fractions for the probability of occurrence for each nucleotide at a sequence
def ProfileWithPseudocounts(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = CountWithPseudocounts(motifs)
    sub= CountWithPseudocounts(motifs) #setting this twice to make it less confusing for the below math
    for symbol in "ACGT":
        for y in range(k):
            profile[symbol][y]=float(sub[symbol][y])/(t+4)
    return profile

#Identify the most common motif from a list of motifs
#motif is a list of strings
#returns a string
def Consensus(motifs):
    k = len(motifs[0])
    count = CountWithPseudocounts(motifs) #get a dictionary for how often nucleotide occurs at each position
    consensus = ""
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol
    return consensus

#Score how different a motif is from the consensus sequence
#motifs is a list of strings
#returns an integer, the higher, the more different from the consensus
def Score(motifs):
    con = Consensus(motifs)
    k = len(motifs[0])
    j=len(motifs)
    value=0
    for x in range(k):
        for y in range(j):
            if motifs[y][x]!=con[x]:
                value+=1
    return value

#A slightly faster way to find good motifs, but not without its problems (misses things not in the sequence)
#dna is a list of strings
#k is an integer for the length of the motif of interest
#t is the number of dna strings (I was asked to write it like this, it's really not necessary)
#output is a list of strings, all of length k
def GreedyMotifSearch(dna, k, t):
    best_motifs=[]
    for i in range(0, t):
        best_motifs.append(dna[i][0:k]) #sets best_motifs to the first three k-mers in each row
    n=len(dna[0]) #how long each row is
    for i in range(n-k+1): #iterate through each row as a k-mer
        motifs=[] #set up an empty motif matrix
        motifs.append(dna[0][i:i+k]) #fill motif matrix with k-mer
        for j in range(1,t): #go through each row in the matrix
            p=ProfileWithPseudocounts(motifs[0:j]) #give a score for each nucleotide for the motif
            motifs.append(ProfileMostProbableKmer(dna[j],k,p)) #For the j row, get the motif with the highest score using the profile scores
        if Score(motifs) < Score(best_motifs): #If the new kmers score better than the artbitrary ones chosen earlier, set them as the best_motifs
            best_motifs=motifs
    return best_motifs

sample=open(folder+"/dataset_159_5.txt", 'r').read().splitlines()
dna = ['GTACAAGTGAACGGCTGTTTACAGTCGGAGTAGCAGGTCTTCTGTTATTGGCGATGGCCGGTATATACGAAAGGCGAAGCGTGGTGCGATCTCTTTGCCTGTGCTCTCCGAATCCCTAGGGAGCGCAATCATAGGCTCTCCCTGCAAATGGAGCGT',
 'TAAACATAAGCCACCGGCTCCCTCTTTAAGCACAGCCCACCGCACAGAACCGATGCGTCCGATGAAATCAGTGACTGGTGGAAGTTAAGAGCTCTAGGCTGTCGACGCGCCTCGCTCTCGGTCCCCCCTTGAAGGGCCTCCCTGTAGTCAGGACAA',
 'AAACTCGCTGAGTTTAACTGCGCCGTCGTCTCTGGTTCAGGTGTGCTAAAGCCAGGTCCATTTGGCACTTAGACATGCACTCACATTTTCCTTCTGCGGACCTCCCTGTGGCACTTATAGCAAGAGGATTCACCCGAAGGACCGTCCGGAAAAATA',
 'CGGCACCTACACTTTACCCTGGGATGAGTTGCCTCATAATGTGTTGATAATGGCATATAGTCCAAACTGCCCAGGCGCTCCATGCAACTAATCGTGCATCTTTCCCTGTCTCTCCGTGCACGATAGTGTGGAACGGACAGGACGATGCGCTCGAGC',
 'GTAGCCCTCCAGAAGACTCCGTTCCAAGTTATCCGATTACCCAGAATTGTACGATAGCTACTCTGCTGCTCCCTTTTCGCGTCCAATGCGATATCGAATGAACCCGGTGGCTACATTACCTAGCCCTCCGAGTGTTTGTTACGTGGGATCTCCGTG',
 'TGAAGTAAAAGGCGGGTCTCCTTGACCGATTTAATGGTACTGTAGACGCAAGGGTGAAGCCCTGTGGCCGTAAGCACCCTACGCCAGACGCGGCATTTTGGCGGTTTTCGTTTAGGTTCTAAAACGCAACTAGCTGTTACTAGCTTCATTACAGTG',
 'GGGGGTTTCTATCTCTTCAAGTGGCCCTGTTACTCTCCAGAGGTGCGGCAAAGTCATTAACGGGACGGGCTGCATTTCCGTCTGGCTAGCATCTGGTACGGGTGCGGCTACTCAGTGGTCATAGGTTGGGTGCAATTGACCTTTAGGGCCTCCTTG',
 'TGCCTATCGGCGCTTCTGGGGAACGGGTTCTCCCTGGCTAATCTGCAGAATGAAAGGTTGCTTTTACAGAATCTTCAAGCCGGGTGCTCCAGCGACAGAAATTTCCGTGAATATCCCATATGATGTACGGTAGTAATGTTCCCCTACTTTGTTCAA',
 'TACGATCCGCACCGTTTGGGACCGATCCTGGGTGAGAGCAGGGCTAGTACGCTGCTACCAGGTCCCTTATTCGGCCTCGTCGAACACACTAACAACTAACGTGCCCTACGCCGACTCCTCCATACGGTGTACAGGCTCTCCTTGCTTCCGGAAATC',
 'GCCCGAACATCTTACTTAGCAGGTTATCAGCAAAATCGATCACTATAGAAGTTCCTGGTTAAATAGGTTAAAAGGTTCGAGGATAAGGCTTGTAAGAGTTTGCTGCTTCGGAACTCCTTGGGAGAACAGCCGCTCACGGCAAAATAACAGCGAGAT',
 'AACTTTCGATTCGATAAGTCCCAGTGTTCGATCGAGCTTGTTTTCAGCGGGAGCTCCTTGGCCCGCAGTGTAAAGCGGTGAGACACCAATCTGTTCCAGCCGAAAAGATGGTAGAAGCGGACTAGGATTGTAGCTTTTCCCATTGAAAATTCGCTT',
 'ACTGCGTTTGCCACACCCGATTTATTAGACGGCGGGTAGCAACCGTGGACGAGAAATGAATATATTTTGTCAGAGTAGGGGGCTACTCGCAAATTCTCGCTAAGGCGATGTTCAGTTCCAGTGTCTAGCATCCGTTACCGCAGACGGATCTCCGTG',
 'AATATCGCACAACTAGAGCCTCGCTTTAAAAGTCAGTAGAGGACTCCCGGGGCCTCCATGGTCTTAAAGGGATCGTTAGAGTAAAGCACGGTACGCTATAGTGCCCAGGGGGATGGTTTTAAGGGAAAGCAACCCATTTGAGGGTACTAACCGAGC',
 'TTTGCTCCATGGCGATGTCTCCACTCGCTTCTCTCTCGACCCTAGGTCCGAAGAGCGGTCTATTCCAAGGAGTGGAGCTCCTTGGTTTCTGGCGCACTTTCTCAAACAAAATATAAGTGCCCTCCCCCGGGCATACGTGAATCTTCTAGAGTATTG',
 'TGGGCCTCCATGATTCCGCTTGGCACCGTCAGTGGACTTGAAAGCCAGTGGATCAAGGGGCAAAAAGGTAAAAAGCGTATATGTGGTGAGATATTTTCTTAAAACACACTTTTGTCTACAAAGACTACGGAAATTATCTGGGTCCTCTGGAGATCA',
 'TTCTTCAGGGGAATCGAGGCGCACCGGTCCTCCTTGGTTTCGCTAGACGCGACGGGTTGTCCCTTGAGGTGCAGAGAGCATCCACTCCAGGACGACCCCGGTTTAGGCTAGTATAAACTTGAGAAGCAATGACATACATGACAAGCGTAGCGCCTG',
 'GGGTCCTCCCTGGTCTTGCACGAAGTACCAGTTATGTGCTAACGTGGCAATAGGTCCGGGCCATACGTGCTTAACATATTAATTAGGGATGTTAACGGTACGAACAAGCGAGCAGCATATGACGCCTTCTTCGGGAACTTTTATACAGTGCAGCTG',
 'CGGTCCTCCCTGTTGGGTGTCGTTCTGCTGGATTAGGCATGGACGTCCCATCGCAAGCGCAGGCTATGGTCTTCACTCGGTAAGGTACCGGCTCCCAACACTAGACAGAGTCACCTAGAGTAGAATCGCGACCGTTGAAGCGCCGAGTTTTTTACT',
 'GATTGCCAAATGGGTTGGCCCTTATTGAGTGCAAAGATACCCATCACTTTCAGTATGCGTTGGATCTCCTTGATGACTTAGAGAGAGCTTGGCTAACGTTACTTAACTGTAAACGATCACCCTGGGACCCAAATCATGTACGAAATATACGCGATA',
 'GTAAAAGAGAGCGGACGCCGTATGTTCATCATCGCGGGGACCTCCGTGTAATGTAAGCGTAAAAAGGTCGGGGTTCCCCTTCAAGCTGTATTTTCATGGAGAGCCTTCGGTATGGGGTTATGACGAGTACGGAGCATCCGAATGCTTAATGGAGTC',
 'AGGAACACTGCGAGATCGTCTATTAGGACCGGACTGTACCTAGGGATAGTTTGTATGATACAGAGGGTCCTAGACGGATATGGGTCTGGCAATCCTAGGCTCTCCGTGGGTCACCGTGGTATAATCTTTAGGCATGACATCCCTTCCGATCCCGCC',
 'GCTTGTGCACGACCGTTACGTGAATAGCCGAGTGTGTGTGGCACTGTTATAGCTTATGTACAATTCCTTATTTATACAGCCGCACGACTAGCAGACGTCCTTCCATACAGGGGCTCCGTGGGTGAGTGTTAGGAAGCTGACATAAGGGTAGTATCT',
 'TTATTACTTCGCATGGTCTCAATGGAGACTGGAAACGAGAACCGCATTCAAACTTAGTTTCACGGTGGAGTGTAACCCCACTAAGTGCATCCTGTCATATATGCATATCGGACCTCCATGCCGTCAGCACTTGTCGGACAGGATTTAGTGGTAGCC',
 'GGGGTCTCCTTGTGGGGATAAACCATAGGCAGTCGAGCAACATAACCTCTTGGCTATGATATATGTTTAGACAACTGTTGCCCCTAATGCTGGTCCCATTCAAAAGTGTATGAATCGGCTAAAATTAGTCAGGCTGTCGATTACTACCCGTAGAGA',
 'CTCGCCTCTGTCTAGACCCCTTATCGTCCTACTATAGGCATCTCGATCTCCTGGAGGAGAGAATTACTATTCTATGTATCGGGAGGATAAAGGTAATTCACCACACTTCGTCAATTGGACGGGCACTCCCTGATGGTTTGAACCACACTCATGTGA']

text = GreedyMotifSearch(dna, 12, 25)

print(*text, sep=" ")

print(*list(sequences.keys()), sep=" ")

sample = open(folder+"/dataset_160_9.txt", 'r').read().splitlines()
text = ['GAGCCCCATTCCGCGCGATGCTCCTTACCGAGGGTTTCGCGAAAGGCTCACATGCTAGTTAGTTCGGCGAAGACAGCTGTGGGAGTTATCCGGCCCCTGTGCCCTAAGCAGTTTTGGATTACTGTCCGATTGATGCGCAGCTCAGACTCTTAGTCT',
 'TCACGGCCTGTGTTGATCGGTCCTAGAGTACAAGTACACCATACGAAAGAATACGTCCGTGCAGGCGGAACGCCATGTGTTAAATGCCCTGAAGATGTAATAAGGGAGGAGACTCTAGTTCCCTGCCACATACCCAACGTAACCAGTAGTTCCAAA',
 'ATTCCCCGAGTCGTTTCTAAAGACCGCACGTTAATTCGCGAGGGGAGAACTAGATAATGCAGCGAAACGGCGAACGCGCTAGAAGTTGAAGGTACCCCAGCCAGCGGTAACATATAAAGAATGTAATATCGTCTTACTCTCAGGCAAAGGCTAGTT',
 'GAGAAAAGAATATGACGCCCCTTACCGGGCGCAATATTGACCGTGTATGAAAAAACCCTTTTCGTGGTGAAGAATAATCTAGTTCCCCCACTGGTTGGACCGCTTGTGGTATGTGACTTGGAACTGTTGTAAGCGAAGTTTGAGGATTGCAGCTGG',
 'AAACCCTATCAAGGAAGTTGATACGCGACGAAGGCGATCGGGTAATCAAGACCGCGCCTTCGCTCATGATATTAACGCCGAAGCCAATCGGGTAGCGGTCACAATGTACGCACTTGGTGAGACACACTAGTTTTGAGTCAACTATAGTCGGTCTTA',
 'CTGAAAGCTCGTCACAACCTAGTTGGAGTATAATATACTGGAGCCGACGGGTCTGTAAGGTTGGTACGTCCCGCGACGCCAGCTTCATTCGGACTTTCCTTCGAACTAACATTACACTTCAGGGGAATTAGCCAGGAATTCGCGACGGGAGCACGA',
 'TATAACCTAGTTGTTGCTTTACAGTTGAGGCCCTTTAACCCTAAGCCACCGGTGCACGTGTGAATGTTGCACTTCGTTTGAGGTCGCCTACGAAACCATGGTTAGAATGAGAAGAGAAACGAGTAGTGTTGGGCTTACTCGGTGCTGAGGTAGGGG',
 'AGCCTTATGTGTCACGCTCACGTCTAAAATTTCCGTACAGGCATCGCGCTGTACTGAAACCTTGCAGACAACAACACACTAGTTTATGCAGGCTGCCGCAGCTCTCCATGGTGGCGCTTGACTAGTGCCATGTCACGGGGTATGGTCCTCCTATGT',
 'AACAGGCAGTCTCTTTCTCTATCGTCAAGCTGAAGACCAACTGGGCCGCCAGTTACAAGATTATGCGGGCTCGTGTTGACATTAGGATAGTAACCGTTCTACGTAGGTCAAAGCCTAGTTGCAATCACATTCGCCGCACCCTATCACTGTGGGCGA',
 'TCATATGTCTCCGCGATCGCCGGCTTGTCGTCCTACGAGAAAAACTTTGAAAACCTAGTTGGCAGAACGCAGAGGTATCACATGGGTACCACCAAGAATAGACGGTTCAGAAAGCGGCAGAAGGGGTCTTATGCAATGTAAAACGGCGAGGTTTTG',
 'TAATGCGGAAGTCTGGTTTACCAACGTAGTCGTCGCCCGACCGCCGCGCCGTTACTACAAGAGAGTCTAGTTATAACGCTGACTTTCACCGAACTATCTGGCTGTGGAATGTCTCAGCCTTCCTCGGGAGTGTCTCGACCTTTAGGGAACGACTAA',
 'GCCCGCGTGCCAACCAGAGGTTTGTCGCGGTCCTCGCAATTGAGGAATACTGCATACTACTACATACTAGTTTGCTCTCCTCTAGGAAAGCTTGGCTTTCAGGTGCGCTAAGCGCCACTGAGTCATTAACCCAGCCACGAAGACCCTTTAAGGTGG',
 'AGACACTATTATCTGCACTCAAGTTGGAAATACTGGGCTCGAACTAACAAAACTCTAGTTTGTCACGTGTTTTCGGCCTAGCGGATCATTCAAGATTACGAATCAAGTATCTGATATGCCCAATGGTTATCTCAGTCAAAACATGACGAGTGGCAG',
 'TCACTAGAGGCTTAGATTCTAGTTACTTAGTAGTAAAGTGAGCCTATAGGGTAAGGTCCAGCCGCATTGTTTGTTCGAACCGAGCACGCTCTCGAGTATGGGTGGACGTCCTCACGGTCGCTGTGTCGAACAAAAGCAGGTAATCTCGATTTCCCA',
 'GGACTTGGAACTTTTACTGCGGCAATGTAATATATGTGGATTTCGGATAGTATGACAAAGGTACATTCGATAGAAAAGCTAGTTGACTAAAACATTTTGGCAGCGATATAGGATCTTAGATTCGCTTTAACTATATAACAAGCGTATCTTTTATGT',
 'GCTTGGGAGGGCTATATCTAAAAACCCATATGGTTTCATATACTAGTTGAGTATCCTAAGCGGTTTTTCTCATCATGGTTCTTTTGGCCCTTAAACCGTGGTGGTCAGTACTCGGATAAATCTTACGAATGACCGCGCTCCCAGATAACGCGAGTG',
 'AAGACGCTAGTTGGATTAAGCACACTTCTTGGAGCTCGCCCAACGCGATCAGAGCGTTCCCTAGCGTACGAGGGGACAGACGCTGCCGGCGAGTCGGAATTGTCCCCTTTGTAAAACCCAATCGACCTCTACACGGCGCCTAAGAGCGTTCCTACG',
 'GAGCATACTGTTGCCCCACTAGAACCTTGATAACCGCCACAAATTATGGGCGGCTGAGTGCGAGTCTGGGGAAGTAACCAGCCTCCACATCGCTGTGACAGTCTAGTTGAGCGTCTGTGGCCATATTGGTTGGTCTGGTGCGCGAGAGACACCGTG',
 'GAGTGGACACAAACTTCCGTTCCAACGACGAGATGGGTGTCTGGCACCTAGATGCTAGTTGTAGCTCTACGGATTACCCGTCCTACTAGCTCTTAATTAATTCCAGTGCGGCCACTGATGACACTAGGATGTACGCGGTTGACGGTACAAGCTATA',
 'TGTGATACTTCCATTGTACAGAATAAGTAGTTAGTTTGAGCCACTTCGTGGAGTTACGTAAAAAGCCTAGTTACTTAGGGCGGCAGAGTCCAGAGTATCATCGTTGATATATTAGAACGTGTTACGATTCCTCCGATGAGACAGGTCTTGGAACGG',
 'CTCCCGTCGAGGTCCGCGAGGTCGGGATCGAGCAAGCAACTGCCTAACACAATAGTATCGACACGAGATAACTGTTCCCGTAAAATCCCAACAGAACGAATCCCCTGGCGCTCAGCGTCCGTCGACGACCTGAGTCTAGTCCTAAACAACCTAGTT',
 'CCGGTCACGAGTAGCGATAAAGCGTTAATCTGGTAATCATCTTCAAGTACTCACCGTCCTGACATTCTAGTTTGCGGCAAGTCAGCTACCCGAAGGAAATGTACTATGGGAAGATGCCGATATTCCGGTCGAAACACCGCCAGTACTATCACTCAA',
 'TGCGGACATGGTAGCCCACAGAAACAACCTGACCTGCCTAGCTTGACTGGAAATTGTCATCAAAACCTAGTTGTTCGATATGGATTCAACAAGTAGGAGGTGAGAGGAGACCCGAAAAGCCCGGTTCCAACGTCGTTGGACAAACAATCGGTGCGT',
 'CACACCACTGTACGTTATCCCTGCCTGTTCACATAGGCGTATACCCGATCGAGCCACTCAGTCGCTCCCGGGCGAATGGTCGGTGATGACCACACCGTCGATGCAGTAAGTGGAACTGACCGAAGGTGGAGGAAGAAACTAGTTGGAAGTCGCATG',
 'TAGTACTCCGTAGTCGACGAGCCGGAGACTCTAGTTTTCCAACTAGAATTCCTGCGGAGTCCACAATCATTTGGTACCAGAATCTAGAGCAGACCTCCATTAAAATCTAACGGATCGCATATGAATGGTTTGCGTACCGCTAACAACCATATAGGG']

res = GreedyMotifSearch(text, 12, 25)
print(*res, sep=" ")

text = ["GGCGTTCAGGCA",
"AAGAATCAGTCA",
"CAAGGAGTTCGC",
"CACGTCAATCAC",
"CAATAATATTCG"]

GreedyMotifSearch(text, 3, 5)


#Find the best motif according to a profile and a list of dna strings
#profile is a dictionary with keys A C G T
#dna is a list of strings
#output is a list of the best motifs for each string
def BestProfileMotif(profile, dna):
    kmer=len(profile['A'])
    pattern=dna.copy()
    for j in range(len(dna)):
        score=-1
        for i in range(len(dna[j])-kmer+1):
            score2=Pr(dna[j][i:i+kmer], profile) #probability of occurrence, want this to be as close to 1 as possible
            if score2>score:
                score=score2
                pattern[j]=dna[j][i:i+kmer]
            if score2==score:
                pattern[j]=pattern[j]
    return pattern


DNA = ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"]
test = ["TGA", "GTT", "GAA", "TGT"]
profile = ProfileWithPseudocounts(test)
BestProfileMotif(profile, DNA)


#need random number generator
import random

#Make random motifs of length k
#dna is a list of strings
#k is how long the desired motifs will be (integer)
#t is how many strings are in list dna (instructed to include this, I feel like you could easily leave it out and include it in the function)
#output is a list of strings length k
def RandomMotifs(dna, k, t):
    sub = []
    for i in range(t-1):
        x = random.randint(1, (len(dna[i])-k))
        sub.append(dna[i][x:x+k])
    return sub

#Start with random motifs to get closer to the best motifs
#Randomized way is faster, but not really as reliable...
#dna is a list of strings
#k is length of motif (integer)
#t is how many strings are in the list dna
#output is a string of length k
def RandomizedMotifSearch(dna,k,t):
    m = RandomMotifs(dna, k, t) #generate random motifs firt
    best_motifs = m
    while True:
        profile = ProfileWithPseudocounts(m)
        m = BestProfileMotif(profile, dna) #identify the best motifs out of the random motifs according to profile
        if Score(m) < Score(best_motifs): 
            best_motifs = m
        else:
            return best_motifs 

sample = open(folder+"/dataset_161_5.txt", 'r').read().splitlines()
text = sample[1:len(sample)-1]
text[0]
len(text)
sample[0]

res = RandomizedMotifSearch(text, 15, 20)

print(*res, " ")

prac = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

RandomizedMotifSearch(prac,8,5)


#Randomized motif search gets better over many iterations
N=1000
score1=5000
for i in range(0,N):
    if i%10==0:
        print(i) #This is just so I know where I am over the iterations
    tempmotif=RandomizedMotifSearch(text,15, 20)
    score=Score(tempmotif)
    if score<score1:
        score1=score
        best_motifs=tempmotif


print(*best_motifs, "\n")


p1 = (600-15)/(600-15+1) #probability of not capturing the 15-mer in one string
p2 = p1**10 #10 cases of not capturing the 15-mer
p3 = 1 - p2 #subtract to find the probability of capturing it once
p3

p4 = 1/586 * ((585/586) ** 9)
p3-p4*10


random.uniform(0,1)

#Given the probability of specific sequences occuring, randomly select one according to its probability
#input is a dictionary with a key and its probabiility of occuring (may not be in probability format yet)
#output is one of the keys according to the chance of its probability
def Random(prob_prof):
    new_prof = {}
    total = 0
    for key in prob_prof:
        total += prob_prof[key] #add all the occurences together
    for key in prob_prof:
        new_prof[key] = prob_prof[key]/total #set the new profile keys as their previous values over all counts of occurrences to get a probability value
    sum_prof = {}
    y = 0 #add them gradually together so you can look for ones that are in a range
    for key in new_prof: #need to add them together so when a number is randomly selected, there will be defined ranges
        y += new_prof[key] 
        sum_prof[key] = y
    randomint = random.uniform(0,1)
    for key in sum_prof:
        if randomint <= sum_prof[key]: #if random integer is less than a key value, it's going to be that key
            return key

test = {"seq1": .6,
 "seq2": .3,
 "seq3": .1}

Random(test)

#Generate a sequence of k-length according to the probability profile of it occuring in the dna string
#text is the dna string
#profile is the list of probabilities for each position for keys A C G T
#k is length of the string of interest
def ProfileGeneratedString(text, profile, k):
    n=len(text)
    probabilities={}
    for i in range(0, n-k+1):
        probabilities[text[i:i+k]]=Pr(text[i:i+k], profile)
    return Random(probabilities)

test = "AAACCCAAACCC"
prof = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
k = 2

ProfileGeneratedString(test,prof,k)

#A way of random sampling that will help you get good motifs. However, it can get stuck at motifs that are only locally good, and randomness can be an issue.
#Can also accidentally throw out good strings. Overall, I feel like this doesn't work super well, at least not how I wrote it.
#It is a bit faster than iterating through everything, though.
#dna is a list of strings
#k is the length of the motif of interest (integer)
#t is how many strings in dna (integer)
#N is how many iterations to run for Gibbs Sampler
#output is a list of the "best performing" motif strings, all length k
def GibbsSampler(dna, k, t, N):
    score1 = 5000
    bestmotif = []
    for i in range(0,70): #set randomized motif search for a few iterations to narrow down on well-performing motifs
        tempmotif=RandomizedMotifSearch(dna,k,t)
        score=Score(tempmotif)
        if score<score1:
            score1=score
            bestmotif=tempmotif 
    for j in range(0, N): #how many iterations we're doing
        i = random.randint(0, t-1) #going to choose a random string to drop
        reducedmotif=[]
        for b in range(0,t):
            if b!= i: #going to drop the randomly selected motif from analysis
                reducedmotif.append(bestmotif[b])
        prof = ProfileWithPseudocounts(reducedmotif) #generate a profile with these motifs
        mi = ProfileGeneratedString(dna[i], prof, k) #make a new string in place of the old one based on this profile
        reducedmotif.insert(i, mi)
        if Score(reducedmotif) < Score(bestmotif): #see if this new set of strings performs better
            bestmotif=reducedmotif #if it's better, use it for the next round
    return bestmotif

#Gibbs Sampler Problem
sample = open(folder+"/dataset_163_4.txt", 'r').read().splitlines()
sample = open(folder+"/gibbs.txt", 'r').read().splitlines()
sample[0]
text = sample[1:len(sample)-1]
text[len(text)-1]

text =  ['CTCTTGATTAAGGAGATGTAAAACTCTTTCCGGACATTAACTTGTCGATTGGTTCGTTTTATGATTGTTAGCCCATACAACGAGTGCTACTTTCGACGATTACCTGGCAACAATAGACAAGTCAGGGCCGCGGAAGACTGATCCCCTATACAGACCGTTATCATGCTACGAGAACGGTTGTCTAGCAACTCTTAGCTACGTGTGACGTCCACCGGCGTCGAGCCTGGCGACTATTAAATTCGCATGCGCTAAAAGCACCTGTTATAAACGGCTGTCAGCGATGTTCGGCCGATATGCGCATCTTCGTTTCCTCTTGATTAAGGAG',
 'ATGTAAAACTCTTTCCGGACATTAACTTGTCGATTGGTTCGTTTTATGATTGTTAGCCCATACAACGAGTGCTACTTTCGACGATTACCTGGCAACAATAGACAAGTCAGGGCCGCGGAAGACTGATCCCCTATACAGACCGTTATCATGCTACGAGAACGGTTGTCTAGCAACTCTTAGCTACGTGTGACGTCCACCGGCGTCGAGCCTGGCGACTATTAAATTCGCATGCGCTAAAAGCACCTGTTATAAACGGCTGTCAGCGATGTTCGGCCGATATGCGCATCTTCGTAAGCGCACCGGGGTGTTCCTCTTGATTAAGGAG',
 'GAGATGATAGGTTGGCCGGTTCGCCTCGATACGGTCCACGCCTGCTGGAATCTAGCTAGACAATTGCTTAGTGGATTCATTCTCCTCACCCCTGTAATTTACCCTTACCGGGGTGGGGAGGAAATACTCCACGTAGAACACGTTTACGAGCCTAAGGGCCGAGAATCACATAAGGCGTCTAACTATTAAGTGCCTTTGGTATCGATTATTGTGTTTTTCCCCATGCCCGCAGTCCTCCACTTAATAGACTGCTATCAACTATGGTAAATCAATTTCCACGATCGGGCTCTCGAACTTCTGTGTTATCCGATACGTCGCCGAAATC',
 'GCCTAATTGAATTATAAAGTATTTCGTCCGACATATCGCCATGTTGACTGTATGCGCATGGAATTCGCTTCGAGAAGTTCCTCGGGGTGAGGCACGTTTTGAAGAACCCGGAAGCTCCTTCGGTTGAGCCTAAGTTTACTCTATAGGCAATCTCACCATCCGCGTCCACCCAATCGCGTGAGGTAAGATCTAAGTCCGGCTGCAAGTATCCATAAGGCCCCTTGCGGATGGTCACGTCTCTTAGCAAGGAGTCAATGAGATCGGCCCTCCCTACCCTTAGTCTATGTTTTGGCATAAGCATTGGGAATTGTGTAGGATATGTGAG',
 'CGTTTCATCTACATGACATTGCTGCTACGACATGCGTGTCGCCCTCCTGGAGCCCAGTGTTGATCACCGTGGGAACGTTCCTAATAGCTGAAGTGAGGACTGGGAATTCGTTCACTTGACGTCTCACCTGTCGATTTATGCATTTGAAGCTCAATTTGGGGGTAAATTGGAATGAGAGCGAAGAGACGTTTACCTATCCTTCTAATAGGAAACTTCTAGTTGGATGATGAGATAAGTTTTATGGGGTGTATATTGGGCGTCAATGAACCCTCGCCAGTGTAAACACCAATTTCCATTGAGGTTGGGTGGTAGAGTCCGCGGGACA',
 'TAGACTAACCCACACGTAACCAATTGGTTTTTCGGACAGGGTGAAGGGATGTGTGCATCGAAAGTTTTTAGCTACGACTGTAATATCCACTTCACCTCTGTCCACCAGTACAATCCAGGTAATAAATCTCCTCTGGCTGGTGCTTTAAAGGGAGTCTGTTTCACGATCCTTGAACAGGTGCGTCTCACGAGGACGTGTATGAATTTTCATAATAGACGTGTTCCCGAGCCACCAACAGGAGCGTGCCTGATTCGGAAGATGCAAAGCCAATTGCATACCACCTGCACAGGAGGAGGCATGGATCGCAAGTTTACCGGGTGCAAGG',
 'CCTACTTGACAAGCGTAGGCGCGGTACGCAAGTGTTGCGTTCTCCCTCGCAACACCCGTCAGTGCTACGGGGACGGGTTTTACGACTTGACGCTCTTCCGGCCACCTGCATTAACTCGACGGAATGAGCACGGCTCGGTAGGCGATCGAGTATGCGTCATGGGAAAATAGGAATCGGACGCCCCTCGGGCATATTAAGCCTGCGTTCGTGTTGTCCTTACGATATTAGCCTACCAAGTTTCGAGGGGTGCCAAGCTCAAGTGATCCGGAACTTTGCTTTACCACCACCGCCATCCAGGGCATTATACATCGCTCCCTTGTGACCT',
 'TAATACACATCCTCGGACTCCACATGACGATACCACTAAAAAATCAACGACCTTTCGGCCGCATGATAGGTCATGAGGGGGCAGTTTATTCTCGGTTCCTGTTTACCGGGGTATGGTAAATCTGCAGGGTTGCACACCCGATCAGCTTGTAGGCTTTCGTGCTTTCAGATTTCTAACAATACGTTAAAGATTTTTGAGTTAGAGAAAGAGCGTCGAACATACTGTCGTACCAATTTACTCTTTACGATCATTCGCCCGCAGCATTCCGGTGCAATCGATTATTCGCATAGTCATTCCCCTGTTCCGTGGCTATTCTTCGTACCTT',
 'AATGGGATTGCTGAACAAGAAGGCGGCTTAGACTGTCTATGGCTTCCGATCGGACTAACGGCGAATAATAGTAAGATTACGGATCCCTGACAGCTTCAGTCCGCAAACGACACCACAGGCTCCTGTAGTAAAACAGACAGCCACTATAGCGCGATTGTTGGCCCCCCCTTAAGTTGCTCGGGGTGGTCCAACAGTCCCCAGAAGACATACGACGGGATGTATATAATGAAATTCGCCTTCTTTAAGAAGATGCTCTGGCAGTTTCATATAGGGGCCCGCTGTTGAAAATCGGATGAGTGAGGATACATGCGTTTGCGTTCGTGTC',
 'GATACTCCTATCGCGCAGTGACCTCCCTGCGTTCATATTTAGCCCTACTTTGACGAGACAGATAGCTGGGAAAGCCTATTCGACATATATACTGCGATGACTCCGGAACGTAAAAGAGTAAATCGACATATTTAGTGGCTTGGATTTGAGTAGTATCGCAACCTACGCCGATGCGGAAAATTAAACATACCGGGGTGTCCCATATGAGGGGGGCGAAATCTCCGAGGATTGAGTACTCGTGCCCCCGACTTTTTTTCGACTCGCGGCAATGAAAACCGAAGGAGGCACGAAGTGGTACATGTGTACCCCTCTTTGGTTACTCATG',
 'CGCAGGCTCATTCGTTCACGAACACACGGAACTACCCAGCGCGTTGATGCTCCAAAACGAGGCCACGTTCACAGAACCGAAACACCGATAAAAGCGCGCCAACAACCCGACGACGCACAGGGTGAAATGGCACTTACGGCTCTTTCATGATCTTCGACCGAAGGAATGGAGGGGGTCACCTGGCCCGGCCCGGTGAGTGCTTGTATAGGCGTTTGTACTGAGGTACCAGGACCGGGCGCTGCACAAGCTGCCATTCTAGCGTATTCTCATATCCAAATGGCTCGCAAGTTTAGGAGGGTGGGGCTCCCGCCAGGCCGTCATATCC',
 'ACGTTTCGCAGCTGAGGTAAGGAAACCGGGGTGGAATCACCCTCGAAGCTGGTCGCGCCGGCATCTATTGTTGAGCAGGTCATCACAATTCCTCTATTTCTATGATACAACTTCGACGATCCACGGGATATGTAACGCCGGAACACAGGAGTAAATGTGATTGACAGGGGCTCATCCGTCTGCCCAAACGGCATCTACGCAATGACTGCATAGGTTTTGTGTAAAAGAGTTTGTCATCTACCCAACCAGGACAAGTCAGCCCGCGCAAACGGCCCACGCGCACATATCAAGCCCGTCAGGCGCCCGCAGAAACAGATCCTAAGTT',
 'GTGGCTGTGCGTAACCGTCTAAATGTAAAAGCGCACATGAGGTAAGTTTACACAGGTGACCCAAGTGATCCTGATCGAGATGGGTAACCGCATTTCTGTGAGTCGGGACACTGGGTGTTACCAGTTGCCAGAAATTCGGCGGGGAGTGAGTTCGGTCGGTATTTATGACTAGGTCATTGGGCTGCAGCGCTCCGCAACAGTCCATGGTTTATAGTTGGAACAGACCGGGGTGATTCATTAAAAGAACATTCATCTGCTTAGAAAAATAGATTTACGTTCCGTAGAACCGTAAGAAATTACTGGCTAACCCAACATAAAAGCTGAG',
 'TCGAATCCGCCACATGCAAGGCTCAATGTTGACAACTCTTGTGGAGAAGACATTGCAAGACAGCTTGAAGGAGGTCCGCTAGAGCTAGTCTACGCTCCGTGTCAAAGCCTGGAGAACATACGATAATGAGTTAGACCGGGGCCAATTAGTTTACCGGGGATCGCTGAACAACCGGTCCGTGACCATACACTTAGTTGGGTAGCAATACATCTGGCCCGGTCAATTTCATCTAAGGCACCCGATATGAGGACGTGTGCAATACACATATTTTCGGTGCTGTCATGTCCTGTGAGGTTTGCATGGCTGACCGTACTAGTATTAACAG',
 'GGCCTTAATGAATCGCTCTGTCATGCATGCATTGGGATGGGGACCCCTCCGTTAGCTGTGATGGGTCGAGACGTACGATGTACCGCCCCTTTTACCGGGGTGAGCGATTCTCGTGCGAAATGTTCTCCCACTTTGTCCGGCCGTGCGCGCAGCATACTGGGCAGCCTGCGTTCCCCGCCCCCCCACATGACCACGACTGGGTTCGCCATCGTCAGCTTAAAATCCGTATGGTTAGGGAAGATAGCGTCCAAATGGGAAGCATGCACGTAATTCAGACTGAGTCCCTCTGTATTCTGTCTTGGACGTAATGAAACTCTATAAAACT',
 'CCCTAGCAAAAGCCCTCTTCAATCACTTGCAATTGTTCTTTACCCCCTTAACTCAGCTTGACCATCATCAATCCAAACCGAAGCTTCGGCCACATCCAGTATGCCGAACAAAGGCAGTAGATTATGCGATCATTCTGTTCTATAACTTTCTTTCTACCCTCACCGATCCACATATTAGCTGTGATTTCGAGTCATTCTGATTGACTATTCATGACTTCCGGCTAAAGACAGGTACATGAAGTGAGCCGGGGTGATACAGGAGTGGGATGCTTGGGCAGCGTTCAATTGGAGAAATCGGAGAGTTGCTACCATTCCTGCTGTCTGG',
 'GCGTGACGGCTCTATAAGAGAACTACGACCAATAGTACCGGTGGTCCTCAGCCTTAAATATAGTGTAAAGTCGTCCGGGGTGATTCAAATGGGTGTCTTTAAACTTATTTAGAAGTACATTGTGCCTAGGTTTCCGGGACTTGCCATAATTGAGAGTCCCTCATTCTCGGTGAGGAGCGCCGAAGTCCCGTTAATCTGGCGTGTCCCGTGATGCATCATCTAAGTTATCAGTCAGTCGCACGCACTCCTACATGACTGAACCAGTGCGCGCTGAGATGGTACGCGTGCTCACTGTCCAAGGAGACGGACACGTATCAACTGGCGC',
 'GCCTATGGAATTGCAATGGAGTTATGTCCAGTACAGAGGTGAAAGTTTACCGGACAGAGATTACCAACCCCGGGATTAGGGGAGATCGAGCTGCGGGCTCGTGGGCCAAGTATTCAACGAACAAAGCTTAAGTAAAGCAGCGAAACGCCTACCGGTACAACAGGCGGTTCAGGTGACTACCAATAAAGTAAATGTTCGGACGCAGACGTCTAAGCAAGTGACGGCCTAGGAGTTTACGCCCTACAACCCACCAGCCACCAACGGCAAATAAGTCCCTACTGACCGCGGCATTTTGCGCACCGAACTAGCCGTCAACCACATCACG',
 'CGGTCCATGTCTCTAGCGCAAATGGATAGGTTCTGTATATACGGCACCTGGCCCAGCACGTCTTTACACAATAAACAATAACCCGAGTGGTGTTAGGTGAGACTTACTAAGGGACCCGCGCAACAACGGGTCCAAGGTGACGGATTTTAATCGTTGCGTGTCGATATCTCGCAGCATCTAAGACTGAGAATGGCGGGATTCACTCCTCGGACTAGGACATCTTCCCAAAGTTTACCAATGTGAGAGACAGAGGTGCACCACTAGGCACGGATGTATGCGCGAGCAATTGAACAATACGGTCCTGTTATACAGTTCACGTTAACAC',
 'GTCGGTTCAGAGCAACTTTACATAGAGGAACGGAAAGAGCAACATTCTTCCCAAGTTTACCGTCATGGTTTGCGAGTACAGCGGCCGGCACTACTGGCGGAGTGAGCCACATCGTTGGCTGGGACCGAGAAACTGCGAGTCTTTAAACGGACCCGCGCCCCAGACACTAGTGTTTCCTATGCGCGCGCATAAAAAGCCAGTCCCGGTAACTGGAGTTCAGGACCAAGGAGTTTGGACAAGCTTGCTAATCGAAATACCATTTGTGTTGCGATCTTGGAGCGTGCGTAGCGCTTACGGTCGAAACGTACCCCGCAGTATTATACCC']

res = GibbsSampler(text, 15, 20, 2000)
Score(res)


print(*res, sep = "\n")
