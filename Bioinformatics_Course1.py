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

len(pattern)

def PatternCount(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if pattern == text[i:(i+len(pattern))]:
            count += 1
    return count

text = "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT"

PatternCount(text, "TTT")

PatternCount(text, pattern)

def FrequentWords(text, k):
    freq_patt = {}
    count = 0
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        counti = PatternCount(text, pattern)
        freq_patt[pattern] = counti
    pat = []
    m = max(freq_patt.values())
    for key in freq_patt:
        if freq_patt[key]==m:
            pat += [key]
    return pat

text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"

text[0:4]

FrequentWords(text, 4)

bob = {"george":4,
       "kevin":3,
       "craig":8}

bob.keys()
bob["george"] += 1

bob["larry"] = 1

sample = open(folder+'/dataset_2_10.txt', 'r').read().splitlines()
text = sample[0]
k = sample[1]

FrequentWords(text, 12)

text[0]

text[::-1]

def ReverseComplement(text1):
    output = ""
    text = text1[::-1]
    for i in range(len(text)):
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

def PatternMatching(pattern, text):
    pos = []
    for i in range(len(text)):
        if text[i:i+len(pattern)]==pattern:
            pos.append((i))
    print(*pos, sep=" ")


pattern = "ATAT"
text = "GATATATGCATATACTT"

text[0:4]

PatternMatching(pattern, text)

sample = open(folder+'/dataset_3_5.txt', 'r').read().splitlines()
pattern = sample[0]
text = sample[1]

#Find k-mers that are repeated at least t times in a l length clump

def Clump(genome, k, l, t):
    dic = {}
    sequences = {}
    for i in range(l-k+1):
        if genome[i:(i+k)] in dic:
            dic[genome[i:(i+k)]] += 1
#            if dic[genome[i:(i+k)]] >= t:
 #               sequences.append(genome[i:(i+k)])
        else:
            dic[genome[i:(i+k)]] = 1
        #print(dic)
    for i in range(l, len(genome)):
        for key in dic:
            if dic[key] >= t:
                sequences[key] = ""
                #print(sequences)
                #print(dic)
        dic[genome[i-l:(i-l+k)]] -= 1
        if genome[i-k+1:(i+1)] in dic:
            dic[genome[i-k+1:(i+1)]] += 1        
        else:
            dic[genome[i-k+1:(i+1)]] = 1
    print(*list(sequences.keys()), sep=" ")
            
text = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
5 50 4

text2 = "CCACGCGGTGTACGCTGCAAAAAGCCTTGCTGAATCAAATAAGGTTCCAGCACATCCTCAATGGTTTCACGTTCTTCGCCAATGGCTGCCGCCAGGTTATCCAGACCTACAGGTCCACCAAAGAACTTATCGATTACCGCCAGCAACAATTTGCGGTCCATATAATCGAAACCTTCAGCATCGACATTCAACATATCCAGCG"

Clump(text, 8, 28, 3)

BetterClumpFinding(text, 8, 28, 3)

sample = open(folder+'/dataset_4_5.txt', 'r').read().splitlines()
text = sample[0]
sample[1]
ecoli = open(folder+'/E_coli.txt', 'r').read().splitlines()
ecoli1=ecoli[0]
ecoli1[0:30]
Clump(ecoli1, 9, 500, 3)

PatternCount("ACTGTACGATGATGTGTGTCAAAG", "TGT")
FrequentWords("CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA",3)
PatternMatching("CGC","ATGACTTCGCTGTTACGCGC")

def PatternToNumber(pattern):
    pattern = pattern[::-1]
    q = len(pattern) - 1
    number = 0
    while q >= 0 :
        if pattern[q] == 'A':
            num = 0
        if pattern[q] == 'C':
            num = 1
        if pattern[q] == 'G':
            num = 2
        if pattern[q] == 'T':
            num = 3
        number += (num * (4**q))
        q -= 1
    return(number)

pattern[::-1]

912/(4)
228/4
57%4
14.25%4
3.5625%4
0.890625%4

pattern = "ATGCAA"
print(PatternToNumber(pattern))

def NumberToPattern(index, k):
    num = index
    rem = 0
    pattern = ""
    q = k - 1
    while q >= 0:
        rem = num%4
        num = num/4
        if rem < 1:
            pattern += "A"
        if rem <2 and rem >=1:
            pattern += "C"
        if rem <3 and rem >= 2:
            pattern += "G"
        if rem < 4 and rem >=3:
            pattern += "T"
        q -= 1
    pattern2 = pattern[::-1]
    return pattern2

5437%4

NumberToPattern(5437, 7)
NumberToPattern(5437, 8)


def ComputingFrequencies(text, k):
    freq = {}
    pattern = ""
    list = []
    for i in range(4**k):
        freq[i] = 0
    #   print(freq)
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        j = PatternToNumber(pattern)
        freq[j]+=1
    for i in range(4**k):
        list.append(freq[i])
    return freq

x = ComputingFrequencies("ACGCGGCTCTGAAA", 2)

x = [3, 4, 1, 3]

max(x)

max(x.values())

sample = open(folder+"/dataset_2994_5.txt", 'r').read().splitlines()
text = sample[0]
k = sample[1]
k
x = ComputingFrequencies(text, 5)
f = open(folder+"/answer.txt", 'w')
f.write()

def FasterFrequentWords(text, k):
    freq_patt = []
    freq_array = ComputingFrequencies(text, k)
    max_count = max(freq_array.values())
    #print(max_count)
    for i in range(4**k):
        if freq_array[i] == max_count:
            pattern = NumberToPattern(i, k)
            freq_patt.append(pattern)
    return freq_patt

FrequentWords(text, 4)
FasterFrequentWords(text, 4)

def LastSymbol(pattern):
    pattern2=pattern[len(pattern)-1]
    return pattern2

def Prefix(pattern):
    pattern2 = pattern[0:(len(pattern)-1)]
    return pattern2

def SymbolToNumber(symbol):
    if symbol == "A":
        return 0
    if symbol == "C":
        return 1
    if symbol == "G":
        return 2
    if symbol == "T":
        return 3


def PatternToNumber2(pattern):
    if pattern is not "ACGT":
        return 0
    symbol = LastSymbol(pattern)
    prefix = Prefix(pattern)
    return 4 * PatternToNumber2(Prefix) + SymbolToNumber(symbol)

PatternToNumber(pattern)

sample = open(folder+"/dataset_3010_2.txt", 'r').read().splitlines()
pattern = sample[0]

NumberToPattern(8308, 11)

sample = open(folder+"/dataset_3010_5.txt", 'r').read().splitlines()

z = {1: 3, 20: 2, 4: 2}

sorted(z)

max(z)

z[1]+=1
z[1]

z[2]+=1

max(z.values())

def FindingFrequentWordsBySorting(text, k):
    freq_patt = []
    count = {}
    for i in range(0, len(text)-k):
        pattern = text[i:i+k]
        Index = PatternToNumber(pattern)
        if Index in list(count.keys()):
            count[Index] += 1
        else:
            count[Index]=1
    SortedIndex = sorted(count)
#    print(count)
    #print(SortedIndex)
#    for i in range(1, len(SortedIndex)):
#        if SortedIndex[i] == SortedIndex[i-1]:
#            count[SortedIndex[i]] = count[SortedIndex[i-1]] +1 
#            print(count[SortedIndex[i]])
    max_count = max(count.values())
    #print(max_count)
    for key in count:
        if count[key]==max_count:
            pattern = NumberToPattern(key, k)
            freq_patt.append(pattern)
    return freq_patt
    
FindingFrequentWordsBySorting(text, 4)    

def BetterClumpFinding(genome, k, l, t):
    freq_patt = []
    clump = {}
    for i in range(0, 4**k-1):
        clump[i]=0
    text = genome[0:l]
    freq_array=ComputingFrequencies(text,k)
    for key in clump:
        if freq_array[key] >= t:
            clump[key] += 1
    for i in range(1, len(genome)-l):
        first_patt = genome[i-1:i-1+k]
        index = PatternToNumber(first_patt)
        freq_array[index] -= 1
        last_patt = genome[i+l-k:i+l]
        index = PatternToNumber(last_patt)
        freq_array[index] += 1
        if freq_array[index] >= t:
            clump[index] = 1
    for key in clump:
        if clump[key] ==1:
            pattern = NumberToPattern(key, k)
            freq_patt.append(pattern)
    return freq_patt
            

    
def Skew(genome):
    skew = 0
    G = 0
    C = 0
    result = [0]
    for i in range(len(genome)):
        if genome[i] == "G":
            G +=1
        if genome[i] == "C":
            C +=1
        skew = G - C
        result.append(skew)
    return result

result = Skew("GAGCCACCGCGATA")
print(*result, sep=" ")


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

def HammingDistance(p, q):
    dist = 0
    if len(p)==len(q):
        for i in range(len(p)):
            if p[i]==q[i]:
                score+=0
            else:
                dist+=1
    return dist

p="GGGCCGTTGGT"
q="GGACCGTTGAC"

sample=open(folder+"/dataset_9_3.txt", 'r').read().splitlines()
p = sample[0]
q=sample[1]

HammingDistance(p,q)

def ApproxPatternMatch(pattern, text, d):
    result = []
    n = len(pattern)
    for i in range(len(text)-n+1):
        x = HammingDistance(pattern, text[i:i+n])
        #print(x)
        if x <= d:
            result.append(i)
    x = len(result)
    return x

def ApproximatePatternMatching(text, pattern, d):
    result = []
    n = len(pattern)
    for i in range(len(text)-n+1):
        x = HammingDistance(pattern, text[i:i+n])
        #print(x)
        if x <= d:
            result.append(i)
    x = len(result)
    return result

ApproxPatternMatch("ATTCTGGA","CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3)

ApproxPatternMatch("CCC", "CATGCCATTCGCATTGTCCCAGTGA", 2)

pattern = "ATTCTGGA"
text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"

sample=open(folder+"/dataset_9_4.txt", 'r').read().splitlines()
pattern = sample[0]
text=sample[1]
sample[2]

x = ApproxPatternMatch(pattern, text, 6)
print(*x, sep = " ")

def Count2(text, pattern):
    x = len(ApproxPatternMatch(pattern, text, 2))
    return x

Count2("CATGCCATTCGCATTGTCCCAGTGA", "CCC")

ApproxPatternMatch("GAGG","TTTAGAGCCTTCAGAGG", 2)

sample = open(folder+"/dataset_9_6.txt", 'r').read().splitlines()
pattern= sample[0]
text=sample[1]
n=sample[2]

ApproxPatternMatch(pattern, text, 2)

n = {"a":1, "b":2, "c":3}
list(n.keys())

import itertools

bases=['A','T','G','C']
x = [''.join(p) for p in itertools.product(bases, repeat=4)]
len(x[0])

def FrequentMismatch(text, k, d):
    bases=['A','T','G','C']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    dic = {}        
    for j in range(len(kmers)):
        for i in range(len(text)-k+1):
            x = HammingDistance(kmers[j], text[i:i+k])

            if x <= d:
                if kmers[j] in list(dic.keys()):
                    dic[kmers[j]]+=1
                else:
                    dic[kmers[j]]=1
    keys = []
    for key in dic:
        if dic[key]==max(dic.values()):
            keys.append(key)
    return keys

text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
   4 1
   
FrequentMismatch(text, 4, 1)

sample = open(folder+"/dataset_9_7.txt", 'r').read().splitlines()
text=sample[0]
sample[1]
sample[2]


x = FrequentMismatch(text, 6, 3)
print(*x, sep=" ")

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

HammingDistance("CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG", "ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT")

sample=open(folder+"/dataset_9_8.txt", 'r').read().splitlines()
text=sample[0]
sample[1]
FrequentMismatchRC(text, 5, 2)
            
#    Neighbors(pattern, d)
#        if d = 0
#            return {pattern}
#        if |pattern| = 1 
#            return {A, C, G, T}
#        neighborhood ← an empty set
#        suffix_neighbors ← Neighbors(Suffix(pattern), d)
#        for each string text from suffix_neighbors
#            if HammingDistance(Suffix(pattern), text) < d
#                for each nucleotide x
#                    add x • text to neighborhood
#            else
#                add FirstSymbol(pattern) • text to neighborhood
#        return neighborhood
def Neighbors(pattern, d):
    if d == 0:
        return pattern
    elif len(pattern) == 1: 
        return {'A', 'C', 'G', 'T'}
    neighborhood = []
    suffix_neighbors = Neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        nucleotides = {'A', 'C', 'G', 'T'}
        if HammingDistance(pattern[1:], text) < d:
            for x in nucleotides:
                neighborhood.append(x + text)
        else:
            neighborhood.append(pattern[0] + text)
    return neighborhood

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

def MotifEnumeration(dna, k, d):
    motifs = []
    for i in range(len(dna)):
        if i == 0:
            for j in range(len(dna[i])-k+1):
                pattern = (dna[i][j:j+k])
                pattern_p=Neighbors(pattern, d)
                if isinstance(pattern_p, str)==True:  #If there's only one value, it's a string not a list.
                    motifs.append(pattern_p)
                else:
                    motifs.extend(pattern_p)
    motifs = list(set(motifs))
    temp = {}
    removemotif = []
    for j in range(1, len(dna)):
        for motif in motifs:
            #print(motif)
            for i in range(len(dna[j])-k+1):
                if HammingDistance(motif, dna[j][i:i+k]) <= d:
                    if motif in list(temp.keys()):
                        temp[motif]+=1
                    else:
                        temp[motif]=1
                else:
                    if motif not in list(temp.keys()):
                        temp[motif]=0
            if temp[motif] == 0:
                removemotif.append(motif)
        #print(motifs)
        temp = {}
    for motif in list(set(removemotif)):
        motifs.remove(motif)
    return motifs

test = ["ACGT", "ACGT", "ACGT"]

Neighbors("ACG", 0)

MotifEnumeration(test, 3, 0)

#cycle through the strand
#when moving to the next strand, ask if motif value is less than 0. If so, then remove from motifs list
#then move to the next strand

3 1
text = ["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]
MotifEnumeration(text, 3, 1)

text = open(folder+"/dataset_156_8.txt", 'r').read().splitlines()
thing = text[1:len(text)-1]
text[0]

MotifEnumeration(thing, 5, 1)

{'TTG': 1, 'ACT': 1, 'ATA': 1, 'GGT': 0, 'TTC': 2, 'AGT': 0, 'TGA': 2, 'ATC': 0, 'GTT': 1, 'TTT': 2, 'TGC': 1, 'ATG': 0, 'TCT': 1, 'CGG': 0, 'AGC': 1, 'AAT': 0, 'CGC': 1, 'TTA': 1, 'GGG': 0, 'GCC': 1, 'CTT': 2, 'TGG': 1, 'TAT': 0, 'AGG': 0, 'ATT': 1}

10/4
1 2 3 4 5 6 7 8 9 10 
A A A G C T G C T G  
7*15

import math

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
    # your code here
    return count


def Profile(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = {}
    profile = Count(motifs)
    sub= Count(motifs)
    for x in range(t): #for key in position x...
        for y in range(k): #for list position in that key...
            symbol=motifs[x][y]
            profile[symbol][y]=(sub[symbol][y])/t
    # insert your code here
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

test = Profile(motifs)

test["A"][1]

math.log2(test["A"][1])

len(test["A"])

def entropy(motifs):
    prof = Profile(motifs)
    result = []
    result2 = 0
    for i in range(len(prof["A"])):
        temp = 0
        for key in prof:
            if prof[key][i] != 0:
                temp += abs(prof[key][i]*math.log2(prof[key][i]))
            else:
                temp += 0
        result.append(temp)
        result2 += temp
    return result2

0.5*math.log2(0.5)*2

0.25*math.log2(0.25)*4

1*math.log2(1)

0.25*math.log2(0.25)*2+0.5*math.log2(0.5)


entropy(motifs)

def HammingDistance(p, q):
    dist = 0
    if len(p)==len(q):
        for i in range(len(p)):
            if p[i]!=q[i]:
                dist+=1
            else:
                dist+=0
    return dist

dog = {'blue': 5, 'brown': 1, 'zebra': 4}

min(dog)

min(dog.values())

len({})

dog = "meow"

cat = []

cat.append(dog)

def Motifs(pattern, dna):
    motifs = []
    print(pattern)
    k= len(pattern)
    if isinstance(dna, str)==True: #verify if the input is a list or a string
        new_dna = []
        new_dna.append(dna)
        dna = new_dna
    for seq in dna:
        score = k
        temp_motif = {}
        for i in range(len(seq)-k+1):
            new_score = HammingDistance(pattern, seq[i:i+k])
            if new_score < score:
                temp_motif[seq[i:i+k]]=new_score
        if len(temp_motif)==0:
            motifs.append(seq[0:k])
        else:
            low_score = min(temp_motif.values())
            for key in temp_motif:
                if temp_motif[key] == low_score:
                    motifs.append(key)
                    break
    return motifs

dna = ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG"]    

motifs("AAA", dna)

motifs("GATTCTCA", "GCAAAGACGCTGACCAA")

motifs("AAG", "GCAATCCTCAGC")

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

def MedianString(dna, k):
    kmers = []
    median = ""
    if isinstance(dna, str)==True:
        distance = len(dna)
        for i in range(len(dna)-k+1):
            kmers.append(dna[i:i+k])
    else:
        distance = len(dna[0])
        for seq in dna:
            for i in range(len(seq)-k+1):
                kmers.append(seq[i:i+k])
    kmers = list(dict.fromkeys(kmers)) #remove duplicates
    #print(kmers)
    for kmer in kmers:
        if distance > d(kmer, dna):
            distance = d(kmer, dna)
            median = kmer
            print(kmer)
            print(distance)
    return median

MedianString(text, 7)
MedianString2(text, 7)

#Search kmers not in the main text
def MedianString2(dna, k):
    kmers = []
    median = []
    if isinstance(dna, str)==True:
        distance = len(dna)*10
        for i in range(len(dna)-k+1):
            kmers.append(dna[i:i+k])
    else:
        distance = len(dna[0])*len(dna)
        for seq in dna:
            for i in range(len(seq)-k+1):
                kmers.append(seq[i:i+k])
    kmers = list(dict.fromkeys(kmers)) #remove duplicates
    kmers2 = []
    for i in range(len(kmers)):#add more kmers that are not just in the main text
        kmers2.append(kmers[i])
        newkmers = Neighbors(kmers[i], k%2+1) #arbitrarily decided on this, tbh, but the length will likely change
        kmers2.extend(newkmers)
    kmers2 = list(dict.fromkeys(kmers2))
    #print(kmers)
    #print(kmers2)
    for kmer in kmers2:
        if distance > d(kmer, dna):
            distance = d(kmer, dna)
            median = []
            median.append(kmer)
            #print(kmer)
            #print(distance)
        elif distance == d(kmer, dna):
            median.append(kmer)
    return median

test = ["ATA","ACA","AGA","AAT","AAC"]

d("AAA", test)

MedianString2(test, 3)
motifs("GTC", test)
motifs("AAA", test)
d("GTC", test)
d("GAT", test)

motifs("GAT", test)

d("AAA", test)

MedianString2(dna, 3)

3
dna = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"]

motifs("GAC", dna)

d("GAC", dna)

MedianString(dna, 3)

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

def ProfileMostProbableKmer(text, k, profile):
    score=-1
    pattern=text[0:k]
    for i in range(len(text)-k+1):
        score2=Pr(text[i:i+k], profile)
        if score2>score:
            score=score2
            pattern=text[i:i+k]
        if score2==score:
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
    # your code here
    return count

def ProfileWithPseudocounts(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = {}
    profile = CountWithPseudocounts(motifs)
    sub= CountWithPseudocounts(motifs)
    for symbol in "ACGT":
        for y in range(k):
            profile[symbol][y]=float(sub[symbol][y])/(t+4)
    # insert your code here
    return profile

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
        if Score(motifs) < Score(best_motifs):
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

sample = open(folder+"/dataset_5164_1.txt", 'r').read().splitlines()
pattern = sample[0]
text = sample[1].split(" ")

d(pattern, text)

text=["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
      "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
      "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"]

MedianString2(text, 7)

text = {"A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
        "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
        "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
        "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]}

Pr("CAGTGA", text)

test = ["ATA",
"ACA",
"AGA",
"AAT",
"AAC"]

MedianString2(test, 3)

AAA

d("AAA", test)
test = ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"]

motifs("AAA", test)

test=["TCTGCAGAGCGCAAGCGTCTCAATGTTTTTCGCTGCAAGTGTTAGCTGCGTTTGTCGACACACAAATGAAGTGACCACATCAAACCTAATTATG",
      "TCCCGCTCGAGGACAGACACCGGTGCAGCCGAGGGTATTATAGTCTGCTCGTATGGTTTGTATGGAGGATAATAGATGGGGTGACGAAATGAGA", 
      "CCTCCCTGCCAGGTTGGTGAATTTAAGCATAAAGACCTGGAGGACTGACGGGTCTGGTCGACCACCATTCGTCGGTTCTCGAGCGCTGTTTCTA", 
      "GATACAGGTCGCTTGAAATGTCCCTGAGACGTCCGGGCCAAGAGCGAAACAAATCTCAGTCCTGGGTGGCGGCAGCGAAGGGATTAGATTACTT", 
      "ATGGAACCGTGTCACGGCAGTCTACCTTTTACCCAGAGCTGTATATGACGGTTAGTCCCCGAGGCCATCTCGCACCCTACTGAGCACAAATATA", 
      "AATCCCTGCAGGTCAGGCTGTTCTGAAAGAACCCTGCAGGCTGGTGCTATTCCTTAAGACGCCGAGGTTCAGATACTCTCAGCCAGAGAGCAGA", 
      "CCTATTAACTTCTCTCTAGATCGTAAACGACTATGGACTTACACGTCCCGCTTTCTTGTCTTGGGTGCCGTCACACTAGCCCACTACAGTCGGG", 
      "GGCGCACCTTCGATATTAATGCGAGTCAATTGCCGTTTAAGTCCGCGCGTAGCGGATGCACGTAGGGACTTTAAAGCCCCATGCCAACTGTATA", 
      "TCTTGGTTACAATTACCCTCCTGAGCCGATAAGCGCCTCTCAGGCACACTGTGGCAGTTCATCGTGTCCTTCGTCTTATTCCCTGCCCCCAGCA", 
      "ATCTCAATGCGATATATATGGTCCTGCCTAAAAATAGAGCTAGTGAGCGTGGTGGTATCATTGTTGACTACCCGGGACAGGCGACACCGTACGA", 
      "ATCTATAGAAACATCTAGCACAGGAACACCCATAAAAACACCACTAGCCGACCTGTTCGTCGGTTTTCCATGCGAATGGTAATATTACGTTTGT", 
      "ACGCATACGTCAGAAATCTGTGAAAGTACGACGCCCACTGGACCCACGGTCGCTGCGTTCACAATCATTGTGATCCTGGATCCATAACGCCTAA", 
      "GACCAAGGCCCGAAAGTGTGGCGAGATAAGAAAACCCACTCCATATGACCGTGGATCCCGTGGGAGTCCTCCTTTCATCCCTCTAGCGCTAGAA", 
      "ATACCCCCGTAGTTTTCGGACTCAAAAAGTATTTATGTGACTTGCAATAGTAGAGACCGGCGTCTTTGCCGCTGTCACCTGCTGTATGATCACC", 
      "CGATGGTACCCAAGAGTGGTCGGTCTACCAAGCATCGCGACCCGCTTTCGTTCACCAGCCGAGAGTGCCATGTAAGTCAGCAGCTAACACCCTG", 
      "TGTCTCAAGCATTACTATGGCTGGCCTCGTTATTAGATTCTGCCTAAGCCGCATGGGCGCGTGACCGTGGGCACATTGCCCTGAATTACCCGCC", 
      "AGGTGCAATTAACTGTACTAAGCTCAAGATCTCAATCCGCAAAGGCCATTCAGGCGACTCCGTGTGGTGCCTATACCGCACAAACGTAGGAGCG", 
      "CATTAAGACGGATTGGGACCCGCGGGAGAAAGGCGCGTTGGAGCGAGGGGGTCATTTGCACTGCCGAGCAGACGCGTCTAGGTCATCTTCTCTA", 
      "CCGAAAAGCCCAAAGGTTTTATGCAAGTTCCACGAGAATCGTCGAGCGCTGAGGTACGGCCACCGAATCACAATTCCATTCACATTCGCGGAGT", 
      "TTTCTTTGGCAGATGTCCTCCGGACTTCCGGAATGAAGCTTCAGGGCGCGCCTATCGCGAAACCCCGACGTGACGGAACGTCCTGTCCTGGACT", 
      "AGTTGCCCAAGGTTTGGACCTGCCCGGTTTTCTGGCCGCGTATGAAGTATTCCAGTAGATCCCTGCGATACTTCTCGCTACGCTAGGTCCGATA", 
      "CTCGTAAATAAACGGGATAACAAAGTGTTATCGATTATGTATAGCGGTCAGGGTGTCTGCACTTATCGGAGTCGGGAGGAGACCCTTCCTAAGA", 
      "CGCTTATGTCGTAAGTATAGCAACTAAGGGGCCATATGAACAAGTGTAGTACACGCGGCACCAGGTACCCTTAGGCGTTGTTCTTCGGAGGTTT", 
      "CTAAAGACCCCGCGTGGAATTAGAGGGAGCTCCAACGATAGAAATGAAGGAGGGGGCTTGTATAAATCCAACTGAAGGTGAACGAAGCTATAGA", 
      "CGTTTCCCTACACGGCCAATGAACAGCGGGCGCATAACATCTCTTAGCCCTGACTTCTGCGACGCATATATTTACCGGGACATGATCACAGTGT", 
      "TTACAGCACAGCAAGAAATTCCGTGCTTGGGGGATCTACCGTGCATACTCCTGAGCCCAGATATTGGCCCGTAAGAATCATGCATTGAAAATAA", 
      "TGGTCTCACCAGTGGGCACCGCGTATTGTCAGGCAAGAGCAATTAGTGATTCGTACATAGAGGTCTTCAATCTCATGTATAGAGTCGCTGGTTC", 
      "ATGCACCTAGAAGAACATAGTTAAGAAGCGTCATGGAGCTACAGAGTTCCGTGGCAATGTACAGTAGCTGTTGCCAATTTGTTGGGACTCAGCC", 
      "AGGTCTCTTCGAAATAGGTAGTGCCGGAACCATACTGAGAGTCTGCAAACGTTGCTTCGGAATGTACCTGTTGAGGACACTACTCTGACGGGCG", 
      "GCGATCACCCGCGACTCGTTCACTGGTACTTCCAATGCACTGCGCCCAATGCGCAATACGCGAATTCTCTTAAATGATCTCAAATATATGAGCT"]

d("AATTGGG", test)

test = ["AGGCGGCACATCATTATCGATAACGATTCGCCGCATTGCC",
"ATCCGTCATCGAATAACTGACACCTGCTCTGGCACCGCTC",
"AAGCGTCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT",
"AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG",
"AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT",
"AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT",
"AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG",
"AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"]
GreedyMotifSearch(test, 5,8)

test = ["GCACATCATTAAACGATTCGCCGCATTGCCTCGATAGGCG",
"TCATAACTGACACCTGCTCTGGCACCGCTCATCCGTCGAA",
"AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTCGGC",
"AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG",
"AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT",
"AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT",
"AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTGAGGCG",
"AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCG"]

GreedyMotifSearch(test, 5,8)

import random

def RandomMotifs(dna, k, t):
    sub = []
    for i in range(t-1):
        x = random.randint(1, (len(dna[i])-k))
        sub.append(dna[i][x:x+k])
    return sub


def RandomizedMotifSearch(dna,k,t):
    m = RandomMotifs(dna, k, t)
    best_motifs = m
    while True:
        profile = ProfileWithPseudocounts(m)
        m = Motifs(profile, dna)
        if Score(m) < Score(best_motifs):
            best_motifs = m
        else:
            return best_motifs 

DNA = ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"]
test = ["TGA", "GTT", "GAA", "TGT"]
profile = ProfileWithPseudocounts(test)
Motifs(profile, DNA)

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

len(test)

len(prac)

N=1000
score1=5000
for i in range(0,N):
    if i%10==0:
        print(i)
    tempmotif=RandomizedMotifSearch(text,15, 20)
    score=Score(tempmotif)
    if score<score1:
        score1=score
        best_motifs=tempmotif

print(*best_motifs, "\n")


def Motifs(profile, dna):
    kmer=len(profile['A'])
    pattern=dna.copy()
    for j in range(len(dna)):
        score=5000
        for i in range(len(dna[j])-kmer+1):
            score2=Pr(dna[j][i:i+kmer], profile)
            if score2<score:
                score=score2
                pattern[j]=dna[j][i:i+kmer]
            if score2==score:
                pattern[j]=pattern[j]
    return pattern

def HammingDistance(p, q):
    dist = 0
    if len(p)==len(q):
        for i in range(len(p)):
            print(p, q)
            if p[i]!=q[i]:
                dist+=1
            else:
                dist+=0
    return dist

p1 = (600-15)/(600-15+1) #probability of not capturing the 15-mer in one string
p2 = p1**10 #10 cases of not capturing the 15-mer
p3 = 1 - p2 #subtract to find the probability of capturing it once
p3

p4 = 1/586 * ((585/586) ** 9)
p3-p4*10

test = [4, 2, 1]

random.uniform(0,1)

def Random(prob_prof):
    new_prof = {}
    total = 0
    for key in prob_prof:
        total += prob_prof[key]
    for key in prob_prof:
        new_prof[key] = prob_prof[key]/total
    sum_prof = {}
    y = 0 #add them gradually together so you can look for ones that are in a range
    for key in new_prof:
        y += new_prof[key]
        sum_prof[key] = y
    randomint = random.uniform(0,1)
    for key in sum_prof:
        if randomint <= sum_prof[key]:
            #print(randomint, key)
            return key

Random(test)

test.keys()

test["seq1"]

test = {"seq1": .6,
 "seq2": .3,
 "seq3": .1}

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

# first, import the random package
import random

for i in range(0,N):
    if i%10==0:
        print(i)
    tempmotif=RandomizedMotifSearch(text,15, 20)
    score=Score(tempmotif)
    if score<score1:
        score1=score
        best_motifs=tempmotif

# Input:  Integers k, t, and N, followed by a collection of strings dna
# Output: GibbsSampler(dna, k, t, N)
def GibbsSampler(dna, k, t, N):
    score1 = 5000
    bestmotif = []
    for i in range(0,70):
        tempmotif=RandomizedMotifSearch(dna,k,t)
        score=Score(tempmotif)
        if score<score1:
            score1=score
            bestmotif=tempmotif
    for j in range(0, N):
        i = random.randint(0, t-1)
        reducedmotif=[]
        for b in range(0,t):
            if b!= i:
                reducedmotif.append(bestmotif[b])
        prof = ProfileWithPseudocounts(reducedmotif)
        mi = ProfileGeneratedString(dna[i], prof, k)
        reducedmotif.insert(i, mi)
        if Score(reducedmotif) < Score(bestmotif):
            bestmotif=reducedmotif
    return bestmotif

# place all subroutines needed for GibbsSampler below this line
def RandomizedMotifSearch(dna,k,t):
    m = RandomMotifs(dna, k, t)
    best_motifs = m
    while True:
        profile = ProfileWithPseudocounts(m)
        m = Motifs(profile, dna)
        if Score(m) < Score(best_motifs):
            best_motifs = m
        else:
            return best_motifs 


def Normalize(Probabilities):
    temp = Probabilities.copy()
    val=0
    for key in Probabilities.keys():
        val+=Probabilities[key]
    for key in Probabilities.keys():
        temp[key]=temp[key]/val
    return temp


def WeightedDie(Probabilities):
    kmer=''
    p=random.uniform(0,1)
    y = 0
    newdic = Probabilities.copy()
    for key in Probabilities.keys():
        y += Probabilities[key]
        newdic[key]= y
    for key in Probabilities.keys():
        if p <= newdic[key]:
            kmer = key
            return kmer

        
def Pr(text, profile):
    p=1
    for i in range(len(text)):
        p=p*profile[text[i]][i]
    return p  

def ProfileGeneratedString(text, profile, k):
    n=len(text)
    probabilities={}
    for i in range(0, n-k+1):
        probabilities[text[i:i+k]]=Pr(text[i:i+k], profile)
    probabilities=Normalize(probabilities)
    return WeightedDie(probabilities)

def RandomMotifs(dna, k, t):
    sub = []
    for i in range(t):
        x = random.randint(1, (len(dna[i])-k))
        sub.append(dna[i][x:x+k])
    return sub


def Motifs(profile, dna):
    kmer=len(profile['A'])
    pattern=dna.copy()
    for j in range(len(dna)):
        score=-1
        for i in range(len(dna[j])-kmer+1):
            score2=Pr(dna[j][i:i+kmer], profile)
            if score2<score:
                score=score2
                pattern[j]=dna[j][i:i+kmer]
            if score2==score:
                pattern[j]=pattern[j]
    return pattern

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
            count[symbol][j] += 1#Add one every time this position occurs in the matrix at this position
    # your code here
    return count

def ProfileWithPseudocounts(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = {}
    profile = CountWithPseudocounts(motifs)
    sub= CountWithPseudocounts(motifs)
    for symbol in "ACGT":
        for y in range(k):
            profile[symbol][y]=float(sub[symbol][y])/(t+4)
    # insert your code here
    return profile

def ProfileMostProbableKmer(text, k, profile):
    score=-1
    pattern=text[0:k]
    for i in range(len(text)-k+1):
        score2=Pr(text[i:i+k], profile)
        if score2>score:
            score=score2
            pattern=text[i:i+k]
        if score2==score:
            pattern=pattern
    return pattern

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

def Pr(text, profile):
    p=1
    for i in range(len(text)):
        p=p*profile[text[i]][i]
    return p

def Consensus(motifs):
    k = len(motifs[0])
    count = CountWithPseudocounts(motifs)
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

sample = open(folder+"/upstream250.txt", 'r').readlines()

" for python3"
from itertools import groupby

def fasta_iter(fasta_name):
"""
modified from Brent Pedersen
Correct Way To Parse A Fasta File In Python
given a fasta file. yield tuples of header, sequence
"""
"first open the file outside "
    fh = open(fasta_name)

# ditch the boolean (x[0]) and just keep the header or sequence since
# we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
    # drop the ">"
        headerStr = header.__next__()[1:].strip()

    # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

    fiter = fasta_iter('testset.fas')

    for ff in fiter:
        headerStr, seq = ff
        print(headerStr)