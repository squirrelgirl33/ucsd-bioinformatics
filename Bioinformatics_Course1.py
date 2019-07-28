# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 22:13:08 2019

@author: Kathryn
"""

sample = open('dataset_2_7.txt', 'r').read().split('\n')
sample

Text = sample[0]
Pattern = sample[1]

len(Pattern)

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Pattern == Text[i:(i+len(Pattern))]:
            count += 1
    return count

text = "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT"

PatternCount(text, "TTT")

PatternCount(Text, Pattern)

def FrequentWords(Text, k):
    FrequentPatterns = {}
    count = 0
    for i in range(len(Text)-k+1):
        pattern = Text[i:i+k]
        counti = PatternCount(Text, pattern)
        FrequentPatterns[pattern] = counti
    pat = []
    m = max(FrequentPatterns.values())
    for key in FrequentPatterns:
        if FrequentPatterns[key]==m:
            pat += [key]
    return pat

Text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"

Text[0:4]

FrequentWords(Text, 4)

bob = {"george":4,
       "kevin":3,
       "craig":8}

bob.keys()
bob["george"] += 1

bob["larry"] = 1

sample = open('dataset_2_10.txt', 'r').read().split('\n')
Text = sample[0]
k = sample[1]

FrequentWords(Text, 12)

Text[0]

Text[::-1]

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

sample = open('dataset_3_2.txt', 'r').read().split('\n')
Text = sample[0]

ReverseComplement(Text)

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

sample = open('dataset_3_5.txt', 'r').read().split('\n')
pattern = sample[0]
text = sample[1]

#Find k-mers that are repeated at least t times in a L length clump

def Clump(genome, k, L, t):
    dic = {}
    sequences = {}
    for i in range(L-k+1):
        if genome[i:(i+k)] in dic:
            dic[genome[i:(i+k)]] += 1
#            if dic[genome[i:(i+k)]] >= t:
 #               sequences.append(genome[i:(i+k)])
        else:
            dic[genome[i:(i+k)]] = 1
        #print(dic)
    for i in range(L, len(genome)):
        for key in dic:
            if dic[key] >= t:
                sequences[key] = ""
                #print(sequences)
                #print(dic)
        dic[genome[i-L:(i-L+k)]] -= 1
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

sample = open('dataset_4_5.txt', 'r').read().split('\n')
text = sample[0]
sample[1]
ecoli = open('E_coli.txt', 'r').read().split('\n')
ecoli1=ecoli[0]
ecoli1[0:30]
Clump(ecoli1, 9, 500, 3)

PatternCount("ACTGTACGATGATGTGTGTCAAAG", "TGT")
FrequentWords("CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA",3)
PatternMatching("CGC","ATGACTTCGCTGTTACGCGC")

def PatternToNumber(Pattern):
    Pattern = Pattern[::-1]
    q = len(Pattern) - 1
    Number = 0
    while q >= 0 :
        if Pattern[q] == 'A':
            num = 0
        if Pattern[q] == 'C':
            num = 1
        if Pattern[q] == 'G':
            num = 2
        if Pattern[q] == 'T':
            num = 3
        Number += (num * (4**q))
        q -= 1
    return(Number)

Pattern[::-1]

912/(4)
228/4
57%4
14.25%4
3.5625%4
0.890625%4

Pattern = "ATGCAA"
print(PatternToNumber(Pattern))

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


def ComputingFrequencies(Text, k):
    freq = {}
    pattern = ""
    list = []
    for i in range(4**k):
        freq[i] = 0
    #   print(freq)
    for i in range(len(Text)-k+1):
        pattern = Text[i:i+k]
        j = PatternToNumber(pattern)
        freq[j]+=1
    for i in range(4**k):
        list.append(freq[i])
    return freq

x = ComputingFrequencies("ACGCGGCTCTGAAA", 2)

x = [3, 4, 1, 3]

max(x)

max(x.values())

sample = open("dataset_2994_5.txt", 'r').read().split("\n")
Text = sample[0]
k = sample[1]
k
x = ComputingFrequencies(Text, 5)
f = open("answer.txt", 'w')
f.write()

def FasterFrequentWords(Text, k):
    FrequentPatterns = []
    FrequencyArray = ComputingFrequencies(Text, k)
    maxCount = max(FrequencyArray.values())
    #print(maxCount)
    for i in range(4**k):
        if FrequencyArray[i] == maxCount:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns

FrequentWords(Text, 4)
FasterFrequentWords(Text, 4)

def LastSymbol(Pattern):
    pattern2=Pattern[len(Pattern)-1]
    return pattern2

def Prefix(Pattern):
    pattern2 = Pattern[0:(len(Pattern)-1)]
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


def PatternToNumber2(Pattern):
    if Pattern is not "ACGT":
        return 0
    symbol = LastSymbol(Pattern)
    prefix = Prefix(Pattern)
    return 4 * PatternToNumber2(Prefix) + SymbolToNumber(symbol)

PatternToNumber(pattern)

sample = open("dataset_3010_2.txt", 'r').read().split("\n")
pattern = sample[0]

NumberToPattern(8308, 11)

sample = open("dataset_3010_5.txt", 'r').read().split("\n")

z = {1: 3, 20: 2, 4: 2}

sorted(z)

max(z)

z[1]+=1
z[1]

z[2]+=1

max(z.values())

def FindingFrequentWordsBySorting(Text, k):
    FrequentPatterns = []
    Count = {}
    for i in range(0, len(Text)-k):
        Pattern = Text[i:i+k]
        Index = PatternToNumber(Pattern)
        if Index in list(Count.keys()):
            Count[Index] += 1
        else:
            Count[Index]=1
    SortedIndex = sorted(Count)
#    print(Count)
    #print(SortedIndex)
#    for i in range(1, len(SortedIndex)):
#        if SortedIndex[i] == SortedIndex[i-1]:
#            Count[SortedIndex[i]] = Count[SortedIndex[i-1]] +1 
#            print(Count[SortedIndex[i]])
    maxCount = max(Count.values())
    #print(maxCount)
    for key in Count:
        if Count[key]==maxCount:
            Pattern = NumberToPattern(key, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
    
FindingFrequentWordsBySorting(Text, 4)    

def BetterClumpFinding(Genome, k, L, t):
    FrequentPatterns = []
    Clump = {}
    for i in range(0, 4**k-1):
        Clump[i]=0
    Text = Genome[0:L]
    FrequencyArray=ComputingFrequencies(Text,k)
    for key in Clump:
        if FrequencyArray[key] >= t:
            Clump[key] += 1
    for i in range(1, len(Genome)-L):
        FirstPattern = Genome[i-1:i-1+k]
        index = PatternToNumber(FirstPattern)
        FrequencyArray[index] -= 1
        LastPattern = Genome[i+L-k:i+L]
        index = PatternToNumber(LastPattern)
        FrequencyArray[index] += 1
        if FrequencyArray[index] >= t:
            Clump[index] = 1
    for key in Clump:
        if Clump[key] ==1:
            Pattern = NumberToPattern(key, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
            

    
def Skew(Genome):
    skew = 0
    G = 0
    C = 0
    result = [0]
    for i in range(len(Genome)):
        if Genome[i] == "G":
            G +=1
        if Genome[i] == "C":
            C +=1
        skew = G - C
        result.append(skew)
    return result

result = Skew("GAGCCACCGCGATA")
print(*result, sep=" ")


def MinSkew(Genome):
    skew = Skew(Genome)
    low = min(skew)
    pos = []    
    for i in range(len(skew)):
        if skew[i]==low:
            pos.append(i)
    return pos

MinSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")


sample = open("dataset_7_6.txt", 'r').read().split("\n")
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

sample=open("dataset_9_3.txt", 'r').read().split("\n")
p = sample[0]
q=sample[1]

HammingDistance(p,q)

def ApproxPatternMatch(Pattern, Text, d):
    result = []
    n = len(Pattern)
    for i in range(len(Text)-n+1):
        x = HammingDistance(Pattern, Text[i:i+n])
        #print(x)
        if x <= d:
            result.append(i)
    x = len(result)
    return x

def ApproximatePatternMatching(Text, Pattern, d):
    result = []
    n = len(Pattern)
    for i in range(len(Text)-n+1):
        x = HammingDistance(Pattern, Text[i:i+n])
        #print(x)
        if x <= d:
            result.append(i)
    x = len(result)
    return result

ApproxPatternMatch("ATTCTGGA","CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3)

ApproxPatternMatch("CCC", "CATGCCATTCGCATTGTCCCAGTGA", 2)

Pattern = "ATTCTGGA"
Text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"

sample=open("dataset_9_4.txt", 'r').read().split("\n")
Pattern = sample[0]
Text=sample[1]
sample[2]

x = ApproxPatternMatch(Pattern, Text, 6)
print(*x, sep = " ")

def Count2(Text, Pattern):
    x = len(ApproxPatternMatch(Pattern, Text, 2))
    return x

Count2("CATGCCATTCGCATTGTCCCAGTGA", "CCC")

ApproxPatternMatch("GAGG","TTTAGAGCCTTCAGAGG", 2)

sample = open("dataset_9_6.txt", 'r').read().split("\n")
Pattern= sample[0]
Text=sample[1]
n=sample[2]

ApproxPatternMatch(Pattern, Text, 2)

n = {"a":1, "b":2, "c":3}
list(n.keys())

import itertools

bases=['A','T','G','C']
x = [''.join(p) for p in itertools.product(bases, repeat=4)]
len(x[0])

def FrequentMismatch(Text, k, d):
    bases=['A','T','G','C']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    dic = {}        
    for j in range(len(kmers)):
        for i in range(len(Text)-k+1):
            x = HammingDistance(kmers[j], Text[i:i+k])

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

Text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
   4 1
   
FrequentMismatch(Text, 4, 1)

sample = open("dataset_9_7.txt", 'r').read().split('\n')
Text=sample[0]
sample[1]
sample[2]


x = FrequentMismatch(Text, 6, 3)
print(*x, sep=" ")

def FrequentMismatchRC(Text, k, d):
    bases=['A','T','G','C']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    dic = {}        
    RCText = ReverseComplement(Text)
    for j in range(len(kmers)):
        for i in range(len(Text)-k+1):
            x = HammingDistance(kmers[j], Text[i:i+k])
            y = HammingDistance(kmers[j], Text[i:i+k])
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

sample=open("dataset_9_8.txt", 'r').read().split("\n")
text=sample[0]
sample[1]
FrequentMismatchRC(text, 5, 2)
            
#    Neighbors(Pattern, d)
#        if d = 0
#            return {Pattern}
#        if |Pattern| = 1 
#            return {A, C, G, T}
#        Neighborhood ← an empty set
#        SuffixNeighbors ← Neighbors(Suffix(Pattern), d)
#        for each string Text from SuffixNeighbors
#            if HammingDistance(Suffix(Pattern), Text) < d
#                for each nucleotide x
#                    add x • Text to Neighborhood
#            else
#                add FirstSymbol(Pattern) • Text to Neighborhood
#        return Neighborhood
def Neighbors(Pattern, d):
    if d == 0:
        return Pattern
    elif len(Pattern) == 1: 
        return {'A', 'C', 'G', 'T'}
    Neighborhood = []
    SuffixNeighbors = Neighbors(Pattern[1:], d)
    for Text in SuffixNeighbors:
        nucleotides = {'A', 'C', 'G', 'T'}
        if HammingDistance(Pattern[1:], Text) < d:
            for x in nucleotides:
                Neighborhood.append(x + Text)
        else:
            Neighborhood.append(Pattern[0] + Text)
    return Neighborhood

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

def MotifEnumeration(Dna, k, d):
    motifs = []
    for i in range(len(Dna)):
        if i == 0:
            for j in range(len(Dna[i])-k+1):
                Pattern = (Dna[i][j:j+k])
                Patternp=Neighbors(Pattern, d)
                if isinstance(Patternp, str)==True:  #If there's only one value, it's a string not a list.
                    motifs.append(Patternp)
                else:
                    motifs.extend(Patternp)
    motifs = list(set(motifs))
    temp = {}
    removemotif = []
    for j in range(1, len(Dna)):
        for motif in motifs:
            #print(motif)
            for i in range(len(Dna[j])-k+1):
                if HammingDistance(motif, Dna[j][i:i+k]) <= d:
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

text = open("dataset_156_8.txt", 'r').read().split("\n")
thing = text[1:len(text)-1]
text[0]

MotifEnumeration(thing, 5, 1)

{'TTG': 1, 'ACT': 1, 'ATA': 1, 'GGT': 0, 'TTC': 2, 'AGT': 0, 'TGA': 2, 'ATC': 0, 'GTT': 1, 'TTT': 2, 'TGC': 1, 'ATG': 0, 'TCT': 1, 'CGG': 0, 'AGC': 1, 'AAT': 0, 'CGC': 1, 'TTA': 1, 'GGG': 0, 'GCC': 1, 'CTT': 2, 'TGG': 1, 'TAT': 0, 'AGG': 0, 'ATT': 1}

10/4
1 2 3 4 5 6 7 8 9 10 
A A A G C T G C T G  
7*15

import math

def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0]) #get the length of each motif row (x axis)
    for symbol in "ACGT": #Make A C G T a key
        count[symbol] = [] #make an empty list for each key
        for j in range(k): #for the length of the motif...
            count[symbol].append(0) #...for the key, put a zero for each position in the list for now
    t = len(Motifs) #get the number of keys in the dictionary (y axis)
    for i in range(t): #for every key...
        for j in range(k): #for list position in each key...
            symbol = Motifs[i][j] #for this specific position for a key...
            count[symbol][j] += 1 #Add one every time this position occurs in the matrix at this position
    # your code here
    return count


def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = Count(Motifs)
    sub= Count(Motifs)
    for x in range(t): #for key in position x...
        for y in range(k): #for list position in that key...
            symbol=Motifs[x][y]
            profile[symbol][y]=(sub[symbol][y])/t
    # insert your code here
    return profile

Motifs = [
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

test = Profile(Motifs)

test["A"][1]

math.log2(test["A"][1])

len(test["A"])

def entropy(Motifs):
    prof = Profile(Motifs)
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


entropy(Motifs)

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

def Motifs(Pattern, Dna):
    motifs = []
    print(Pattern)
    k= len(Pattern)
    if isinstance(Dna, str)==True: #verify if the input is a list or a string
        DNA = []
        DNA.append(Dna)
        Dna = DNA
    for seq in Dna:
        score = k
        temp_motif = {}
        for i in range(len(seq)-k+1):
            new_score = HammingDistance(Pattern, seq[i:i+k])
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

Dna = ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG"]    

Motifs("AAA", Dna)

Motifs("GATTCTCA", "GCAAAGACGCTGACCAA")

Motifs("AAG", "GCAATCCTCAGC")

def d(Pattern, Dna):
    motifs = Motifs(Pattern, Dna)
    if len(motifs)==0:
        return len(Pattern)
    score = 0
    for motif in motifs:
        x = HammingDistance(Pattern, motif)
        score += x
    return score

d("AAA", Dna)

d("GAC", Dna)

def MedianString(Dna, k):
    kmers = []
    Median = ""
    if isinstance(Dna, str)==True:
        distance = len(Dna)
        for i in range(len(Dna)-k+1):
            kmers.append(Dna[i:i+k])
    else:
        distance = len(Dna[0])
        for seq in Dna:
            for i in range(len(seq)-k+1):
                kmers.append(seq[i:i+k])
    kmers = list(dict.fromkeys(kmers)) #remove duplicates
    #print(kmers)
    for kmer in kmers:
        if distance > d(kmer, Dna):
            distance = d(kmer, Dna)
            Median = kmer
            print(kmer)
            print(distance)
    return Median

MedianString(text, 7)
MedianString2(text, 7)

#Search kmers not in the main text
def MedianString2(Dna, k):
    kmers = []
    Median = []
    if isinstance(Dna, str)==True:
        distance = len(Dna)*10
        for i in range(len(Dna)-k+1):
            kmers.append(Dna[i:i+k])
    else:
        distance = len(Dna[0])*len(Dna)
        for seq in Dna:
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
        if distance > d(kmer, Dna):
            distance = d(kmer, Dna)
            Median = []
            Median.append(kmer)
            #print(kmer)
            #print(distance)
        elif distance == d(kmer, Dna):
            Median.append(kmer)
    return Median

test = ["ATA","ACA","AGA","AAT","AAC"]

d("AAA", test)

MedianString2(test, 3)
Motifs("GTC", test)
Motifs("AAA", test)
d("GTC", test)
d("GAT", test)

Motifs("GAT", test)

d("AAA", test)

MedianString2(Dna, 3)

3
Dna = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"]

Motifs("GAC", Dna)

d("GAC", Dna)

MedianString(Dna, 3)

sample = open("dataset_158_9.txt", 'r').read().split("\n")
text = sample[0]
text2 = sample[1]
Dna = ['TGTTGGAAATATCACGTGGACCCTCTCGTGCAGTGGATATAA',
 'TGTCCGCGGCAACACCCTAGTGCCGGAGCCGGATCTGTGCAA',
 'TGGGCGTGGACGTACCCTGACCCCCAATTTCGTGAGTATCAA',
 'CTTCATAATTATAACCCTATTTTCAGTCACGAGCACCCTAAG',
 'TACCCTCCAATACAGGACGTACCTAGGTGCGCGCGCCAGACG',
 'TGTATTAACCAAAGATAGTCCCGGTACCCTGCGCGCTCGGCA',
 'AACCCTCCACTCCGGAATAATGTTTACATACCCCCTTTAATT',
 'CTGCTAGCATAGAGCATTGCTTGTTTGGCTGCCCTACACCCT',
 'GTGCAGGGCGTAAGATAGAACCCTGCTACGAGCCGTGGAAAC',
 'CAGAACACCCGCGCAAGAGACCCTTCGGTAGGAATCATGCCT']

MedianString(Dna, 6)

def Pr(Text, Profile):
    p=1
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
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

sample = open("dataset_159_3.txt", 'r').read().split("\n")
text=sample[0]

profile = {
         "A": [0.303, 0.237, 0.25, 0.197, 0.171, 0.25, 0.342, 0.237, 0.092, 0.224, 0.289, 0.303, 0.289],
 'C' : [0.184, 0.316, 0.289, 0.276, 0.329, 0.237, 0.079, 0.303, 0.263, 0.276, 0.25, 0.237, 0.303],
 'G' : [0.289, 0.237, 0.224, 0.25, 0.237, 0.25, 0.355, 0.211, 0.25, 0.276, 0.237, 0.211, 0.263],
 'T' : [0.224, 0.211, 0.237, 0.276, 0.263, 0.263, 0.224, 0.25, 0.395, 0.224, 0.224, 0.25, 0.145]}

ProfileMostProbableKmer(text, 13, profile)

def CountWithPseudocounts(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0]) #get the length of each motif row (x axis)
    for symbol in "ACGT": #Make A C G T a key
        count[symbol] = [] #make an empty list for each key
        for j in range(k): #for the length of the motif...
            count[symbol].append(1) #...for the key, put a one for each position in the list for now
    t = len(Motifs) #get the number of keys in the dictionary (y axis)
    for i in range(t): #for every key...
        for j in range(k): #for list position in each key...
            symbol = Motifs[i][j] #for this specific position for a key...
            count[symbol][j] += 1 #Add one every time this position occurs in the matrix at this position
    # your code here
    return count

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = CountWithPseudocounts(Motifs)
    sub= CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        for y in range(k):
            profile[symbol][y]=float(sub[symbol][y])/(t+4)
    # insert your code here
    return profile

def GreedyMotifSearch(Dna, k, t):
    BestMotifs=[]
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k]) #sets BestMotifs to the first three k-mers in each row
    n=len(Dna[0]) #how long each row is
    for i in range(n-k+1): #iterate through each row as a k-mer
        Motifs=[] #set up an empty motif matrix
        Motifs.append(Dna[0][i:i+k]) #fill Motif matrix with k-mer
        for j in range(1,t): #go through each row in the matrix
            P=ProfileWithPseudocounts(Motifs[0:j]) #give a score for each nucleotide for the motif
            Motifs.append(ProfileMostProbableKmer(Dna[j],k,P)) #For the j row, get the motif with the highest score using the profile scores
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs=Motifs
    return BestMotifs

sample=open("dataset_159_5.txt", 'r').read().split("\n")
Dna = ['GTACAAGTGAACGGCTGTTTACAGTCGGAGTAGCAGGTCTTCTGTTATTGGCGATGGCCGGTATATACGAAAGGCGAAGCGTGGTGCGATCTCTTTGCCTGTGCTCTCCGAATCCCTAGGGAGCGCAATCATAGGCTCTCCCTGCAAATGGAGCGT',
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

text = GreedyMotifSearch(Dna, 12, 25)

print(*text, sep=" ")

print(*list(sequences.keys()), sep=" ")

sample = open("dataset_160_9.txt", 'r').read().split("\n")
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

sample = open("dataset_5164_1.txt", 'r').read().split("\n")
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

Motifs("AAA", test)

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

def RandomMotifs(Dna, k, t):
    sub = []
    for i in range(t-1):
        x = random.randint(1, (len(Dna[i])-k))
        sub.append(Dna[i][x:x+k])
    return sub


def RandomizedMotifSearch(Dna,k,t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 

DNA = ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"]
test = ["TGA", "GTT", "GAA", "TGT"]
profile = ProfileWithPseudocounts(test)
Motifs(profile, DNA)

sample = open("dataset_161_5.txt", 'r').read().split("\n")
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
        BestMotifs=tempmotif

print(*BestMotifs, "\n")


def Motifs(Profile, Dna):
    kmer=len(Profile['A'])
    pattern=Dna.copy()
    for j in range(len(Dna)):
        score=5000
        for i in range(len(Dna[j])-kmer+1):
            score2=Pr(Dna[j][i:i+kmer], Profile)
            if score2<score:
                score=score2
                pattern[j]=Dna[j][i:i+kmer]
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

def Random(ProbProf):
    NewProf = {}
    total = 0
    for key in ProbProf:
        total += ProbProf[key]
    for key in ProbProf:
        NewProf[key] = ProbProf[key]/total
    SumProf = {}
    y = 0 #add them gradually together so you can look for ones that are in a range
    for key in NewProf:
        y += NewProf[key]
        SumProf[key] = y
    randomint = random.uniform(0,1)
    for key in SumProf:
        if randomint <= SumProf[key]:
            #print(randomint, key)
            return key

Random(test)

test.keys()

test["seq1"]

test = {"seq1": .6,
 "seq2": .3,
 "seq3": .1}

#Gibbs Sampler Problem

sample = open("dataset_163_4.txt", 'r').read().split("\n")
sample = open("gibbs.txt", 'r').read().split("\n")
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
        BestMotifs=tempmotif

# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    score1 = 5000
    bestmotif = []
    for i in range(0,70):
        tempmotif=RandomizedMotifSearch(Dna,k,t)
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
        Mi = ProfileGeneratedString(Dna[i], prof, k)
        reducedmotif.insert(i, Mi)
        if Score(reducedmotif) < Score(bestmotif):
            bestmotif=reducedmotif
    return bestmotif

# place all subroutines needed for GibbsSampler below this line
def RandomizedMotifSearch(Dna,k,t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 


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

        
def Pr(Text, Profile):
    p=1
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p  

def ProfileGeneratedString(Text, profile, k):
    n=len(Text)
    probabilities={}
    for i in range(0, n-k+1):
        probabilities[Text[i:i+k]]=Pr(Text[i:i+k], profile)
    probabilities=Normalize(probabilities)
    return WeightedDie(probabilities)

def RandomMotifs(Dna, k, t):
    sub = []
    for i in range(t):
        x = random.randint(1, (len(Dna[i])-k))
        sub.append(Dna[i][x:x+k])
    return sub


def Motifs(Profile, Dna):
    kmer=len(Profile['A'])
    pattern=Dna.copy()
    for j in range(len(Dna)):
        score=-1
        for i in range(len(Dna[j])-kmer+1):
            score2=Pr(Dna[j][i:i+kmer], Profile)
            if score2<score:
                score=score2
                pattern[j]=Dna[j][i:i+kmer]
            if score2==score:
                pattern[j]=pattern[j]
    return pattern

def CountWithPseudocounts(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0]) #get the length of each motif row (x axis)
    for symbol in "ACGT": #Make A C G T a key
        count[symbol] = [] #make an empty list for each key
        for j in range(k): #for the length of the motif...
            count[symbol].append(1) #...for the key, put a one for each position in the list for now
    t = len(Motifs) #get the number of keys in the dictionary (y axis)
    for i in range(t): #for every key...
        for j in range(k): #for list position in each key...
            symbol = Motifs[i][j] #for this specific position for a key...
            count[symbol][j] += 1#Add one every time this position occurs in the matrix at this position
    # your code here
    return count

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = CountWithPseudocounts(Motifs)
    sub= CountWithPseudocounts(Motifs)
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

def Score(Motifs):
    con = Consensus(Motifs)
    k = len(Motifs[0])
    j=len(Motifs)
    value=0
    for x in range(k):
        for y in range(j):
            if Motifs[y][x]!=con[x]:
                value+=1
    return value

def Pr(Text, Profile):
    p=1
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p

def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

sample = open("upstream250.txt", 'r').readlines()

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