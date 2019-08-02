# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 21:12:52 2019

@author: Kathryn
"""

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

    
def FrequencyMap(text, k):
    freq = {}
    n = len(text)
    for i in range(n-k+1):
        pattern = text[i:i+k]
        if pattern in freq:
            freq[pattern]+=1
        else:
            freq[pattern]=1
    return freq

Text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
k = 10

test = FrequencyMap(Text, 3)

test["ATA"]

FrequencyMap(Text, k)

m = max(freq.values())

max(test.values())


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key]==m:
            words += [key]
        # add each key to words whose corresponding frequency value is equal to m
    return words

FrequentWords(Text,k)

def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern

len("GCATG")

def Reverse(Pattern):
    RevPat = ""
    n = len(Pattern)
    for i in range(0,n):
        RevPat = Pattern[i] + RevPat
    return RevPat

Reverse("GCATG")

def Complement(Pattern):
    CompPat = ""
    n = len(Pattern)
    for i in range(0, n):
        if Pattern[i]=="A":
            CompPat = CompPat + "T"
        if Pattern[i]=="T":
            CompPat = CompPat + "A"
        if Pattern[i]=="C":
            CompPat= CompPat + "G"
        if Pattern[i]=="G":
            CompPat = CompPat + "C"
    return CompPat

Complement("GCATG")

def ReverseComplement(Pattern):
    NewSeq = ""
    NewSeq2 = ""
    for char in Pattern:
        NewSeq = char + NewSeq
    for char in NewSeq:
        if char=="A":
            NewSeq2 = NewSeq2 + "T"
        if char=="T":
            NewSeq2 = NewSeq2 + "A"
        if char=="C":
            NewSeq2 = NewSeq2 + "G"
        if char=="G":
            NewSeq2 = NewSeq2 + "C"
    return NewSeq2

def RevComp(Pattern):
    NewSeq = Reverse(Pattern)
    NewSeq = Complement(NewSeq)
    return NewSeq

RevComp("GCATG")

ReverseComplement("GCATG")

test = "GCATG"


range(len("GCATG"))

pat = "ATAT"
gen = "GATATATGCATATACTT"

def PatternMatching(Pattern, Genome):
    sub = []
    ng = len(Genome)
    np = len(Pattern)
    for i in range(ng - np+1):
        if Genome[i:i+np]==Pattern:
            sub.append(i)
    return sub

PatternMatching(pat, gen)

test = []
test.append(2)
test

Text2 = "aactctatacctcctttttgtcgaatttgtgtgatttatagagaaaatcttattaactgaaactaaaatggtaggtttggtggtaggttttgtgtacattttgtagtatctgatttttaattacataccgtatattgtattaaattgacgaacaattgcatggaattgaatatatgcaaaacaaacctaccaccaaactctgtattgaccattttaggacaacttcagggtggtaggtttctgaagctctcatcaatagactattttagtctttacaaacaatattaccgttcagattcaagattctacaacgctgttttaatgggcgttgcagaaaacttaccacctaaaatccagtatccaagccgatttcagagaaacctaccacttacctaccacttacctaccacccgggtggtaagttgcagacattattaaaaacctcatcagaagcttgttcaaaaatttcaatactcgaaacctaccacctgcgtcccctattatttactactactaataatagcagtataattgatctga"

PatternCount(Text2, "ATGATCAAG")
PatternCount(Text2, "CTTGATCAT")

FrequentWords(Text2, 9)

PatternCount("GACCATCAAAACTGATAAACTACTTAAAAATCAGT","AAA")
FrequentWords("TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT", 3)
RevComp("GATTACA")

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

e_coli = open("E_coli.txt", "r").read()

len(e_coli)

SymbolArray(e_coli, "C")

def SymbolArray(Genome, symbol):
    array = []
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        count = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
        array.append(count)
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def SkewArray(Genome):
    Skew=[0]*(len(Genome)+1)
    for i in range(len(Genome)):
        if Genome[i] == "A" or "T":
            Skew[i+1]=Skew[i]
        if Genome[i] == "G":
            Skew[i+1]=Skew[i]+1
        if Genome[i]=="C":
            Skew[i+1]=Skew[i]-1
    return Skew

test5="C"
test4="CATGGGCATCGGCCATACGCC"

SkewArray(test4)

prac=[1,5,6,20,-400, 2]
min(prac)

len(prac[0:])

def MinimumSkew(Genome):
    pos=[]
    x=SkewArray(Genome)
    low=min(x)
    for i in range(len(x)):
        if x[i]==low:
            pos.append(len(x[0:(i)]))
    return pos

test6="TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
test7="AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTC"
min(SkewArray(test7))


MinimumSkew(test6)

SkewArray(test6)

test1 = "GGGCCGTTGGT"
test2 = "GGACCGTTGAC"

test1[1]

def HammingDistance(p, q):
    score=0
    for i in range(len(p)):
        if p[i]==q[i]:
            score+=0
        else:
            score+=1
    return score

HammingDistance(test1, test2)

test1="ATTCTGGA"
test2="CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text) - len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions
    
ApproximatePatternMatching(test2, test1, 3)

def ApproximatePatternCount(Pattern, Text, d):
    count=0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count+=1
    return count
    
ApproximatePatternCount(test1, test2, 4)

test1 = "CATTCCAGTACTTCGATGATGGCGTGAAGA"
MinimumSkew(test1)

x=0
for y in range(0,5):
    x+=y
print(x)

test1 = "TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC" 

test2 = "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"

HammingDistance(test1, test2)

a=list(range(5))
b=a
a[2]=12
b

test="ACGGGGATTACC"
prof={'A':[0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
'C':[0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
'G':[0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
'T':[0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]}

motifs = ["AAAGTCTAG","AGGGCCAAG","CCAGGATAG","CAAGAATAG"]

prof[test[0]]

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

Count(motifs)

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

Profile(motifs)
# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
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

Consensus(motifs)
# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
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
    # Insert code here
Score(motifs)

def Pr(Text, Profile):
    p=1
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p

Pr("TCGGTA",quiz)

Pr(test, prof)


test="ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
k=5
prof={'A':[0.2, 0.2, 0.3, 0.2, 0.3],
'C':[0.4, 0.3, 0.1, 0.5, 0.1],
'G':[0.3, 0.3, 0.5, 0.2, 0.4],
'T':[0.1, 0.2, 0.1, 0.1, 0.2]}


def ProfileMostProbableKmer(text, k, profile):
    score=-1
    pattern=text[0:k]
    for i in range(len(text)-k+1):
        score2=Pr(text[i:i+k], profile)
        if score2<score:
            score=score2
            pattern=text[i:i+k]
        if score2==score:
            pattern=pattern
    return pattern

ProfileMostProbableKmer(test,k,prof)
ProfileMostProbableKmer(test2,k2,prof2)

test="TTACCATGGGACCGCTGACTGATTTCTGGCGTCAGCGTGATGCTGGTGTGGATGACATTCCGGTGCGCTTTGTAAGCAGAGTTTA"
k=12
prof={'A':[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.1, 0.2, 0.3, 0.4, 0.5],
'C':[0.3, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.4, 0.3, 0.2, 0.2,  0.1],
'G':[0.2, 0.1, 0.4, 0.3, 0.1, 0.1, 0.1, 0.3, 0.1, 0.1, 0.2, 0.1],
'T':[0.3, 0.4, 0.1, 0.1, 0.1, 0.1, 0.0, 0.2, 0.4, 0.4, 0.2, 0.3]}

test2="AACCGGTT"
k2=3
prof2={'A':[1.0, 1.0, 1.0],
'C':[0.0, 0.0, 0.0],
'G':[0.0, 0.0, 0.0],
'T':[0.0, 0.0, 0.0]}



test3=["GGCGTTCAGGCA",
    "AAGAATCAGTCA",
    "CAAGGAGTTCGC",
    "CACGTCAATCAC",
    "CAATAATATTCG"]

test4=[]
for i in range(0, len(test3)):
    test4.append(test3[i][0:3])
print(test4)

def GreedyMotifSearch(Dna, k, t):
    BestMotifs=[]
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k]) #sets BestMotifs to the first three k-mers in each row
    n=len(Dna[0]) #how long each row is
    for i in range(n-k+1): #iterate through each row as a k-mer
        Motifs=[] #set up an empty motif matrix
        Motifs.append(Dna[0][i:i+k]) #fill Motif matrix with k-mer
        for j in range(1,t): #go through each row in the matrix
            P=Profile(Motifs[0:j]) #give a score for each nucleotide for the motif
            Motifs.append(ProfileMostProbableKmer(Dna[j],k,P)) #For the j row, get the motif with the highest score using the profile scores
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs=Motifs
    return BestMotifs
            
GreedyMotifSearch(test3,3,5)

import math

a=[0.2,0.2,0.9,0.1,0.1,0.1,0.3]
c=[0.1,0.6,0.4,0.1,0.2,0.4,0.6]
g=[1,1,0.9,0.9,0.1]
t=[0.7,0.2,0.1,0.1,0.5,0.8,0.7,0.3,0.4]
data_list=[a,c,g,t]

H=0.0
for j in data_list:
    for i in j:
        H=H+i*(math.log(i,2))
print (-H)

quiz={'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]}


Pr("GGGGTA", quiz)

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

len(test)

test= ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]

test[1][0]

test2[]

test2 = CountWithPseudocounts(test)
ProfileWithPseudocounts(test)

Count(test)

test2

Profile(test)

len(test)

float(test2["A"][0])/(len(test))

test=["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC",
"TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT","GGCAAGTACC"]
#GreedyMotifSearch with Pseudocounts

1.1 Motif Finding Meets Oliver Cromwell
11 out of 13 steps passed
8 out of 9 points  received
Code Challenge (3 points): Write a function GreedyMotifSearchWithPseudocounts(Dna, k, t) that takes a list of strings Dna followed by integers k and t and returns the result of running GreedyMotifSearch, where each profile matrix is generated with pseudocounts. Then add this function to Motifs.py. (Hint: Ideally, you should only need an extremely small modification to your original GreedyMotifSearch function.)

Click here for this problem's test datasets.

Sample Input:

3 5
GGCGTTCAGGCA
AAGAATCAGTCA
CAAGGAGTTCGC
CACGTCAATCAC
CAATAATATTCG
Sample Output:

TTC
ATC
TTC
ATC
TTC
Code Challenge — Write a program, test using stdin → stdout
 Right.
Now you have access to the Forum of Solutions where you can discuss your solution with others.
You just solved a difficult problem, congratulations! You can help others in the comments below.
1
# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
2
# Output: GreedyMotifSearch(Dna, k, t)
3
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
4
    BestMotifs=[]
5
    for i in range(0, t):
6
        BestMotifs.append(Dna[i][0:k]) #sets BestMotifs to the first three k-mers in each row
7
    n=len(Dna[0]) #how long each row is
8
    for i in range(n-k+1): #iterate through each row as a k-mer
9
        Motifs=[] #set up an empty motif matrix
10
        Motifs.append(Dna[0][i:i+k]) #fill Motif matrix with k-mer
11
        for j in range(1,t): #go through each row in the matrix
12
            P=ProfileWithPseudocounts(Motifs[0:j]) #give a score for each nucleotide for the motif
13
            Motifs.append(ProfileMostProbableKmer(Dna[j],k,P)) #For the j row, get the motif with the highest score using the profile scores
14
        if Score(Motifs) < Score(BestMotifs):
15
            BestMotifs=Motifs
16
    return BestMotifs
17
​
18
# Copy all needed subroutines here.  These subroutines are the same used by GreedyMotifSearch(),
19
def CountWithPseudocounts(Motifs):
20
    count = {} # initializing the count dictionary
21
    k = len(Motifs[0]) #get the length of each motif row (x axis)
22
    for symbol in "ACGT": #Make A C G T a key
23
        count[symbol] = [] #make an empty list for each key
24
        for j in range(k): #for the length of the motif...
25
            count[symbol].append(1) #...for the key, put a one for each position in the list for now
26
    t = len(Motifs) #get the number of keys in the dictionary (y axis)
27
    for i in range(t): #for every key...
28
        for j in range(k): #for list position in each key...
29
            symbol = Motifs[i][j] #for this specific position for a key...
30
            count[symbol][j] += 1 #Add one every time this position occurs in the matrix at this position
31
    # your code here
32
    return count
33
​
34
def ProfileWithPseudocounts(Motifs):
35
    t = len(Motifs)
36
    k = len(Motifs[0])
37
    profile = {}
38
    profile = CountWithPseudocounts(Motifs)
39
    sub= CountWithPseudocounts(Motifs)
40
    for symbol in "ACGT":
41
        for y in range(k):
42
            profile[symbol][y]=float(sub[symbol][y])/(t+4)
43
    # insert your code here
44
    return profile
45
​
46
def ProfileMostProbableKmer(text, k, profile):
47
    score=-1
48
    pattern=text[0:k]
49
    for i in range(len(text)-k+1):
50
        score2=Pr(text[i:i+k], profile)
51
        if score2<score:
52
            score=score2
53
            pattern=text[i:i+k]
54
        if score2==score:
55
            pattern=pattern
56
    return pattern
57
​
58
def Score(Motifs):
59
    con = Consensus(Motifs)
60
    k = len(Motifs[0])
61
    j=len(Motifs)
62
    value=0
63
    for x in range(k):
64
        for y in range(j):
65
            if Motifs[y][x]!=con[x]:
66
                value+=1
67
    return value
68
​
69
def Pr(Text, Profile):
70
    p=1
71
    for i in range(len(Text)):
72
        p=p*Profile[Text[i]][i]
73
    return p
74
​
75
def Consensus(Motifs):
76
    k = len(Motifs[0])
77
    count = CountWithPseudocounts(Motifs)
78
    consensus = ""
79
    for j in range(k):
80
        m = 0
81
        frequentSymbol = ""
82
        for symbol in "ACGT":
83
            if count[symbol][j] > m:
84
                m = count[symbol][j]
85
                frequentSymbol = symbol
86
        consensus += frequentSymbol
87
    return consensus
88
​
89
# except that you should replace Count(Motifs) and Profile(Motifs) with the new functions
90
# CountWithPseudocounts(Motifs) and ProfileWithPseudocounts(Motifs).
Your submissions
Passed: 2183 Correct submissions: 22% You got: 3 points out of 3
 69  20
Step 11
 Comments
 Solutions
  

Submission

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
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

# Copy all needed subroutines here.  These subroutines are the same used by GreedyMotifSearch(),
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

Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

len(Dna[0])
len(Dna)

x=ProfileWithPseudocounts(test)
len(x["A"])


PTC = { 'A': [0.8, 0.0, 0.0, 0.2 ],'C': [ 0.0, 0.6, 0.2, 0.0], 'G': [ 0.2 ,0.2 ,0.8, 0.0], 'T': [ 0.0, 0.2, 0.0, 0.8]}   
Dna=['TTACCTTAAC','GATGTCTGTC','ACGGCGTTAG','CCCTAACGAG','CGTCAGAGGT']
Dna2="TTACCCTAAC"
Dna[0]

def Motifs(Profile, Dna):
    kmer=len(Profile['A'])
    pattern=Dna.copy()
    for j in range(len(Dna)):
        score=-1
        for i in range(len(Dna[j])-kmer+1):
            score2=Pr(Dna[j][i:i+kmer], Profile)
            if score2>score:
                score=score2
                pattern[j]=Dna[j][i:i+kmer]
            if score2==score:
                pattern[j]=pattern[j]
    return pattern

Motifs(PTC, Dna)


    

import random


def RandomMotifs(Dna, k, t):
    sub = []
    for i in range(t):
        x = random.randint(1, (len(Dna[i])-k))
        sub.append(Dna[i][x:x+k])
    return sub

RandomMotifs(Dna, 4, 5)



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
        
test=["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
"CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
"ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
"GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
"GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
"CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
"GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
"GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
"GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
"TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

t=len(Dna)
k=15
N=10

Profile(test)
Profile(test.remove(test[2]))
Profile(test[1])
Pr(test[1], ProfileWithPseudocounts(test))

test
test2 = test.remove(test[2])


import random

RandomizedMotifSearch(test,8,5)

test2=[]
score1=5000
for i in range(0,N):
    tempmotif=RandomizedMotifSearch(dna,k,t)
    score=Score(tempmotif)
    if score<score1:
        score1=score
        BestMotifs=tempmotif

len(test2)

Score(test2)

Dna[0][81:101]

CountWithPseudocounts(test2)

x = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}

len(x)

x['A']

def Normalize(Probabilities):
    temp = Probabilities.copy()
    val=0
    for key in Probabilities.keys():
        val+=Probabilities[key]
    for key in Probabilities.keys():
        temp[key]=temp[key]/val
    return temp

y = Normalize(x)

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
        
WeightedDie(y)
        
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

test = "AAACCCAAACCC"
prof = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
k = 2

ProfileGeneratedString(test,prof,k)



Normalize(x)

def GibbsSampler(Dna, k, t, N):
    tempmotif = RandomMotifs(Dna, k, t) #choose random motifs
    bestmotif=tempmotif
    for j in range(0, N): #repeat N number of times
        i = random.randint(0, t-1)  #generate a random number in the length of the DNA string
        reducedmotif=[]
        for j in range(0,t): #make a new list for every DNA string that's not the one you previously selected
            if j!= i:
                reducedmotif.append(bestmotif[j]) #removed randomly chosen motif
        prof = ProfileWithPseudocounts(reducedmotif) #profile this new set
        Mi = ProfileGeneratedString(Dna[i], prof, k) #For the removed motif, change it
        reducedmotif.insert(i, Mi)
        if Score(reducedmotif) < Score(bestmotif): #compare score of the new set to the old set
            bestmotif=reducedmotif #if new set is better, replace the old set
    return bestmotif



k= 8 
t= 5 
N=1000

Score(Dna2)

Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

GibbsSampler(Dna, 8, 5, 1000)

RandomizedMotifSearch(Dna, 8, 5)

N=1000
score1=5000
for i in range(0,N):
    if i%10==0:
        print(i)
    tempmotif=RandomizedMotifSearch(Dna, 8, 5)
    score=Score(tempmotif)
    if score<score1:
        score1=score
        BestMotifs=tempmotif

print(BestMotifs)

for i in range(10):
    print(GibbsSampler(Dna, 8, 5, 1000))

random.randint(1,t-1)

Dna2 = Dna[0]

GibbsSampler(Dna, k, t, N)
        
import random

def RandomMotifs(Dna, k, t):
    sub = []
    for i in range(t):
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