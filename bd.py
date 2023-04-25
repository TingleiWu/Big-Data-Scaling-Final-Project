#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import time
import multiprocessing as mp


# # Step1

# In[2]:


k = 15  # set the k-mer size to 15
read_kmers = []  # initialize an empty list to store the k-mers
reads = []  # initialize an empty list to store the reads

with open('reads.fq', 'r') as f:
    # loop over the lines in the file, reading four lines at a time
    for i, line in enumerate(f):
        if i % 4 == 1:  # extract the second line of every four lines
            read = line.strip()  # remove any leading/trailing whitespace
            reads.append(read)  # add the read to the list of reads
            for j in range(len(read) - k + 1):
                kmer = read[j:j+k]  # extract a k-mer from the read
                read_kmers.append(kmer)  # add the k-mer to the list of k-mers


# In[3]:


print(len(read_kmers))


# # Step2

# In[4]:


k = 15  # set the k-mer size to 15
bin_kmers = []  # initialize an empty list to store the k-mers
refbinreads = []  # initialize an empty list to store the sequences

with open('reference_chr21_20000000_20050000.fa', 'r') as f:
    # loop over the lines in the file, skipping the first line
    for line in f.readlines()[1:]:
        # extract the sequence from the line and remove any leading/trailing whitespace
        sequence = line.strip().split(',')[2]
        refbinreads.append(sequence)  # add the sequence to the list of sequences
        for j in range(len(sequence) - k + 1):
            kmer = sequence[j:j+k]  # extract a k-mer from the sequence
            bin_kmers.append(kmer)  # add the k-mer to the list of k-mers


# In[5]:


print(len(bin_kmers))


# # Step3

# In[6]:


all_kmers = read_kmers + bin_kmers
kmer_set = set(all_kmers)
distinct_kmers = len(kmer_set)

#print(len(all_kmers))
print("Number of distinct kmers: ", distinct_kmers)


# # Step4

# In[7]:


encode_result = []
for kmer in kmer_set:
    current_result = []
    current_result.append(kmer)
    for read in reads:
        if kmer in read:
            current_result.append(1)
        else:
            current_result.append(0)
    encode_result.append(current_result)   


# In[8]:


columnnames = []
for i in range(0,2000):
    columnnames.append('Read'+str(i))


# In[9]:


columnnames = ["kmer"] + columnnames


# In[10]:


df = pd.DataFrame(encode_result, columns=columnnames)


# In[11]:


pd.set_option('display.max_columns', None)
print(df.head(5))


# # Step5
# 

# In[12]:


refbinencode_result = []
for kmer in kmer_set:
    current_result = []
    current_result.append(kmer)
    for read in refbinreads:
        if kmer in read:
            current_result.append(1)
        else:
            current_result.append(0)
    refbinencode_result.append(current_result)  


# In[13]:


columnnames = []
for i in range(0,500):
    columnnames.append('bin'+str(i))
columnnames = ["kmer"] + columnnames


# In[14]:


refbindf = pd.DataFrame(refbinencode_result, columns=columnnames)


# In[15]:


pd.set_option('display.max_columns', None)
print(refbindf.head(5))


# # Step6

# In[16]:


def minhash(encode_result):
    readminhash = []
    for index in range(0,1000):
        np.random.seed(index)
        permuation = np.random.permutation(72530)
        curpermutation = []
        for bin in range(1,len(encode_result[0])):
            exist = 0
            for curindex in permuation:
                if encode_result[curindex][bin] == 1:
                    exist +=1
                    curpermutation.append(exist)
                    break
                else:
                    exist +=1
        readminhash.append(curpermutation)
    return readminhash


# In[17]:


def minhashpara(encode_result, loopnumber):
    readminhash = []
    for index in range(0,loopnumber):
     
        permuation = np.random.permutation(72530)
        curpermutation = []
        for bin in range(1,len(encode_result[0])):
            exist = 0
            for curindex in permuation:
                if encode_result[curindex][bin] == 1:
                    exist +=1
                    curpermutation.append(exist)
                    break
                else:
                    exist +=1
        readminhash.append(curpermutation)
    return readminhash


# In[18]:


def minhash_parallel(encode_result):
    def collect_result(result):
        global readminhash
        readminhash.append(result)
    with mp.Pool(20) as pool:
        for i in range(0,20):
            pool.apply_async(minhashpara, args=(encode_result, 50), callback=collect_result)
        pool.close()
        pool.join()
    return readminhash


# In[19]:


start_time = time.time()
readminhash_1 = minhash(encode_result)
end_time = time.time()
print("Sequential method time:", end_time - start_time, "seconds")


# In[20]:

readminhash = []
start_time = time.time()
readminhash_parallel = minhash_parallel(encode_result)
end_time = time.time()
print("Parallel method time:", end_time - start_time, "seconds")
result = []
for i in readminhash_parallel:
  result = result + i
readminhash_parallel = result


# In[21]:


start_time = time.time()
binminhash = minhash(refbinencode_result)
end_time = time.time()
print("Sequential method time:", end_time - start_time, "seconds")


# In[22]:

readminhash = []
start_time = time.time()
binminhash_parallel = minhash_parallel(refbinencode_result)
end_time = time.time()
print("Parallel method time:", end_time - start_time, "seconds")
result = []
for i in binminhash_parallel:
  result = result + i
binminhash_parallel = result

# # Step7

# In[23]:


jaccard = []
for readnum in range(0,2000):
    maxbinsim = 0
    for binnum in range(0,500):
        num = 0
        for rownum in range(0,1000):
            if readminhash_1[rownum][readnum] == binminhash[rownum][binnum]:
                num += 1
        cursim = num / 1000
        if cursim > maxbinsim:
            maxbinsim = cursim
            maxbinname = 'bin' + str(20000000 + 100 * binnum) + "_" + str(20000100 + 100 * binnum)
    jaccard.append(["read" + str(readnum), maxbinname, maxbinsim])


# In[24]:


print(jaccard)


# In[ ]:


jaccard_parallel = []
for readnum in range(0,2000):
    maxbinsim = 0
    for binnum in range(0,500):
        num = 0
        for rownum in range(0,1000):
            if readminhash_parallel[rownum][readnum] == binminhash_parallel[rownum][binnum]:
                num += 1
        cursim = num / 1000
        if cursim > maxbinsim:
            maxbinsim = cursim
            maxbinname = 'bin' + str(20000000 + 100 * binnum) + "_" + str(20000100 + 100 * binnum)
    jaccard_parallel.append(["read" + str(readnum), maxbinname, maxbinsim])


# In[ ]:


print(jaccard_parallel)


# # Step8

# In[25]:


# Read in the CSV file using Pandas and store reference_end in a list
df = pd.read_csv('read_position_benchmark.csv')
reference_end = df['reference_end'].tolist()


# In[26]:


# Extract the numbers after the _ in each list in jaccard
numbers = [int(x[1].split('_')[1]) for x in jaccard]

# Calculate the Pearson correlation between the two lists
correlation, p_value = pearsonr(reference_end, numbers)
print("Pearson correlation: ", correlation)


# In[ ]:


# Extract the numbers after the _ in each list in jaccard
numbers_parallel = [int(x[1].split('_')[1]) for x in jaccard_parallel]

# Calculate the Pearson correlation between the two lists
correlation_parallel, p_value_parallel = pearsonr(reference_end, numbers_parallel)
print("Pearson correlation: ", correlation_parallel)


# In[ ]:


## Draw a scatter plot to show the correlation
#plt.scatter(reference_end, numbers)
#plt.xlabel("Reference end")
#plt.ylabel("Bin end")
#plt.show()


# In[27]:


# Calculate the mean square error
mse = np.square(np.subtract(reference_end, numbers)).mean()
print("Mean square error: ", mse)


# In[ ]:


# Calculate the mean square error
mse_parallel = np.square(np.subtract(reference_end, numbers_parallel)).mean()
print("Mean square error: ", mse_parallel)

