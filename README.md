# Big-Data-Scaling-Final-Project

# Outline
- [Introduction](#introduction)
- [Method](#method)
- [Results](#results)
- [Conclusion](#conclusion)
- [Contributors](#contributors)


# Introduction

Bioinformatics is a field that involves the application of computational techniques and tools to analyze biological data, and one of the common challenges faced in this field is the alignment of reads to the reference genome. This task is particularly difficult due to several factors such as the large size of the reference genome, point mutations, and sequencing errors present in the reads. These issues can result in significant computational costs and time-consuming alignment processes. To address these challenges, researchers have turned to innovative solutions such as minhash. Minhash is a technique that can be used to align second-generation sequencing reads to the reference genome, and it has been found to be effective in reducing computational costs and speeding up the alignment process. The primary objective of this project is to leverage minhash to identify the most similar bin in the reference genome for each read and record the most similar bin for each read as well as the similarity measure.

# Method

- Kmers
- One-Hot Encoding
- Minhash
- Jaccard Similarity
- Parallel Processing

# Results

The reference bin with the highest Jaccard similarity for each read: 

<img width="232" alt="Screen Shot 2023-04-24 at 9 18 24 PM" src="https://user-images.githubusercontent.com/89117508/234158189-d01a43ec-b2b8-483d-8a7d-e58ac7efa276.png">

Comparing the difference in running time between sequential programming and parallel programming: 

<img width="631" alt="Screen Shot 2023-04-24 at 11 20 45 PM" src="https://user-images.githubusercontent.com/89117508/234173710-412b736f-cfb3-453f-b091-83ef32e55e44.png">


By further evaluating our results against a benchmark to ensure their accuracy. To do so, we compared the end position of the reference bin from our approach to the end position of the reference bin from the benchmark dataset. We used the Pearson correlation coefficient to measure the strength of the correlation between the two sets of data, and found a very high correlation coefficient of 0.99999. This means that our approach produced results that are strongly correlated with the benchmark data. We also calculated the mean square error between the two sets of data, which gave us a value of around 898.1. While this value may seem high at first glance, it falls within the range of our expectations.

<img width="585" alt="Screen Shot 2023-04-24 at 9 13 00 PM" src="https://user-images.githubusercontent.com/89117508/234158700-dce1fc8d-38d1-4740-84b2-d3b6478f3c73.png">



# Conclusion

In conclusion, our approach of using minhash to determine the most similar reference bin for each read has shown great promise for future applications. By calculating the Jaccard similarity between the reference bins and each read, we were able to identify the reference bin with the highest similarity score, and our results were strongly correlated with the benchmark data.

Overall, our approach holds great potential for various applications, including bioinformatics, data mining, and information retrieval. The results of this study can be used to guide further research in these areas and contribute to the development of more efficient and accurate algorithms for data analysis.


# Contributors

- Tinglei Wu (tinglei.wu@vanderbilt.edu)
- Minghao Zhao (minghao.zhao@vanderbilt.edu)
- Jingting Xu (jingting.xu@vanderbilt.edu)

