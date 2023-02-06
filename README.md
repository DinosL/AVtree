# Adaptive Indexing in High-Dimensionl Metric Spaces

### Compile
Compile using ```make all distFunc=<option>``` where ```<option>``` can be one of the following distance functions:
   - L1, for Manhattan Distance
   - L2, for Euclidean Distance
   - ED, for Edit Distance
   
Two executable files are generated; ```range``` for range queries and ```knn``` for kNN queries. See examples below for correct usage.
   

### Parameters
We tried to be as consistent as possible when it comes to parameters:
| Parameters | README |
| ------ | ------ |
| -l | Linear scan. When used for kNN queries, use ```-n``` to set the number of neighbors. |
| -s | Standard version. When used for kNN queries, use ```-n``` to set the number of neighbors. |
| -m | Standard version with mediocre cracking. When used for kNN queries, use ```-n``` to set the number of neighbors. |
| -c | Standard version with mediocre cracking and Caching |
| -t | Threshold value to be used |
| -n | Number of neighbors for kNN queries |

Note that for kNN queries only methods ```-l, -s, -m``` are available.

### Files

- ##### tree/AVtree.h and tree/AVtreeStrings.h: 
    
    Contains all the required methods for the tree index structure.  

- ##### main_range.cpp: 
    
    Contains the code used for the range query experiments.
    ##### Example 
    - ###### Linear scan
    
        ```sh
        $ ./range -l datasets/mnist50D_data.txt datasets/mnist50D_queries.txt
        ```
    - ###### Standard mediocre with caching with 128 threshold
    
        ```sh
        $ ./range -c -t 128 datasets/mnist50D_data.txt datasets/mnist50D_queries.txt 
        ```

- ##### main_knn.cpp: 
    
    Contains the code used for the kNN query experiments.
    ##### Example 
    - ###### Linear scan, 5 nearest neighbors
    
        ```sh
        $ ./knn -l -n 5 datasets/mnist50D_data.txt datasets/mnist50D_queries.txt 
        ```
    - ###### Standard mediocre with caching with 128 threshold, 5 nearest neighbors
    
        ```sh
        $ ./knn -m -t 128 -n 5 datasets/mnist50D_data.txt datasets/mnist50D_queries.txt
        ```
