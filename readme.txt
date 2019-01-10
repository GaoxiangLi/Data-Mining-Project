coresets.py - This program generates a lightweight coreset using the algorithm from Bachem, Lucic, Krause "Scalable k-Means Clustering via Lightweight Coresets"

To run this code, python version must be greater than Python 3.6

usage: python coresets.py dataset_filename m
	dataset_filename: 	filename of dataset to be converted
	m: 			number of subsamples in lightweight coreset
example: python coresets.py bio_train.dat 1000

In the experiments done in the paper, the values of m used were: {1000, 2000, 5000, 10000, 20000}

Once the lightweight coresets are computed, they are exported to a file called "export.dat" in the same directory.
The resulting coresets are evaluated in a separate k-means++ algorithm and compared for accuracy and performance