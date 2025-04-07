question = """
Write a function that takes a list of gene expression vectors and
normalizes each gene independently using z-score normalization.

Each inner list represents one gene, and each number is the gene’s 
expression level in a specific sample. If a gene has zero standard 
deviation (i.e., all values are identical), return a vector of zeros for that gene. 
If the input list is empty, return an empty list.

This step is a direct pre-processing step for clustering algorithms like K-means. 
Note: Outputs assume exact arithmetic; small floating-point deviations are acceptable.
"""

def normalize_expression(data: list[list[float]]) -> list[list[float]]:
    # Your code here
    if not data:
        return []
    rows = []
    for row in data:
        n = len(row)
        s = sum(row)
        mean = s/n
    
        std = 0
        for ele in row:
            std += (ele-mean)**2
        rows.append([])
        std = (std/(n-1)) ** 0.5
        if std == 0:
            rows[-1].extend([0 for _ in range(len(row))])
        else:
            for ele in row:
                rows[-1].append((ele-mean)/std)
    return rows

question_hard = """
Implement a basic version of the K-Means clustering algorithm 
to group genes by similarity in their expression patterns.

You must not use libraries like sklearn — write it from scratch using 
just Python’s built-in types and the math module. 
Initialize centroids by selecting the first k genes from the data. 
If the input data is empty, return an empty list. 
If k exceeds the number of genes, raise a ValueError. 
For best results, consider preprocessing the data with z-score normalization (e.g., from the Easy Challenge). 
The algorithm should terminate after max_iters iterations even if clusters haven’t fully converged.

Tip: Use Euclidean distance to measure similarity between genes and centroids. 
Cluster assignments may vary slightly due to initialization or iteration limits, 
but the grouping should reflect expression similarity.
"""
def euclidean_distance(gene, centroid):
    s = 0.0
    for (g, c) in zip(gene, centroid):
        s += (g-c)**2
    s = s ** 0.5
    return s

def get_closest_center(data, centroids):
    nearest_centroid = []
    for gene in data:
        nearest_centroid.append(-1)
        min_distance = float('inf')
        for ind, centroid in enumerate(centroids):
            distance = euclidean_distance(gene, centroid)
            if distance < min_distance:
                min_distance = distance
                nearest_centroid[-1] = ind
    return nearest_centroid

def calc_centroids(data, clusters, k):
    centroids = [[0 for _ in range(len(data))] for _ in range(k)]
    ns = [clusters.count(i) for i in range(k)]
    for gene, cluster_id in zip(data, clusters):
        for ind, g in enumerate(gene):
            centroids[cluster_id][ind] += g/ns[cluster_id]
    return centroids

def kmeans_gene_clustering(data: list[list[float]], k: int, max_iters: int = 100) -> list[int]:
    # Your code here
    if not data:
        return []
    if k > len(data):
        raise ValueError("k cannot be greater than number of genes")
    elif k == len(data):
        return [i for i in range(k)]
    else:
        data = normalize_expression(data)
        centroids = data[:k]
        for iter_num in range(max_iters):
            clusters = get_closest_center(data, centroids)
            centroids = calc_centroids(data, clusters, k)
            
        clusters = get_closest_center(data, centroids)
        return clusters

arr = [[5.1, 4.9, 5.0], [6.0, 6.1, 5.9], [1.0, 1.1, 0.9], [5.2, 5.0, 5.1], [0.8, 1.2, 1.0]]
arr = normalize_expression(arr)
from sklearn.cluster import KMeans

model = KMeans(n_clusters=2, max_iter=100)

clusters = model.fit_predict(arr)
print(clusters)
print(kmeans_gene_clustering(arr, 2))