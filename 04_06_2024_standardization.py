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