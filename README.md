# Classification of Acyclic Fake Surfaces
This repository contains the data for the classification of acyclic cellular fake surfaces of complexity 1-4 and a partial classification of complexity 5. It also contains the  code used to generate the classification and check contractibility. 

## Fake Surface Data

All acyclic cellular fake surfaces of complexity 1-4 can be found in ``fakesurfaces.csv.`` As described in [], the convention used is as follows: 

Each row represents an acyclic cellular fake surface. The two columns represents the complexity followed the number of the 1-skeleton. Then, all the disk attaching maps are listed, where a negative entry denotes going along the edge in the opposite way. Following each disk, the two `Y/N` columns are the answers to whether the disk is embedded and has a trivial $T$-bundle, respectively.  

Thus, our notation of fake surfaces relies on ordering the vertices and edges within a 1-skeleton and also ordering the 1-skeletons of a given complexity. We do this in the following way.

Viewing the 1-skeletons as 4-regular multigraphs, consider the map
$\mathsf{D}:M_{n}(\mathbb{Z})\to \mathbb{Z}$ that gives a decimal representation of a matrix by simply concatenating the rows from top to bottom into an integer with $n^2$ digits. An example makes this clear: 
$$ \mathsf{D}\left(\begin{bmatrix}1&2\\3&4\end{bmatrix}\right)=1234.$$
Explicitly, if $A=(a_{i,j})$, $$ \mathsf{D}(A)=\sum_{i, j=1}^n 10^{n(n-i)+n-j} a_{i,j}$$ as an integer.

For a given 1-skeleton of complexity $n$, we choose the adjacency matrix $A=(a_{i,j})$ and therefore an ordering of the vertices that maximizes $\mathsf{D}(A)$. Note this also provides the canonical ordering of the 1-skeleta for a given complexity: from largest value of $\mathsf{D}(A)$ to smallest. 

We then label the edges in the order they appear in the upper-right triangle of the chosen adjacency matrix, read left-to-right and top-to-bottom. To be precise, using our ordering of the vertices, let $(i_1,j_1)$ and $(i_2,j_2)$ be edges satisfying $i_1\leq j_1$ and $i_2\leq j_2$. Then $(i_1,j_1)<(i_2,j_2)$ if $i_1<i_2$ or $i_1=i_2$ and $j_1<j_2$. We then label edges from 1 to $2n$ in accordance with this ordering. The order of multiple edges between two vertices can be chosen arbitrarily. 

We provide an example to make this clear. Consider the graph $\Gamma$ given by the adjacency matrix 
$$A=
\begin{bmatrix}
2&1&1\\
1&0&3\\
1&3&0
\end{bmatrix}.
$$
Observe that this ordering of the vertices is the one that maximizes $\mathsf{D}(A)$. Then the self-loop is first edge, the two other edges at that vertex are the next two edges, and the three edges connecting the other two vertices are the fourth, fifth, and sixth edges. Thus we label the edges by any of the 12 valid labelings. One such labeling is shown in Fig.~\ref{fig:example}:
\begin{figure}[h]
\caption{Valid Edge Labeling for $\Gamma$}
\label{fig:example}
\centering
\begin{tikzpicture}[
    node distance=1.5cm, 
    every node/.style={circle, draw, fill=black, inner sep=3pt},
    edge_label/.style={midway, draw=none, inner sep=2pt, fill=none} 
] 
    \node (A) {};
    \node (B) [below right of=A, yshift=-0.6cm] {};
    \node (C) [below left of=A, yshift=-0.6cm] {};

    \draw (A) to node[edge_label, above right] {3} (B);
    \draw (B) to[bend left=40] node[edge_label, below] {6} (C); 
    \draw (B) to node[edge_label] {5} (C); 
    \draw (B) to[bend right=40] node[edge_label, above] {4} (C);
    \draw (A) to node[edge_label, above left] {2} (C); 
    \draw (A) to [out=135, in=45, looseness=10] node[edge_label, above] {1} (A); 
\end{tikzpicture}
\end{figure}


## Classification Code
The code that generated the classification can be found in ``fakesurfaces_cla_6.py.`` In the filename, `cla` refers to the fact that it takes as input a command line argument corresponding to the 1-skeleton, and `6` that all 1-skeleta up to complexity 6 are written into the file.

## Checking Contractibility
For a given surface and maximal tree, we can put them at the top of this file as shown and then run ``sage check_contractibility.sage.`` It is also easy to modify the code to check multiple surfaces at once.  

## Generating One-Skeletons
We provide the code used to generate the one-skeleta in ``generate_one_skeleta.py``. Here we specify a complexity at the top, called ``SIZE_OF_MATRIX`` and then use the functions as named. Sample usage is provided at the bottom of the script.   
