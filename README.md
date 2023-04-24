# Advertisement-Management-System
This is a simple advertisement system where we parallely identify the influencers which are connected to maximum number of k-truss groups from a social network presented in the form of a graph. We have implemented the truss-computation parallely by dividing the workload among different processers to improve efficiency. Refer to COL380 Assignment 2 for more explanation on truss decomposition and influencers.

## Task 1
Given a graph G, find the k-truss of G for a range of k in (startK, endK). A k-truss of G is a subgraph of G such that every edge in the subgraph is contained in at least k-2 triangles in G. A triangle is a set of three nodes that are connected by three edges. A k-truss is a maximal k-truss, i.e., a k-truss that cannot be extended to a larger k-truss.

## Task 2
Given a graph G, find the maximum k-truss of the graph for a given k. After that, we find the influencers in the graph. An influencer is a node that is connected to the atleast p number of k-truss groups in the graph. 

## Build Instructions 
1. Clone the repository
2. Change the Makefile to point to the correct directories for the following variables:
    * specify inputpath for the input file
    * specify outputpath for the output file
    * specify the headerpath for the header file
    * specify verbose to ensure the type of output. Refer to problem statement for more details.
3. Run the following command to build the project
    * make
4. Run the following command to run the project
    * make run

