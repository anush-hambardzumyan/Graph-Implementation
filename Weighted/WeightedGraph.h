#ifndef WEIGHTEDGRAPH_H
#define WEIGHTEDGRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <limits>
#include <utility>
#include <iomanip>
#include "DisjointSet.h"

class WeightedGraph {
    private:
    std::vector<std::vector<std::pair<int, int>>> adjList;
    std::vector<std::vector<int>> adjMatrix;

    public:
    WeightedGraph(int x = 0): adjList(x , std::vector<std::pair<int,int>>()), adjMatrix(x,std::vector<int>(x, std::numeric_limits<int>::infinity())) {}
    WeightedGraph(const WeightedGraph& rhv) : adjList(rhv.adjList), adjMatrix(rhv.adjMatrix) {}
    WeightedGraph(std::vector<std::vector<std::pair<int, int>>> List, std::vector<std::vector<int>> Matrix) : adjList(List), adjMatrix(Matrix) {}

    void addVertex();                                                  
    void addEdge(int vertex1, int vertex2, int weight, bool undirected);
    void printList() const;
    void printMatrix() const;
    void removeVertex(int vertex);                                     
    void removeEdge(int vertex1 ,int vertex2, bool undirected);        
    void getNeighbours(int vertex) const;                                    
    bool hasVertex(int vertex) const;                                           
    int vertexCount() const;                                                 
    int edgeCount() const ;
    std::vector<int> dfs(int src, bool recursive) const;
    std::vector<int> bfs(int src) const;
    std::vector<std::vector<int>> bfsExtraCase() const;
    std::vector<std::vector<int>> dfsExtraCase() const;
    int connectedComponentsCount() const;
    void transpose();
    WeightedGraph getTransposedCopy() const;
    std::vector<std::pair<int, std::vector<int>>> getAllPaths(int src, int dest) const; 
    std::vector<int> getNodesOfKthLevel(int k, bool dfs = true) const;  //grac chi
    bool isCycledUndirected() const;
    bool isCycledDirected() const;
    std::vector<int> topologicalSort() const;
    std::list<int> KahnsAlgorithm() const;
    std::vector<std::vector<int>> KosarajusAlgorithm() const; //test
    std::vector<int> SSSPTopSort(int src) const;  //test????
    std::vector<int> DijkstrasAlgorithm(int src) const;
    int PrimsAlgorithm(int src) const;
    int KruskalsAlgorithm() const;

    private:
    void dfsRecursive(int src, std::vector<bool>& visited, std::vector<int>& result) const;
    void dfsIterative(int src, std::vector<int>& result) const;
    std::vector<int> bfsHelper(int src, std::vector<bool>& visited) const;
    void allPathHelper(int src, int dest, std::pair<int, std::vector<int>>& currPath, std::vector<bool>& visited, std::vector<std::pair<int, std::vector<int>>>& allPaths) const;
    bool isCycledUndirectedHelper(int src, int parent,std::vector<bool>& visited) const;
    bool isCycledDirectedHelper(int src, std::vector<bool>& onStack, std::vector<bool>& visited) const;
    void topSortHelper(int src, std::vector<bool>& visited, std::stack<int>& res) const;
    void fillOrder(int src, std::vector<bool>& visited, std::stack<int>& stack) const;
};


#include "WeightedGraph.hpp"
#endif 