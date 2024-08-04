#ifndef UNWEIGHTEDGRAPH_H
#define UNWEIGHTEDGRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <list>


class UnweightedGraph {
    private:
    std::vector<std::vector<int>> adjList;
    std::vector<std::vector<bool>> adjMatrix;

    public:
    UnweightedGraph(int x = 0): adjList(x , std::vector<int>()), adjMatrix(x,std::vector<bool>(x, 0)) {}
    UnweightedGraph(const UnweightedGraph& rhv) : adjList(rhv.adjList), adjMatrix(rhv.adjMatrix) {}
    UnweightedGraph(std::vector<std::vector<int>> List, std::vector<std::vector<bool>> Matrix) : adjList(List), adjMatrix(Matrix) {}

    void addVertex();                                                  
    void addEdge(int vertex1 ,int vertex2, bool undirected);           
    void printList() const;                                                   
    void printMatrix() const;                                                 
    void removeVertex(int vertex);                                     
    void removeEdge(int vertex1 ,int vertex2, bool undirected);        
    void getNeighbours(int vertex) const;                                    
    bool hasVertex(int vertex) const;                                           
    int vertexCount() const;                                                 
    int edgeCount() const ;                                                   
    void dfs(int src, bool recursive) const;
    void dfsExtraCase(int vertex) const;                  //dfs for more than one components
    void bfs(int src) const;
    int countConnectedComponents(bool printDFS = false) const;
    void transpose();
    std::vector<std::vector<int>> getTransposedAdjList() const;
    std::vector<std::vector<bool>> getTransposedAdjMatrix() const;
    std::vector<int> getShortestPath(int src, int dest) const; 
    std::vector<std::vector<int>> getAllPaths(int src, int dest) const; 
    bool isCycledUndirected() const;
    bool isCycledDirected() const;
    std::list<int> topologicalSort() const;
    std::list<int> KahnsAlgorithm() const;
    std::vector<std::vector<int>> KosarajusAlgorithm() const;
    std::vector<std::vector<int>> TarjansAlgorithm() const;

    std::vector<int> getNodesOfKthLevel(int k, bool dfs = true) const; //??????

    private:
    void dfsInVector(int src, std::vector<bool>& visited, std::vector<int>& result) const; //for Kosaraju
    void TarjanHelper(int src, std::vector<int>& ids, std::vector<int>& lowlink, std::stack<int>& st, std::vector<bool>& onStack, std::vector<std::vector<int>>& SCCs) const;
    void dfsRecursive(int src, std::vector<bool>& visited, bool print = true) const;
    void dfsIterative(int src) const;
    void allPathHelper(int src, int dest, std::vector<int>& currPath, std::vector<bool>& visited, std::vector<std::vector<int>>& allPaths) const;
    bool isCycledUndirectedHelper(int src, int parent,std::vector<bool>& visited) const;
    bool isCycledDirectedHelper(int src, std::vector<bool>& onStack, std::vector<bool>& visited) const;
    void getNodesOfLevelDFSHelper(int src, int level, std::vector<int>& res, std::vector<bool>& visited) const;
    void topSortHelper(int src, std::vector<bool>& visited, std::list<int>& res) const;
    void fillOrder(int src, std::vector<bool>& visited, std::stack<int>& stack) const;
};

#include "UnweightedGraph.hpp"
#endif     