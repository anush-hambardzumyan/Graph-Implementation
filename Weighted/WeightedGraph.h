#ifndef WEIGHTEDGRAPH_H
#define WEIGHTEDGRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <limits>
#include <utility>

class WeightedGraph {
    private:
    std::vector<std::vector<std::pair<int, double>>> adjList;
    std::vector<std::vector<double>> adjMatrix;

    public:
    WeightedGraph(int x = 0): adjList(x , std::vector<std::pair<int,double>>()), adjMatrix(x,std::vector<double>(x, std::numeric_limits<double>::infinity())) {}
    WeightedGraph(const WeightedGraph& rhv) : adjList(rhv.adjList), adjMatrix(rhv.adjMatrix) {}
    WeightedGraph(std::vector<std::vector<std::pair<int, double>>> List, std::vector<std::vector<double>> Matrix) : adjList(List), adjMatrix(Matrix) {}

    void addVertex();                                                  
    void addEdge(int vertex1, int vertex2, double weight, bool undirected);
    void printList() const;
    void printMatrix() const;
    void removeVertex(int vertex);                                     
    void removeEdge(int vertex1 ,int vertex2, bool undirected);        
    void getNeighbours(int vertex) const;                                    
    bool hasVertex(int vertex) const;                                           
    int vertexCount() const;                                                 
    int edgeCount() const ;
    
};
#endif 