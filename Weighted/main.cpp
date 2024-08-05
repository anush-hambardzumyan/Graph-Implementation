#include "WeightedGraph.h"

int main() 
{
    WeightedGraph graph(3);
    graph.addEdge(0, 1, 10, false);
    graph.addEdge(0, 2, 7, false);
    graph.addEdge(1, 2, 3, false);
    graph.printList();
    graph.printMatrix();
    std::cout << graph.edgeCount();
}