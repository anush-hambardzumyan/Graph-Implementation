#include "WeightedGraph.h"

int main() 
{

    int V = 5;
    WeightedGraph graph(V);

    graph.addEdge(0, 1, 2);
    graph.addEdge(0, 3, 3);
    graph.addEdge(1, 3, 4);
    graph.addEdge(1, 4, 8);
    graph.addEdge(3, 4, 7);

    std::cout << "Prim's Algorithm MST Weight: " << graph.PrimsAlgorithm(0) << std::endl;
    std::cout << "Kruskal's Algorithm MST Weight: " << graph.KruskalsAlgorithm() << std::endl;
    // WeightedGraph graph(7);
    // graph.addEdge(1,4,10, false);
    // graph.addEdge(1,2,50,0);
    // graph.addEdge(1,3,45,0);
    // graph.addEdge(2, 4, 15,0);
    // graph.addEdge(2, 3, 10, 0);
    // graph.addEdge(3, 5, 30, 0);
    // graph.addEdge(4, 1, 10, 0);
    // graph.addEdge(4, 5, 15, 0);
    // graph.addEdge(5, 2, 10, 0);
    // graph.addEdge(5, 3, 35, 0);
    // graph.addEdge(6, 5, 3, 0);

    // graph.printList();
    // graph.printMatrix();
    // auto res = graph.DijkstrasAlgorithm(1);
    // for (auto elem : res) {
    //     std::cout << elem << " ";
    // }
    
}