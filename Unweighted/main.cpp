#include <iostream>
#include "UnweightedGraph.h"

int main() 
{
    // Create an instance of UnweightedGraph
    UnweightedGraph graph(5);
    graph.addEdge(0, 2, false);
    graph.addEdge(0, 3, false);
    graph.addEdge(1, 0, false);
    graph.addEdge(2, 1, false);
    graph.addEdge(3, 4, false);

    auto res = graph.TarjansAlgorithm();
    for (int i = 0; i < res.size(); ++i) {
        for (int elem : res[i]) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    // UnweightedGraph graph(8);
    // graph.addEdge(0, 1, false);
    // graph.addEdge(1, 2, false);
    // graph.addEdge(2, 0, false);
    // graph.addEdge(2, 3, false);
    // graph.addEdge(3, 4, false);
    // graph.addEdge(4, 7, false);
    // graph.addEdge(6, 7, false);
    // graph.addEdge(6, 4, false);
    // graph.addEdge(4, 5, false);
    // graph.addEdge(5, 6, false);

    // auto res = graph.KosarajusAlgorithm();
    // for (int i = 0; i < res.size(); ++i) {
    //     for (int elem : res[i]) {
    //         std::cout << elem << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // UnweightedGraph graph(5);
    // graph.addEdge(0, 1, false);
    // graph.addEdge(2, 1, false);
    // graph.addEdge(2, 0, false);
    // graph.addEdge(1, 3, false);
    // graph.addEdge(3, 4, false);

    // std::cout << std::boolalpha << graph.isCycledDirected();
    //UnweightedGraph graph(6);

    // graph.addEdge(5, 0, false);
    // graph.addEdge(4, 0, false);
    // graph.addEdge(5, 2, false);
    // graph.addEdge(2, 3, false);
    // graph.addEdge(3, 1 ,false);

    // auto res = graph.topologicalSort();

    // for(auto it : res) {
    //     std::cout << it << " ";
    // }
    // graph.addEdge(0, 9);
    // graph.addEdge(0, 3);
    // graph.addEdge(0, 2);
    // graph.addEdge(0, 1);
    // graph.addEdge(1, 2);
    // graph.addEdge(1, 8);
    // graph.addEdge(2, 4);
    // graph.addEdge(3, 5);
    // graph.addEdge(4, 5);
    // graph.addEdge(5, 6);
    // graph.addEdge(7, 4);
    // graph.addEdge(7, 8);

    // graph.printList();
    //  auto path = graph.getAllPaths(0,5);

    // for(int i = 0; i < path.size(); ++i) {
    //     for (int j = 0 ; j < path[i].size(); ++j) {
    //         std::cout << path[i][j] << " -> ";
    //     } 
    //     std::cout << std::endl;
    // }

//    std::cout<<std::boolalpha << graph.isCycledUndirected();

    // auto path = graph.getNodesOfKthLevel(0);
    // std::cout << path; 
    //for (auto elem : path) std::cout << elem << " ";

    //std::cout << graph.count_connected_components();


    // // Add edges
    // graph.addEdge(0, 1);
    // graph.addEdge(0, 2);
    // graph.addEdge(1, 3);
    // graph.addEdge(1, 4);
    // graph.addEdge(2, 5);

    // // Print adjacency list
    // std::cout << "Adjacency List:" << std::endl;
    // graph.printList();
    // std::cout << std::endl;

    // // Print adjacency matrix
    // std::cout << "Adjacency Matrix:" << std::endl;
    // graph.printMatrix();
    // std::cout << std::endl;

    // // Perform DFS recursively
    // std::cout << "DFS Recursive (starting from vertex 0): ";
    // graph.dfs(0, true);

    // //Perform DFS iteratively
    // std::cout << "DFS Iterative (starting from vertex 0): ";
    // graph.dfs(0, false);

    //Perform BFS (only iterative)
    // std::cout << "BFS Iterative (starting from vertex 0): ";
    // graph.bfs(0);

    // auto vec = graph.getNodesOfKthLevel(2);
    // for (auto elem : vec) std::cout << elem << " ";
    // std::cout << std::endl;
    // return 0;
}
