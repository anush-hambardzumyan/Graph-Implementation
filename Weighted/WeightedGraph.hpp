#include "WeightedGraph.h"

void WeightedGraph::addVertex()
{
    adjList.push_back(std::vector<std::pair<int, double>>());
    adjMatrix.push_back(std::vector<double>(adjList.size(), std::numeric_limits<double>::infinity()));
    for (int i = 0; i < adjMatrix.size() - 1; ++i) {
        adjMatrix[i].push_back(std::numeric_limits<double>::infinity());
    }
}

void WeightedGraph::addEdge(int vertex1, int vertex2, double weight, bool undirected = true)
{
    if(vertex1 < 0 || vertex1 >= adjList.size() || vertex2 < 0 || vertex2 >= adjList.size())
    {
        std::cout << "Unvalid vertices to add edge: " << std::endl;
        exit(0);
        return;
    }

    auto it1 = std::find(adjList[vertex1].begin() , adjList[vertex1].end() , vertex2 );
    if(it1 == adjList[vertex1].end())
    {
        adjList[vertex1].push_back(vertex2);
        adjMatrix[vertex1][vertex2] = true;
    } 

    
    if (undirected) {
        auto it2 = std::find(adjList[vertex2].begin() , adjList[vertex2].end() , vertex1 );
        if (it2 == adjList[vertex2].end()) {
            adjList[vertex2].push_back(vertex1);
            adjMatrix[vertex2][vertex1] = true;
        }
    }
}
}