#include "WeightedGraph.h"

void WeightedGraph::addVertex()
{
    adjList.push_back(std::vector<std::pair<int, double>>());
    adjMatrix.push_back(std::vector<double>(adjList.size(), std::numeric_limits<double>::infinity()));
    for (int i = 0; i < adjMatrix.size() - 1; ++i) {
        adjMatrix[i].push_back(std::numeric_limits<double>::infinity());
    }
}

void WeightedGraph::addEdge(int vertex1, int vertex2, double weight, bool undirected)
{
    if (vertex1 < 0 || vertex1 >= adjList.size() || vertex2 < 0 || vertex2 >= adjList.size()) {
        std::cout << "Invalid vertices to add edge: " << std::endl;
        exit(0);
    }

    auto it1 = std::find_if(adjList[vertex1].begin(), adjList[vertex1].end(),
                            [vertex2](const std::pair<int, double>& edge) { return edge.first == vertex2; });
    if (it1 == adjList[vertex1].end()) {
        adjList[vertex1].push_back(std::make_pair(vertex2, weight));
    } else {
        it1 -> second = weight;
    }
    adjMatrix[vertex1][vertex2] = weight;

    if (undirected)
    {
        auto it2 = std::find_if(adjList[vertex2].begin(), adjList[vertex2].end(),
                                [vertex1](const std::pair<int, double>& edge) { return edge.first == vertex1; });
        if (it2 == adjList[vertex2].end())
        {
            adjList[vertex2].push_back(std::make_pair(vertex1, weight));
        } else {
            it2 -> second = weight;
        }
        adjMatrix[vertex2][vertex1] = weight;
    }
}

void WeightedGraph::printList() const
{
    for (int i = 0; i < adjList.size(); ++i) {
        std::cout << i << " -> ";
        for (int j = 0; j < adjList[i].size(); ++j) {
            std::cout << "{ " <<adjList[i][j].first << ", " <<adjList[i][j].second << " }  ";
        }
        std::cout << std::endl;
    }
}

void WeightedGraph::printMatrix() const
{
    const char* infinity = "inf";
    const int width = 6; 

    std::cout << std::setw(width) << " ";
    for (int i = 0; i < adjMatrix.size(); ++i) {
        std::cout << std::setw(width) << i;
    }
    std::cout << std::endl;

    for (int i = 0; i < adjMatrix.size(); ++i) {
        std::cout << std::setw(width) << i; 
        for (int j = 0; j < adjMatrix[i].size(); ++j) {
            if (adjMatrix[i][j] == std::numeric_limits<double>::infinity()) {
                std::cout << std::setw(width) << infinity;
            } else {
                std::cout << std::setw(width) << adjMatrix[i][j];
            }
        }
        std::cout << std::endl;
    }
}

void WeightedGraph::removeVertex(int vertex)
{
    if(vertex < 0 || vertex > adjList.size() - 1)
    {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
        return;
    }

    adjList.erase(adjList.begin() + vertex);
    adjMatrix.erase(adjMatrix.begin() + vertex);
    for (int i = 0; i < adjMatrix.size(); ++i) {
        adjMatrix[i].erase(adjMatrix[i].begin() + vertex);
    }

    for (int i = 0; i < adjList.size(); ++i) 
    {
        for (int j = 0; j < adjList[i].size(); ++j) 
        {
            if (adjList[i][j].first == vertex) 
            {
                adjList[i].erase(adjList[i].begin() + j);
            }
        }   
    }   
}

void WeightedGraph::removeEdge(int vertex1 ,int vertex2, bool undirected)
{
    if (vertex1 < 0 || vertex1 >= adjList.size() || vertex2 < 0 || vertex2 >= adjList.size())
    {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
        return;
    }

    for (int i = 0; i < adjList[vertex1].size(); ++i)
    {
        if (adjList[vertex1][i].first == vertex2)
        {
            adjList[vertex1].erase(adjList[vertex1].begin() + i);
            return;
        }
    }
    adjMatrix[vertex1][vertex2] = std::numeric_limits<int>::infinity();

    if(undirected) {
        for (int i = 0; i < adjList[vertex2].size(); ++i)
        {
            if (adjList[vertex2][i].first == vertex1)
            {
                adjList[vertex2].erase(adjList[vertex2].begin() + i);
                return;
            }
        }
        adjMatrix[vertex2][vertex1] = std::numeric_limits<int>::infinity();
    }
}        

void WeightedGraph::getNeighbours(int vertex) const
{
    if(vertex < 0 || vertex >= adjList.size())
    {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
        return; 
    }

    std::cout << "Neighbours of vertex are: "; 
    for(int i = 0; i < adjList[vertex].size(); ++i)
    {
        std::cout << adjList[vertex][i].first << " ";
    }

    std::cout << std::endl;
}

bool WeightedGraph::hasVertex(int vertex) const
{
    return !(vertex < 0 || vertex >= adjList.size());
}

int WeightedGraph::vertexCount() const
{
    return adjList.size();
}

int WeightedGraph::edgeCount() const 
{
    if(!adjList.size())
    {
        return 0;
    }

    int count = 0;
    for(int i = 0; i < adjList.size(); ++i)
    {
        count += adjList[i].size();
    }
    return count;
}

std::vector<int> WeightedGraph::dfs(int src, bool recursive) const
{
    if (src < 0 || src >= vertexCount()) {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
    }

    std::vector<int> result;    
    if (recursive) {
        std::vector<bool> visited(vertexCount(), false);
        dfsRecursive(src, visited, result);
    } else dfsIterative(src, result);

    return result;
}

void WeightedGraph::dfsRecursive(int src, std::vector<bool>& visited, std::vector<int>& result) const
{
    visited[src] = true;
    result.push_back(src);

    
    for (auto elem : adjList[src]) {
        if (!visited[elem.first]) {
            dfsRecursive(elem.first, visited, result);
        }
    }
}

void WeightedGraph::dfsIterative(int src, std::vector<int>& result) const
{
    std::vector<bool> visited(adjList.size(), false);
    std::stack<int> st;

    visited[src] = true;
    st.push(src);

    while(!st.empty()) {
        int elem = st.top();
        result.push_back(elem);
        st.pop();

        for (auto neighbor : adjList[elem]) {
            if (!visited[neighbor.first]) {
                visited[neighbor.first] = true;
                st.push(neighbor.first);
            }
        }
    }

}

std::vector<int> WeightedGraph::bfs(int src) const
{
    if (src < 0 || src >= vertexCount()) {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
    }

    std::vector<bool> visited(adjList.size(), false);
    
    return bfsHelper(src, visited);
}

std::vector<int> WeightedGraph::bfsHelper(int src, std::vector<bool>& visited) const
{
    std::queue<int> qu;
    std::vector<int> result;
    visited[src] = true;
    qu.push(src);

    while (!qu.empty()) {
        int elem = qu.front();
        result.push_back(elem);
        qu.pop();

        for (auto neighbor : adjList[elem]) {
            if (!visited[neighbor.first]) {
                visited[neighbor.first] = true;
                qu.push(neighbor.first);
            }
        }
    }
    return result;
}

std::vector<std::vector<int>> WeightedGraph::dfsExtraCase() const
{
    std::vector<std::vector<int>> result; 
    std::vector<bool> visited(vertexCount(), false);
    for (int i = 0; i < vertexCount(); ++i) {
        if (!visited[i]) {
            std::vector<int> currComponent;
            dfsRecursive(i, visited, currComponent);
            result.push_back(currComponent);
        }
    }

    return result;
}

std::vector<std::vector<int>> WeightedGraph::bfsExtraCase() const
{
    std::vector<std::vector<int>> result; 
    std::vector<int> currComponent;
    std::vector<bool> visited(vertexCount(), false);
    for (int i = 0; i < visited.size(); ++i) {
        if (!visited[i]) {
            result.push_back(bfsHelper(i, visited));
        }
    }

    return result;
}

int WeightedGraph::connectedComponentsCount() const
{
    return dfsExtraCase().size();
}

