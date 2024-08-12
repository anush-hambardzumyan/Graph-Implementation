#include "WeightedGraph.h"

void WeightedGraph::addVertex()
{
    adjList.push_back(std::vector<std::pair<int, int>>());
    adjMatrix.push_back(std::vector<int>(adjList.size(), std::numeric_limits<int>::max()));
    for (int i = 0; i < adjMatrix.size() - 1; ++i) {
        adjMatrix[i].push_back(std::numeric_limits<int>::max());
    }
}

void WeightedGraph::addEdge(int vertex1, int vertex2, int weight, bool undirected = false)
{
    if (vertex1 < 0 || vertex1 >= adjList.size() || vertex2 < 0 || vertex2 >= adjList.size()) {
        std::cout << "Invalid vertices to add edge: " << std::endl;
        exit(0);
    }

    auto it1 = std::find_if(adjList[vertex1].begin(), adjList[vertex1].end(),
                            [vertex2](const std::pair<int, int>& edge) { return edge.first == vertex2; });
    if (it1 == adjList[vertex1].end()) {
        adjList[vertex1].push_back(std::make_pair(vertex2, weight));
    } else {
        it1 -> second = weight;
    }
    adjMatrix[vertex1][vertex2] = weight;

    if (undirected)
    {
        auto it2 = std::find_if(adjList[vertex2].begin(), adjList[vertex2].end(),
                                [vertex1](const std::pair<int, int>& edge) { return edge.first == vertex1; });
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
            if (adjMatrix[i][j] == std::numeric_limits<int>::max()) {
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
    adjMatrix[vertex1][vertex2] = std::numeric_limits<int>::max();

    if(undirected) {
        for (int i = 0; i < adjList[vertex2].size(); ++i)
        {
            if (adjList[vertex2][i].first == vertex1)
            {
                adjList[vertex2].erase(adjList[vertex2].begin() + i);
                return;
            }
        }
        adjMatrix[vertex2][vertex1] = std::numeric_limits<int>::max();
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

std::vector<int> WeightedGraph::dfs(int src, bool recursive = true) const
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

void WeightedGraph::transpose()
{
    std::vector<std::vector<std::pair<int, int>>> transposedList(adjList.size());
    for (int i = 0; i < adjList.size(); ++i) {
        for (auto edge : adjList[i]) {
            int neighbor = edge.first;
            int weight = edge.second;
            transposedList[neighbor].push_back(std::make_pair(i, weight));
        }
    }
    adjList = transposedList;

    std::vector<std::vector<int>> transposedMatrix (adjMatrix.size(), std::vector<int>(adjMatrix.size(), std::numeric_limits<int>::max()));
    for (int i = 0; i < adjMatrix.size(); ++i) {
        for (int j = 0; j < adjMatrix[i].size(); ++j) {
            transposedMatrix[j][i] = adjMatrix[i][j];
        }
    }

    adjMatrix = transposedMatrix;
}

WeightedGraph WeightedGraph::getTransposedCopy() const
{
    WeightedGraph transposed(*this);
    transposed.transpose();
    return transposed;
}

std::vector<std::pair<int, std::vector<int>>> WeightedGraph::getAllPaths(int src, int dest) const
{
    if (src < 0 || src > adjList.size() - 1 || dest < 0 || dest > adjList.size() - 1) {
        std::cerr << "Error: Source or destination index is out of bounds." << std::endl;
        return {};
    }

    int n = adjList.size();
    std::vector<bool> visited(n, false);
    std::vector<std::pair<int, std::vector<int>>> allPaths;
    std::pair<int, std::vector<int>> currPath;
    allPathHelper(src, dest, currPath, visited, allPaths);
    return allPaths;
} 

void  WeightedGraph::allPathHelper(int src, int dest, std::pair<int, std::vector<int>>& currPath, std::vector<bool>& visited, std::vector<std::pair<int, std::vector<int>>>& allPaths) const
{
    visited[src] = true;
    currPath.second.push_back(src);

    if(src == dest) {
        allPaths.push_back(currPath);
    } else {
        for (auto i : adjList[src]) {
            if (!visited[i.first]) {
                currPath.first += i.second;
                allPathHelper(i.first, dest, currPath, visited, allPaths);  
                currPath.first -= i.second;
            }
        }
    }
    
    currPath.second.pop_back();
    visited[src] = false;
}

bool WeightedGraph::isCycledUndirected() const
{
    std::vector<bool> visited(adjList.size(), false);
    int parent = -1;
    for (int i = 0; i < adjList.size(); ++i) {
        if (!visited[i]) {
            if (isCycledUndirectedHelper(i, parent, visited)) {
                return true;
            }
        }
    }
    return false;
}

bool WeightedGraph::isCycledUndirectedHelper(int src, int parent, std::vector<bool>& visited) const
{
    visited[src] = true;

    for (const auto& neighbor : adjList[src]) {
        if (!visited[neighbor.first]) {
            if (isCycledUndirectedHelper(neighbor.first, src, visited)) {
                return true;
            }
        } else if (parent != neighbor.first) {
            return true;
        }
    }

    return false;
}

bool WeightedGraph::isCycledDirected() const
{
    std::vector<bool> onStack(adjList.size(), false);
    std::vector<bool> visited(adjList.size(), false);
    for (int i = 0; i < adjList.size(); ++i) {
        if (!visited[i]) {
            if (isCycledDirectedHelper(i, onStack, visited)) {
                return true;
            }
        }
    }
    return false;
}

bool WeightedGraph::isCycledDirectedHelper(int src, std::vector<bool>& onStack, std::vector<bool>& visited) const
{
    visited[src] = true;
    onStack[src] = true;

    for (const auto& neighbor : adjList[src]) {
        if (!visited[neighbor.first]) {
            if (isCycledDirectedHelper(neighbor.first, onStack, visited)) {
                return true;
            }
        } else if (onStack[neighbor.first]) {
            return true;
        }
    }

    onStack[src] = false;
    return false;
}

std::vector<int> WeightedGraph::topologicalSort() const
{
    if (isCycledUndirected()) {
        std::cerr << "Graph is cycled. Topological Sort is impossible. " << std::endl;
        return {};
    }
    std::vector<bool> visited(adjList.size(), false);
    std::stack<int> st;
    std::vector<int> result;
    for (int i = 0; i < adjList.size(); ++i) {
        if (!visited[i]) {
            topSortHelper(i, visited, st);
        }
    }

    while (!st.empty()) {
        result.push_back(st.top());
        st.pop();
    }

    std::reverse(result.begin(), result.end());
    return result;
}

void WeightedGraph::topSortHelper(int src, std::vector<bool>& visited, std::stack<int>& res) const
{
    visited[src] = true;

    for (const auto& neighbor : adjList[src]) {
        if (!visited[neighbor.first]) {
            topSortHelper(neighbor.first, visited, res);
        }
    }

    res.push(src);
}

std::list<int> WeightedGraph::KahnsAlgorithm() const
{
    std::vector<int> inDegree(adjList.size(), 0);
    for (int i = 0; i < adjList.size(); ++i) {
        for (const auto& elem : adjList[i]) {
            ++inDegree[elem.first];
        }
    }

    std::queue<int> q;
   
    for (int i = 0; i < inDegree.size(); ++i) {
        if (inDegree[i] == 0) {
            q.push(i);
        }
    } 

    std::list<int> result;
    while (!q.empty()) {
        int elem = q.front();
        q.pop();
        result.push_front(elem);

        for (const auto& neighbor : adjList[elem]) {
            --inDegree[neighbor.first];
            if (inDegree[neighbor.first] == 0) {
                q.push(neighbor.first);
            }
        }
    }

    if (result.size() != adjList.size()) {
        return {}; //cycled graph
    }

    return result;
}

void WeightedGraph::fillOrder(int src, std::vector<bool>& visited, std::stack<int>& stack) const
{
    topSortHelper(src, visited, stack);
}

std::vector<std::vector<int>> WeightedGraph::KosarajusAlgorithm() const
{
    std::vector<bool> visited(adjList.size(), false);
    std::stack<int> stack;
    for (int i = 0; i < visited.size(); ++i) {
        if (!visited[i]) {
            fillOrder(i, visited, stack);
        }
    }

    WeightedGraph transposedGraph = this -> getTransposedCopy();

    visited = std::vector<bool>(adjList.size(), false);
    std::vector<std::vector<int>> SCCs;
    
    while (!stack.empty()) {
        int elem = stack.top();
        stack.pop();
        if (!visited[elem]) {
            SCCs.push_back(transposedGraph.dfs(elem));
        }   
    }

    return SCCs;
}

std::vector<int> WeightedGraph::SSSPTopSort(int src) const {
    std::stack<int> topologicalsorting;
    std::vector<bool> visited(adjList.size(), false);
    
    for (int i = 0; i < adjList.size(); ++i) {
        if (!visited[i]) {
            topSortHelper(i, visited, topologicalsorting);
        }
    }

    std::vector<int> distances(adjList.size(), std::numeric_limits<int>::max());
    distances[src] = 0; // as a startpoint

    while (!topologicalsorting.empty()) {
        int elem = topologicalsorting.top();
        topologicalsorting.pop();

        if (distances[elem] != std::numeric_limits<int>::max()) { //for not analyzing unnecessary vertices
            for (const auto& neighbor : adjList[elem]) {
                int neighborIndex = neighbor.first;
                int edgeWeight = neighbor.second;

                if (distances[neighborIndex] > distances[elem] + edgeWeight) {
                    distances[neighborIndex] = distances[elem] + edgeWeight;
                }
            }
        }
    }

    return distances;
}

std::vector<int> WeightedGraph::DijkstrasAlgorithm(int src) const {
    int n = adjList.size();
    std::vector<int> dist(n, std::numeric_limits<int>::max());
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<double, int>>> pq;

    dist[src] = 0;
    pq.push({0.0, src});

    while (!pq.empty()) {
        int u = pq.top().second;
        double d = pq.top().first;
        pq.pop();

        if (d > dist[u]) continue;

        for (const auto& neighbor : adjList[u]) {
            int v = neighbor.first;
            double weight = neighbor.second;

            if (dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pq.push({dist[v], v});
            }
        }
    }

    return dist;
}

int WeightedGraph::PrimsAlgorithm(int src) const
{
    if (src < 0 || src >= adjList.size()) {
        std::cerr << "Invalid source to start: " << std::endl;
        exit(1);
    }

    int n = adjList.size();
    std::vector<bool> inMST (n, false);
    std::vector<int> dist (n, std::numeric_limits<int>::max());

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    dist[src] = 0;
    pq.push({0, src});
    int MSTWeight = 0;

    while (!pq.empty()) {
        auto [weight, u] = pq.top();
        pq.pop();
        if (inMST[u]) continue;

        inMST[u] = true;
        MSTWeight += weight;

        for (const auto& [v, nextWeight] : adjList[u]) {
            if (!inMST[v] && nextWeight < dist[v]) {
                dist[v] = nextWeight;
                pq.push({dist[v], v});
            }
        }
    }

    return MSTWeight;
}


int WeightedGraph::KruskalsAlgorithm() const
{
    std::vector<std::pair<int, std::pair<int, int>>> edges;
    int n = adjList.size();
    for (int u = 0; u < n; ++u) {
        for (const auto& [v, weight] : adjList[u]) {
            edges.push_back({weight, {u, v}});
        }
    }

    //sorting edges by weight
    std::sort(edges.begin(), edges.end());

    DisjointSet uf(n);
    int totalWeight = 0;

    for (const auto& [weight, edge] : edges) {
        int u = edge.first;
        int v = edge.second;

        if (uf.findRepresentive(u) != uf.findRepresentive(v)) {
            uf.unify(u, v);
            totalWeight += weight;
        }
    }

    return totalWeight; 

}