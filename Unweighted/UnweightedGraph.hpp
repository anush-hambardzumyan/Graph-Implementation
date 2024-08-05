#include "UnweightedGraph.h"

void UnweightedGraph::addVertex()
{
    adjList.push_back(std::vector<int>());
    adjMatrix.push_back(std::vector<bool>(adjList.size(), 0));
    for (int i = 0; i < adjMatrix.size() - 1; ++i) {
        adjMatrix[i].push_back(0);
    }
}

void UnweightedGraph::addEdge(int vertex1, int vertex2, bool undirected = true)
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

void UnweightedGraph::printList() const
{
    for(int i = 0 ; i < adjList.size() ; i++)
    {
        std::cout << i << " -> ";
        for(int j = 0 ; j < adjList[i].size() ; j++)
        {
            std::cout << adjList[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void UnweightedGraph::printMatrix() const
{
    const int width = 4;

    std::cout << std::setw(width) << " ";
    for (int i = 0; i < adjMatrix.size(); ++i) {
        std::cout << std::setw(width) << i;
    }
    std::cout << std::endl;

    for (int i = 0; i < adjMatrix.size(); ++i) {
        std::cout << std::setw(width) << i; 
        for (int j = 0; j < adjMatrix[i].size(); ++j) {
            std::cout << std::setw(width) << adjMatrix[i][j];
        }
        std::cout << std::endl;
    }
}

void UnweightedGraph::removeVertex(int vertex)
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
            if (adjList[i][j] == vertex) 
            {
                adjList[i].erase(adjList[i].begin() + j);
            }
        }   
    }   
}

void UnweightedGraph::removeEdge(int vertex1, int vertex2, bool undirected = true)
{
    if (vertex1 < 0 || vertex1 >= adjList.size() || vertex2 < 0 || vertex2 >= adjList.size())
    {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
        return;
    }

    for (int i = 0; i < adjList[vertex1].size(); ++i)
    {
        if (adjList[vertex1][i] == vertex2)
        {
            adjList[vertex1].erase(adjList[vertex1].begin() + i);
            return;
        }
    }
    adjMatrix[vertex1][vertex2] = false;

    if(undirected) {
        for (int i = 0; i < adjList[vertex2].size(); ++i)
        {
            if (adjList[vertex2][i] == vertex1)
            {
                adjList[vertex2].erase(adjList[vertex2].begin() + i);
                return;
            }
        }
        adjMatrix[vertex2][vertex1] = false;
    }
    
}

void UnweightedGraph::getNeighbours(int vertex) const
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
        std::cout << adjList[vertex][i] << " ";
    }

    std::cout << std::endl;
}

bool UnweightedGraph::hasVertex(int vertex) const
{
    return !(vertex < 0 || vertex >= adjList.size());
}

int UnweightedGraph::vertexCount() const
{
    return adjList.size();
}

int UnweightedGraph::edgeCount() const 
{
    if(!adjList.size())
    {
        return 0;
    }

    int count = 0;
    for(int i = 0; i < adjList.size(); ++i)
    {
        for(int j = 0; j < adjList[i].size(); ++j)
        {
            count++;
        }
    }
    return count;
}

void UnweightedGraph::dfs(int src,bool recursive = true) const
{
    if (src < 0 || src >= vertexCount()) {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
    }

    if (recursive) {
        std::vector<bool> visited(vertexCount(), false);
        dfsRecursive(src, visited);
        std::cout << std::endl;
    }
    
    else dfsIterative(src);
}

void UnweightedGraph::dfsRecursive(int src, std::vector<bool>& visited, bool print) const
{
    if (src < 0 || src >= vertexCount()) {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
    }

    visited[src] = true;
    if (print) { 
        std::cout << src << " ";
    }
        

    for (auto elem : adjList[src]) {
        if (!visited[elem]) {
            dfsRecursive(elem, visited, print);
        }
    }
}

void UnweightedGraph::dfsIterative(int src) const
{
    if (src < 0 || src >= vertexCount()) {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
    }

    std::vector<bool> visited(vertexCount(), false);
    std::stack<int> stack;
    stack.push(src);

    while(!stack.empty()) {
        int val = stack.top();
        stack.pop();
        if (!visited[val]) {
            std::cout << val << " ";
            visited[val] = true;
        }
        
        for (auto it = adjList[val].rbegin(); it != adjList[val].rend(); ++it) {
            if (!visited[*it]) 
                stack.push(*it);
            
        }
    }
    std::cout << std::endl;
}

void UnweightedGraph::bfs(int src) const
{
    if (src < 0 || src >= vertexCount()) {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
    }

    std::vector<bool> visited(vertexCount(), false);
    std::queue<int> queue;
    queue.push(src);
    visited[src] = true;

    while (!queue.empty()) {
        int val = queue.front();
        queue.pop();

        std::cout << val << " ";
        for (auto elem : adjList[val]) {
            if (!visited[elem]) {
                visited[elem] = true;
                queue.push(elem);
            }
        }
    }
    std::cout << std::endl;
}

void UnweightedGraph::dfsExtraCase(int src) const
{
    if (src < 0 || src >= vertexCount()) {
        std::cout << "Invalid operation: " << std::endl;
        exit(0);
    }

    countConnectedComponents(true);
}

int UnweightedGraph::countConnectedComponents(bool printDFS) const
{
    std::vector<bool> visited(vertexCount(), false);  
    int count = 0;

    for (int v = 0; v < vertexCount(); ++v) {
        if (!visited[v]) {
            dfsRecursive(v, visited, printDFS);
            ++count;
        }
    }

    return count;
}

void UnweightedGraph::transpose()
{
    std::vector<std::vector<int>> transposedList(adjList.size());
    for (int i = 0; i < adjList.size(); ++i) {
        for (int j = 0; j < adjList[i].size(); ++j) {
            transposedList[adjList[i][j]].push_back(i);
        }
    }

    adjList = transposedList;

    std::vector<std::vector<bool>> transposedMatrix (adjMatrix.size(), std::vector<bool>(adjMatrix.size(), false));
    for (int i = 0; i < adjMatrix.size(); ++i) {
        for (int j = 0; j < adjMatrix[i].size(); ++j) {
            transposedMatrix[j][i] = adjMatrix[i][j];
        }
    }

    adjMatrix = transposedMatrix;
}

std::vector<std::vector<int>> UnweightedGraph::getTransposedAdjList() const
{
    std::vector<std::vector<int>> transposedList(adjList.size());
    for (int i = 0; i < adjList.size(); ++i) {
        for (int neighbor : adjList[i]) {
            transposedList[neighbor].push_back(i);
        }
    }
    return transposedList;
}

std::vector<std::vector<bool>> UnweightedGraph::getTransposedAdjMatrix() const
{
    std::vector<std::vector<bool>> transposedMatrix (adjMatrix.size(), std::vector<bool>(adjMatrix.size(), false));
    for (int i = 0; i < adjMatrix.size(); ++i) {
        for (int j = 0; j < adjMatrix[i].size(); ++j) {
            transposedMatrix[j][i] = adjMatrix[i][j];
        }
    }

    return transposedMatrix;
}

std::vector<int> UnweightedGraph::getShortestPath(int src, int dest) const
{
    int n = adjList.size();
    std::vector<bool> visited(n, false);
    std::vector<int> parent(n, -1);
    std::queue<int> queue;

    visited[src] = true;
    queue.push(src);

    while(!queue.empty()) {
        int current = queue.front();
        queue.pop();

        if (current == dest) {
            std::vector<int> path;
            for (int v = dest; v != -1; v = parent[v]) {
                path.push_back(v);
            }
            std::reverse(path.begin(), path.end());
            return path;
        }

        for (int i = 0; i < adjList[current].size(); ++i) {
            if (!visited[adjList[current][i]]) {
                visited[adjList[current][i]] = true;
                parent[adjList[current][i]] = current;
                queue.push(adjList[current][i]);
            }
        }
    }
    return {};
}

std::vector<std::vector<int>> UnweightedGraph::getAllPaths(int src, int dest) const
{
    if (src < 0 || src > adjList.size() - 1 || dest < 0 || dest > adjList.size() - 1) {
        std::cerr << "Error: Source or destination index is out of bounds." << std::endl;
        return {};
    }

    int n = adjList.size();
    std::vector<bool> visited(n, false);
    std::vector<std::vector<int>> allPaths;
    std::vector<int> currPath;
    allPathHelper(src, dest, currPath, visited, allPaths);
    return allPaths;

}

void UnweightedGraph::allPathHelper(int src, int dest, std::vector<int>& currPath, std::vector<bool>& visited, std::vector<std::vector<int>>& allPaths) const
{
    visited[src] = true;
    currPath.push_back(src);
    
    if (src == dest) {
        allPaths.push_back(currPath);
    } else {
        for (auto i : adjList[src]) {
            if (!visited[i]) {
                allPathHelper(i, dest, currPath, visited, allPaths);  
            }
        }
    }
    
    currPath.pop_back();
    visited[src] = false;
}

bool UnweightedGraph::isCycledUndirected() const
{
    std::vector<bool> visited(adjList.size(), false);
    int parent = -1;
    for (int i = 0; i < adjList.size(); ++i) {
        if (!visited[i]) {
            if(isCycledUndirectedHelper(i, parent, visited)) {
                return true;
            }
        }
    }
    return false;
}

bool UnweightedGraph::isCycledUndirectedHelper(int src, int parent,std::vector<bool>& visited) const
{
    visited[src] = true;
    for (int neighbor : adjList[src]) {
        if (!visited[neighbor]) {
            if (isCycledUndirectedHelper(neighbor, src, visited)) {
                return true;
            }
        } else if (parent != neighbor) {
            return true;
        }
    }
    return false;
} 

bool UnweightedGraph::isCycledDirected() const
{
    std::vector<bool> visited(adjList.size(), false);
    std::vector<bool> onStack(adjList.size(), false);

    for (int i = 0; i < adjList.size(); ++i) {
        if (!visited[i]) {
            if (isCycledDirectedHelper(i, onStack, visited)) {
                return true;
            }
        }
    }
    return false;
}

bool UnweightedGraph::isCycledDirectedHelper(int src, std::vector<bool>& onStack, std::vector<bool>& visited) const
{
    visited[src] = true;
    onStack[src] = true;
    for (int neighbor : adjList[src]) {
        if (!visited[neighbor]) {
            if (isCycledDirectedHelper(neighbor, onStack, visited)) {
                return true;
            }
        } else if (onStack[neighbor]) {
            return true;
        }
    }

    onStack[src] = false;
    return false;
}

std::vector<int> UnweightedGraph::getNodesOfKthLevel(int k, bool dfs) const
{   
    std::vector<bool> visited(adjList.size(), false);
    std::vector<int> res;
    getNodesOfLevelDFSHelper(0, k, res, visited);
    return res;
}

void UnweightedGraph::getNodesOfLevelDFSHelper(int src, int level, std::vector<int>& res, std::vector<bool>& visited) const
{
    visited[src] = true;
    if (level == 1) res.push_back(src);
    else {
        for (int i = 0; i < adjList[src].size(); ++i) {
            if (!visited[adjList[src][i]]) {
                getNodesOfLevelDFSHelper(adjList[src][i], level - 1, res, visited);
            }
        }
    }
}

void UnweightedGraph::topSortHelper(int src, std::vector<bool>& visited, std::list<int>& res) const
{
    visited[src] = true;

    for (int i : adjList[src]) {
        if (!visited[i])
            topSortHelper(i, visited, res);
    }

    res.push_front(src);
}

std::list<int> UnweightedGraph::topologicalSort() const
{
    if (isCycledDirected()) {
        std::cout << "Graph is cycled: " << std::endl;
        exit(1);
    }

    std::vector<bool> visited(adjList.size(), false);
    std::list<int> result;

    for (int i = 0; i < adjList.size(); i++) {
        if (!visited[i])
            topSortHelper(i, visited, result);
    }

    return result;
}

std::list<int> UnweightedGraph::KahnsAlgorithm() const
{
    std::vector<int> inDegree(adjList.size(), 0);
    for (int i = 0; i < adjList.size(); ++i) {
        for (int elem : adjList[i]) {
            ++inDegree[elem];
        }
    }

    std::queue<int> queue;

    for (int i = 0; i < adjList.size(); ++i) {
        if (inDegree[i] == 0) {
            queue.push(i);
        }
    }
    std::list<int> result;
    while (!queue.empty()) {
        int elem = queue.front();
        result.push_front(elem);
        queue.pop();

        for (int i : adjList[elem]) {
            --inDegree[i];
            if (inDegree[i] == 0) {
                queue.push(i);
            }
        }
    } 
    if (result.size() != adjList.size()) {
            return {}; //cycled graph
        }
    return result;
}

void UnweightedGraph::dfsInVector(int src, std::vector<bool>& visited, std::vector<int>& result) const
{
    visited[src] = true;
    result.push_back(src);
    for (auto elem : adjList[src]) {
        if (!visited[elem]) {
            dfsInVector(elem, visited, result);
        }
    }
}

void UnweightedGraph::fillOrder(int src, std::vector<bool>& visited, std::stack<int>& stack) const
{
    visited[src] = true;
    for (int i : adjList[src])
        if (!visited[i])
            fillOrder(i, visited, stack);
    stack.push(src);
}

std::vector<std::vector<int>> UnweightedGraph::KosarajusAlgorithm() const
{
    std::stack<int> Stack;
    std::vector<bool> visited(adjList.size(), false);
    for (int i = 0; i < visited.size(); ++i) {
        if (!visited[i]) {
            fillOrder(i, visited, Stack);
        }
    }

    UnweightedGraph transposed(this -> getTransposedAdjList(), this -> getTransposedAdjMatrix());
    visited = std::vector<bool> (adjList.size(), false);

    std::vector<std::vector<int>> SCCs;    
    while (!Stack.empty()) {
        int v = Stack.top();
        Stack.pop(); 
        if (!visited[v]) {
            std::vector<int> component;
            transposed.dfsInVector(v, visited, component);
            SCCs.push_back(component);
        }
    }
    return SCCs;
}

void UnweightedGraph::TarjanHelper(int src, std::vector<int>& ids, std::vector<int>& lowlink, std::stack<int>& st, std::vector<bool>& onStack, std::vector<std::vector<int>>& SCCs) const
{
    static int visitingTime = 0;
    ids[src] = lowlink[src] = ++visitingTime;
    st.push(src);
    onStack[src] = true;

    for (int v : adjList[src]) {
        if (ids[v] == -1) {
            TarjanHelper(v, ids, lowlink, st, onStack, SCCs);
            lowlink[src] = std::min(lowlink[src], lowlink[v]);
        } else if (onStack[v]) {
            lowlink[src] = std::min(lowlink[src], ids[v]);
        }
    }

    if (lowlink[src] == ids[src]) {
        std::vector<int> currSCC;
        while (st.top() != src) {
            int elem = st.top();
            currSCC.push_back(elem);
            st.pop();
            onStack[elem] = false;
        }

        currSCC.push_back(st.top());
        onStack[st.top()] = false;
        st.pop();
        SCCs.push_back(currSCC);
    }
}

std::vector<std::vector<int>> UnweightedGraph::TarjansAlgorithm() const
{
    std::vector<std::vector<int>> SCCs;
    std::vector<bool> onStack(adjList.size(), false);
    std::vector<int> ids(adjList.size(), -1);
    std::vector<int> lowlink(adjList.size(), -1);
    std::stack<int> st;

    for (int i = 0; i < ids.size(); ++i) {
        if (ids[i] == -1) {
            TarjanHelper(i, ids, lowlink, st, onStack, SCCs);
        }
    }

    return SCCs;
}