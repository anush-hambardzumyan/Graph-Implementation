#include "prototypes.hpp"

template <typename T>
void Graph<T>::add_vertex()
{
    vec.push_back(std::vector<int>());
}

template <typename T>
void Graph<T>::add_edge(T vertex1 , T vertex2)
{
    if(vertex1 < 0 || vertex1 >= vec.size() || vertex2 < 0 || vertex2 >= vec.size())
    {
        return;
    }

    auto it = std::find(vec[vertex1].begin() , vec[vertex1].end() , vertex2 );
    if(it == vec[vertex1].end())
    {
        vec[vertex1].push_back(vertex2);
    }
}

template<typename T>
void Graph<T>::print()
{
    for(int i = 0 ; i < vec.size() ; i++)
    {
        std::cout << i << " -> ";
        for(int j = 0 ; j < vec[i].size() ; j++)
        {
            std::cout << vec[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
void Graph<T>::remove_vertex(T vertex)
{
    if(vertex < 0 || vertex > vec.size() - 1)
    {
        std::cout << "Invalid operation: " << std::endl;
        return;
    }

    vec.erase(vec.begin() + vertex);

    for (int i = 0; i < vec.size(); ++i) 
    {
        for (int j = 0; j < vec[i].size(); ++j) 
        {
            if (vec[i][j] == vertex) 
            {
                vec[i].erase(vec[i].begin() + j);
            }
        }   
    }   
}

template<typename T>
void Graph<T>::remove_edge(T vertex1, T vertex2)
{
    if (vertex1 < 0 || vertex1 >= vec.size() || vertex2 < 0 || vertex2 >= vec.size())
    {
        std::cout << "Invalid operation: " << std::endl;
        return;
    }

    for (int i = 0; i < vec[vertex1].size(); ++i)
    {
        if (vec[vertex1][i] == vertex2)
        {
            vec[vertex1].erase(vec[vertex1].begin() + i);
            return;
        }
    }
    std::cout << "Edge wasn't found: " << std::endl;
}


template<typename T>
void Graph<T>::get_neighbours(T vertex)
{
    if(vertex < 0 || vertex >= vec.size())
    {
        std::cout << "Invalid operation: " << std::endl;
        return; 
    }

    std::cout << "Neighbours of vertex are: "; 
    for(int i = 0; i < vec[vertex].size(); ++i)
    {
        std::cout << vec[vertex][i] << " ";
    }

    std::cout << std::endl;
}

template<typename T>
bool Graph<T>::has_vertex(T vertex)
{
    if(vertex < 0 || vertex >= vec.size())
    {
        return false;
    }
    return true;
}

template<typename T>
int Graph<T>::vertex_count()
{
    return vec.size();
}


template<typename T>
int Graph<T>::edge_count()
{
    if(!vec.size())
    {
        return 0;
    }

    int count = 0;
    for(int i = 0; i < vec.size(); ++i)
    {
        for(int j = 0; j < vec[i].size(); ++j)
        {
            count++;
        }
    }
    return count;
}









