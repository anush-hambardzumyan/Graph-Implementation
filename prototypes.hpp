#ifndef PROTOTYPES_HPP
#define PROTOTYPES_HPP


#include <iostream>
#include <vector>


template <typename T>
class Graph
{
    private:
    std::vector<std::vector<int>> vec;

    public:
    Graph(int x): vec(x , std::vector<T>()) {}
    void add_vertex();                                       //done
    void add_edge(T vertex1 , T vertex2);                    //done
    void print();                                            //done
    void remove_vertex(T vertex);                            //done
    void remove_edge(T vertex1 , T vertex2);                 //done
    void get_neighbours(T vertex);                           //done
    bool has_vertex(T vertex);                               //done   
    int vertex_count();                                      //done
    int edge_count();                                        //done
};

#include "implementations.hpp"
#endif      //PROTOTYPES_HPP