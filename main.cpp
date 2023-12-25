#include "prototypes.hpp"

int main()
{
    Graph<int> obj(5);
    obj.add_vertex();
    obj.add_vertex();
    obj.add_vertex();
    obj.add_edge(0,1);
    obj.add_edge(0,2);
    obj.add_edge(0,1);
    obj.add_edge(1,7);
    obj.add_edge(3,6);
    obj.add_edge(7,5);
    //obj.remove_vertex(1);
    //obj.remove_edge(0,2);
    //obj.remove_edge(6,5);
    //obj.get_neighbours(0);
    //std::cout << obj.has_vertex(8);
    std::cout << obj.edge_count();
    std::cout << std::endl;
    obj.print();
}