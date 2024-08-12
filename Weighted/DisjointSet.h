#ifndef DISJOINT_SET
#define DISJOINT_SET
#include <vector>

class DisjointSet { //Kruskals Algorithm, Union Find
    private:    
        std::vector<int> rank, parent;
    public:
        DisjointSet(int n) : parent(n), rank(n, 0) 
        {
            for (int i = 0; i < n; ++i) {
                parent[i] = i;
            }
        }

        int findRepresentive (int u) 
        {
           if (u != parent[u])  
                return parent[u] = findRepresentive(parent[u]);
            return u;    
        }

        void unify (int first, int second) 
        {
            first = findRepresentive(first);
            second = findRepresentive(second);
            if (first == second) return;
            else if (first >= second) {
                parent[second] = first;
                rank[first] += rank[second];
            } else {
                parent[first] = second;
                rank[second] += rank[first];
            }
        }   
};

#endif //DISJOINT_SET