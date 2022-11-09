#include <iostream>
#include <list>
using namespace std;
class Graph {
    int V; // No. of vertice
    // Pointer to an array containing adjacency lists
    list<int>* adj;
    // A function used by DFS
    void DFSUtil(int v, bool visited[], int labels[], int curLabel);

public:
    int* labelsOutput; // pointer to an array of labels
    Graph(int V); // Constructor
    ~Graph();
    void addEdge(int v, int w);
    void connectedComponents();
};