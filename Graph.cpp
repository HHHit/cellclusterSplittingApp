#include <iostream>
#include <list>
#include "Graph.h"
using namespace std;
// Graph class represents a undirected graph
// using adjacency list representation


// Method to print connected components in an
// undirected graph
void Graph::connectedComponents()
{
    // Mark all the vertices as not visited
    bool* visited = new bool[V];
    int curLabel = 0;
    for (int v = 0; v < V; v++)
        visited[v] = false;

    for (int v = 0; v < V; v++) {
        if (visited[v] == false) {
            // print all reachable vertices
            // from v
            DFSUtil(v, visited, labelsOutput, curLabel);
//             for (int k = 0; k < V; k++) {
//                 cout << labelsOutput[k] << ",";
//             }
//             cout << "\n";
            ++curLabel;
//             cout << curLabel << "\n";

        }
    }
    delete[] visited;

}

void Graph::DFSUtil(int v, bool visited[], int labels[], int curLabel)
{
    // Mark the current node as visited and print it
    visited[v] = true;
    labels[v] = curLabel;
    // Recur for all the vertices
    // adjacent to this vertex
    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i) {
        if (!visited[*i]) {
            
            DFSUtil(*i, visited, labels, curLabel);


        }
    

    }

}

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
    labelsOutput = new int[V];
}

Graph::~Graph() { delete[] adj; }

// method to add an undirected edge
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}

// // Driver code
// int main()
// {
//     // Create a graph given in the above diagram
//     Graph g(5); // 5 vertices numbered from 0 to 4
//     g.addEdge(1, 0);
//     g.addEdge(2, 3);
//     g.addEdge(3, 1);
// 
//     cout << "Following are connected components \n";
//     g.connectedComponents();
//     for (int k = 0; k < 5; k++) {
//         cout << g.labelsOutput[k] << ",";
//     }
//     cout << "\n";
//     return 0;
// }