/*#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono> // for timing
#include "Classes/Graph.h"

using namespace std;
using namespace std::chrono;

int n;

vector<int> best_path;
double best_cost = numeric_limits<double>::infinity();
Graph<int> graph;

Graph<int> readAndParseStadium() {
    ifstream file("../Data/Toy-Graphs/stadiums.csv");
    string line;

    getline(file, line); // Skip the header line

    while (getline(file, line)) {
        stringstream ss(line);
        string src_str, dst_str, dist_str;

        getline(ss, src_str, ',');
        int src = stoi(src_str);

        getline(ss, dst_str, ',');
        int dst = stoi(dst_str);

        getline(ss, dist_str);
        double dist = stod(dist_str);

        graph.addVertex(src);
        graph.addVertex(dst);

        graph.addEdge(src, dst, dist);
    }

    file.close();
    return graph;
}

bool recursive_backtracking(vector<int>& path, int pos, double cost) {
    if (pos == graph.getNumVertex()) {
        for (auto e : graph.findVertex(path[pos - 1])->getAdj()) {
            if (e->getDest() == graph.findVertex(0)) {
                double total_cost = cost + e->getWeight();
                if (total_cost < best_cost) {
                    best_cost = total_cost;
                    best_path = path;
                    best_path.push_back(best_path[0]); // Add the initial vertex at the end to close the cycle
                    return true;
                }
            }
        }
        return false;
    }

    for (int i = 0; i < graph.getNumVertex(); ++i) {
        if (!graph.findVertex(i)->isVisited()) {
            auto v = graph.findVertex(pos - 1);
            for (auto e : v->getAdj()) {
                if (e->getDest() == graph.findVertex(i)) {
                    path[pos] = i;
                    graph.findVertex(i)->setVisited(true);
                    if (recursive_backtracking(path, pos + 1, cost + e->getWeight())) {
                        return true;
                    }
                    graph.findVertex(i)->setVisited(false); // Backtrack
                }
            }
        }
    }
    return false;
}

vector<int> tsp(Graph<int>& graph) {
    int n = graph.getNumVertex();
    vector<int> path(n);
    path[0] = 0; // Start from the initial vertex
    graph.findVertex(0)->setVisited(true); // Mark the initial vertex as visited
    recursive_backtracking(path, 1, 0);
    return best_path;
}

int main() {
    auto start = high_resolution_clock::now(); // Start timing
    Graph<int> graph = readAndParseStadium();

    for (auto v : graph.getVertexSet())
        v->setVisited(false);

    vector<int> solution = tsp(graph);

    auto stop = high_resolution_clock::now(); // Stop timing
    auto duration = duration_cast<milliseconds>(stop - start); // Calculate duration

    cout << "Melhor caminho encontrado: ";
    for (int city : solution) {
        cout << city << " ";
    }
    cout << endl;

    cout << "Custo do melhor caminho: " << best_cost << endl;

    cout << "Tempo de execução: " << duration.count() << " milliseconds" << endl;

    return 0;
}
*/




#include "iostream"
#include "App.h"



int main() {




        App::run();


    return 0;
}









