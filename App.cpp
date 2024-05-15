//
// Created by David Ferreira on 11/05/2024.
//

#include "App.h"


#include <iostream>
#include <cmath>
#include <utility>
#include <set>
#include <chrono>
#include "reader.h"

using namespace std;


void display4_1menu();
void display4_2menu();
void display4_1menuSmallGraph();
void display4_1menuSmallGraphStadium();
void display4_1menuSmallGraphShipping();
void display4_1menuSmallGraphTourism();
void display4_1menuMediumGraph();
void display4_1menuMedium(const std::string &filename);
void display4_2menuSmallGraph();
void display4_2menuMediumGraph();
void display4_2menuSmallGraphStadium();
void display4_2menuSmallGraphShipping();
void display4_2menuSmallGraphTourism();
void display4_2menuMedium(const std::string &filename);

void display_OHmenu();
void display_NNmenu();
void display_NNmenuSmallGraph();
void display_NNmenuMediumGraph();
void getValue_NNmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc);
void getValue_NNmenuSmallGraphStadium();
void getValue_NNmenuSmallGraphShipping();
void getValue_NNmenuSmallGraphTourism();
void getValue_NNmenuMediumGraph(const std::string &filename);

void display_SANmenu();
void display_SANmenuSmallGraph();
void display_SANmenuMediumGraph();
void getValue_SANmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc);
void getValue_SANmenuSmallGraphStadium();
void getValue_SANmenuSmallGraphShipping();
void getValue_SANmenuSmallGraphTourism();
void getValue_SANmenuMediumGraph(const std::string &filename);

void display_LINmenu();
void display_LINmenuSmallGraph();
void display_LINmenuMediumGraph();
void getValue_LINmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc);
void getValue_LINmenuSmallGraphStadium();
void getValue_LINmenuSmallGraphShipping();
void getValue_LINmenuSmallGraphTourism();
void getValue_LINmenuMediumGraph(const std::string &filename);

void display_HeldKarp_menu();
void display_HeldKarp_menuSmallGraph();
void display_HeldKarp_menuMediumGraph();
void getValue_HeldKarp_menuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc);
void getValue_HeldKarp_menuSmallGraphStadium();
void getValue_HeldKarp_menuSmallGraphShipping();
void getValue_HeldKarp_menuSmallGraphTourism();
void getValue_HeldKarp_menuMediumGraph(const std::string &filename);

void display_CLUSTERmenu(int clusterOption);
void display_CLUSTERmenuSmallGraph(int clusterOption);
void display_CLUSTERmenuMediumGraph(int clusterOption);
void getValue_CLUSTERmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc, int clusterOption);
void getValue_CLUSTERmenuSmallGraphStadium(int clusterOption);
void getValue_CLUSTERmenuSmallGraphShipping(int clusterOption);
void getValue_CLUSTERmenuSmallGraphTourism(int clusterOption);
void getValue_CLUSTERmenuMediumGraph(const std::string &filename, int clusterOption);


void display_RWmenu();
void getValue_RWsmallGraph();
void getValue_RWmediumGraph();
void getValue_RWlargeGraph();

std::vector<int> preOrderTraversal(Vertex<int>* root, const std::vector<Edge<int>*>& mst);
void tsp(int currentNode, std::vector<int>& path, double currentCost, int level, Graph<int>& graph, std::vector<int>& bestPath, double& bestCost);
void heldKarp(const Graph<int>& graph, std::vector<int>& bestPath, double& bestCost);
std::vector<int> generateInitialTour(Graph<int> graph);
std::vector<int> generateNeighbor(const std::vector<int>& tour);
bool shouldAccept(double delta, double temperature);
std::vector<int> simulatedAnnealing(Graph<int>& graph, double initialTemperature, double coolingRate, int iterations);


int mainMenu(){
    cout << "Loading...";

    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Main Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option that suits your needs:\n";
        cout << "1. Backtracking Algorithm\n";
        cout << "2. Triangular Approximation Heuristic \n";
        cout <<"3. Other Heuristics \n";
        cout << "4 Tsp in the Real World\n";
        cout << "e. Exit\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display4_1menu();
                break;
            case '2':
                display4_2menu();
                break;
            case '3':
                display_OHmenu();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
    return 0;
}

void display4_1menu() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Backtracing Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Small Graph\n";
        cout << "2. Medium Graph \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display4_1menuSmallGraph();
                break;
            case '2':
                display4_1menuMediumGraph();
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void display4_2menu() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Triangular Approximation Heuristic       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Small Graph\n";
        cout << "2. Medium Graph \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display4_2menuSmallGraph();
                break;
            case '2':
                display4_2menuMediumGraph();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void display4_1menuSmallGraph(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Backtracing SmallGraph Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the graph you want:\n";
        cout << "1. Stadium Graph\n";
        cout << "2. Shipping Graph \n";
        cout <<"3. Tourism Graph \n";
        cout << "e. Back to the backtracing Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display4_1menuSmallGraphStadium();
                break;
            case '2':
                display4_1menuSmallGraphShipping();
                break;
            case '3':
                display4_1menuSmallGraphTourism();
            case 'e':
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void display4_1menuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc) {
    Reader reader;
    Graph<int> graph;
    graph = readAndParseFunc(reader);
    std::vector<int> bestPath;
    double bestCost = std::numeric_limits<double>::infinity();
    std::vector<int> path(1, 0);
    graph.findVertex(0)->setVisited(true);

    auto startTime = std::chrono::high_resolution_clock::now();
    tsp(0, path, 0.0, 1, graph, bestPath, bestCost);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTime = endTime - startTime;

    std::cout << "Execution Time: " << execTime.count() << " seconds\n";
    std::cout << "Best Path Cost: " << bestCost << "\n";
    std::cout << "Best Path: ";
    for (size_t i = 0; i < bestPath.size(); ++i) {
        std::cout << bestPath[i];
        if (i < bestPath.size() - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << "\n";
}

void display4_1menuSmallGraphStadium() {
    auto readAndParseStadium = [](Reader& reader) { return reader.readAndParseStadium(); };
    display4_1menuSmallGraph(readAndParseStadium);
}

void display4_1menuSmallGraphShipping() {
    auto readAndParseShipping = [](Reader& reader) { return reader.readAndParseShipping(); };
    display4_1menuSmallGraph(readAndParseShipping);
}

void display4_1menuSmallGraphTourism() {
    auto readAndParseTourism = [](Reader& reader) { return reader.readAndParseTourism(); };
    display4_1menuSmallGraph(readAndParseTourism);
}

void display4_1menuMediumGraph() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Triangular Approximation Medium Graph Menu       \n";
        cout << "-----------------------------\n";
        cout << "How many edges does the graph have:\n";
        cout << "1. 25\n";
        cout << "2. 50\n";
        cout << "3. 75\n";
        cout << "4. 100\n";
        cout << "5. 200\n";
        cout << "6. 300\n";
        cout << "7. 400\n";
        cout << "8. 500\n";
        cout << "9. 600\n";
        cout << "10. 700\n";
        cout << "11. 800\n";
        cout << "12. 900\n";
        cout << "e. Back to the backtracing Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice == "e") {
            cout << "Exiting menu system...\n";
            exitMenu = true;
        } else {
            // Convert choice to an integer for comparison
            int choiceNum = stoi(choice);
            switch (choiceNum) {
                case 1:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
                    break;
                case 2:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_50.csv");
                    break;
                case 3:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_75.csv");
                    break;
                case 4:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_100.csv");
                    break;
                case 5:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_200.csv");
                    break;
                case 6:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_300.csv");
                    break;
                case 7:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_400.csv");
                    break;
                case 8:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_500.csv");
                    break;
                case 9:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_600.csv");
                    break;
                case 10:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_700.csv");
                    break;
                case 11:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_800.csv");
                    break;
                case 12:
                    display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_900.csv");
                    break;
                default:
                    cout << "Invalid input. Please choose a valid option.\n";
            }
        }
    }
}

void display4_1menuMedium(const std::string &filename){
    Reader reader;
    Graph<int> graph;
    graph = reader.readAndParseExtra_Fully_Connected_Graphs(filename);
    std::vector<int> bestPath;
    double bestCost = std::numeric_limits<double>::infinity();
    std::vector<int> path(1, 0);
    graph.findVertex(0)->setVisited(true);

    auto startTime = std::chrono::high_resolution_clock::now();
    tsp(0, path, 0.0, 1, graph, bestPath, bestCost);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTime = endTime - startTime;

    std::cout << "Execution Time: " << execTime.count() << " seconds\n";
    std::cout << "Best Path Cost: " << bestCost << "\n";
    std::cout << "Best Path: ";
    for (int node : bestPath) {
        std::cout << node << (node == 0 ? "\n" : " -> ");
    }

}

void display4_2menuSmallGraph(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to  Triangular Approximation Heuristic SmallGraph Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the graph you want:\n";
        cout << "1. Stadium Graph\n";
        cout << "2. Shipping Graph \n";
        cout <<"3. Tourism Graph \n";
        cout << "e. Back to the backtracing Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display4_2menuSmallGraphStadium();
                break;
            case '2':
                display4_2menuSmallGraphShipping();
                break;
            case '3':
                display4_2menuSmallGraphTourism();
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void display4_2menuSmall(const std::function<Graph<int>(Reader&)>& readAndParseFunc) {

    Reader reader;
    Graph<int> graph = readAndParseFunc(reader);
    Vertex<int>* startVertexPtr = graph.findVertex(0);

    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    vector<Edge<int>*> mst = graph.primMST(startVertexPtr->getInfo());

    vector<int> preorderList = preOrderTraversal(startVertexPtr, mst);

    vector<int> tspTour;
    set<int> visited;
    for (int vertex : preorderList) {
        if (visited.insert(vertex).second) {
            tspTour.push_back(vertex);
        }
    }
    tspTour.push_back(tspTour.front());

    double totalDistance = 0;
    int previous = tspTour.front();
    for (size_t i = 1; i < tspTour.size(); i++) {
        Edge<int> *edge = graph.findEdge(previous, tspTour[i]);
        if (edge)
            totalDistance += edge->getWeight();
        previous = tspTour[i];
    }

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> execTime = end - start; // Calculate execution time

    cout << "TSP Tour: ";
    for (size_t i = 0; i < tspTour.size(); ++i) {
        cout << tspTour[i] << (i < tspTour.size() - 1 ? " -> " : "");
    }
    cout << endl;
    cout << "Total Approximation Distance: " << totalDistance << "\n";
    cout << "Execution Time: " << execTime.count() << " seconds\n"; // Print execution time
}

void display4_2menuSmallGraphStadium() {
    auto readAndParseStadium = [](Reader& reader) { return reader.readAndParseStadium(); };
    display4_2menuSmall(readAndParseStadium);
}

void display4_2menuSmallGraphShipping() {
    auto readAndParseShipping = [](Reader& reader) { return reader.readAndParseShipping(); };
    display4_2menuSmall(readAndParseShipping);
}

void display4_2menuSmallGraphTourism() {
    auto readAndParseTourism = [](Reader& reader) { return reader.readAndParseTourism(); };
    display4_2menuSmall(readAndParseTourism);
}

void display4_2menuMediumGraph() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Triangular Approximation Medium Graph Menu       \n";
        cout << "-----------------------------\n";
        cout << "How many edges does the graph have:\n";
        cout << "1. 25\n";
        cout << "2. 50\n";
        cout << "3. 75\n";
        cout << "4. 100\n";
        cout << "5. 200\n";
        cout << "6. 300\n";
        cout << "7. 400\n";
        cout << "8. 500\n";
        cout << "9. 600\n";
        cout << "10. 700\n";
        cout << "11. 800\n";
        cout << "12. 900\n";
        cout << "e. Back to the backtracing Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice == "e") {
            cout << "Exiting menu system...\n";
            exitMenu = true;
        } else {
            // Convert choice to an integer for comparison
            int choiceNum = stoi(choice);
            switch (choiceNum) {
                case 1:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
                    break;
                case 2:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_50.csv");
                    break;
                case 3:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_75.csv");
                    break;
                case 4:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_100.csv");
                    break;
                case 5:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_200.csv");
                    break;
                case 6:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_300.csv");
                    break;
                case 7:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_400.csv");
                    break;
                case 8:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_500.csv");
                    break;
                case 9:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_600.csv");
                    break;
                case 10:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_700.csv");
                    break;
                case 11:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_800.csv");
                    break;
                case 12:
                    display4_2menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_900.csv");
                    break;
                default:
                    cout << "Invalid input. Please choose a valid option.\n";
            }
        }
    }
}


void display4_2menuMedium(const std::string &filename){
    Reader reader;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename);


    Vertex<int>* startVertexPtr = graph.findVertex(0);

    auto start = std::chrono::high_resolution_clock::now();

    vector<Edge<int>*> mst = graph.primMST(startVertexPtr->getInfo());


    vector<int> preorderList = preOrderTraversal(startVertexPtr, mst);


    vector<int> tspTour;
    set<int> visited;
    for (int vertex : preorderList) {
        if (visited.insert(vertex).second) {
            tspTour.push_back(vertex);
        }
    }
    tspTour.push_back(tspTour.front());

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTime = end - start;

    double totalDistance = 0;
    int previous = tspTour.front();
    for (size_t i = 1; i < tspTour.size(); i++) {
        Edge<int> *edge = graph.findEdge(previous, tspTour[i]);
        if (edge)
            totalDistance += edge->getWeight();
        previous = tspTour[i];

    }


    cout << "TSP Tour: ";
    for (size_t i = 0; i < tspTour.size(); ++i) {
        cout << tspTour[i] << (i < tspTour.size() - 1 ? " -> " : "");
    }
    cout << endl;
    cout << "Total Approximation Distance: " << totalDistance << "\n";
    cout << "Execution Time: " << execTime.count() << " seconds\n";

}

std::vector<int> preOrderTraversal(Vertex<int>* root, const std::vector<Edge<int>*>& mst) {
    std::vector<int> order;
    std::unordered_set<Vertex<int>*> visited;
    std::unordered_map<Vertex<int>*, std::vector<Vertex<int>*>> tree;

    for (auto edge : mst) {
        tree[edge->getOrig()].push_back(edge->getDest());
        tree[edge->getDest()].push_back(edge->getOrig());
    }

    std::function<void(Vertex<int>*)> dfs = [&](Vertex<int>* vertex) {
        if (!vertex || visited.count(vertex)) return;
        visited.insert(vertex);
        order.push_back(vertex->getInfo());

        for (auto next : tree[vertex]) {
            dfs(next);
        }
    };

    dfs(root);
    return order;
}


void tsp(int currentNode, std::vector<int>& path, double currentCost, int level, Graph<int>& graph, std::vector<int>& bestPath, double& bestCost) {
    if (level == graph.getNumVertex()) {
        Edge<int>* edgeBack = graph.findEdge(currentNode, 0);
        if (edgeBack && currentCost + edgeBack->getWeight() < bestCost) {
            bestCost = currentCost + edgeBack->getWeight();
            bestPath = path;
            bestPath.push_back(0);
        }
        return;
    }

    if (currentCost >= bestCost)
        return;

    Vertex<int>* currentVertex = graph.findVertex(currentNode);
    for (Edge<int>* edge : currentVertex->getAdj()) {
        auto nextNode = edge->getDest();
        if (!nextNode->isVisited()) {
            nextNode->setVisited(true);
            path.push_back(nextNode->getInfo());
            tsp(nextNode->getInfo(), path, currentCost + edge->getWeight(), level + 1, graph, bestPath, bestCost);
            nextNode->setVisited(false);
            path.pop_back();
        }
    }
}

void heldKarp(const Graph<int>& graph, std::vector<int>& bestPath, double& bestCost) {
    int n = graph.getNumVertex();

    // dp[mask][i] will be the minimum cost to reach node i with visited nodes as in mask.
    std::vector<std::vector<double>> dp(1 << n, std::vector<double>(n, INF));
    std::vector<std::vector<int>> parent(1 << n, std::vector<int>(n, -1));
    dp[1][0] = 0; // Starting at node 0.

    // Iterate over all subsets of nodes.
    for (int mask = 1; mask < (1 << n); ++mask) {
        for (int u = 0; u < n; ++u) {
            if (mask & (1 << u)) { // If u is in the subset represented by mask.
                for (int v = 0; v < n; ++v) {
                    if (!(mask & (1 << v))) { // If v is not in the subset.
                        Edge<int>* edge = graph.findEdge(u, v);
                        if (edge) {
                            double newCost = dp[mask][u] + edge->getWeight();
                            if (newCost < dp[mask | (1 << v)][v]) {
                                dp[mask | (1 << v)][v] = newCost;
                                parent[mask | (1 << v)][v] = u;
                            }
                        }
                    }
                }
            }
        }
    }

    // Reconstruct the minimum cost to complete the tour and the path.
    bestCost = INF;
    int lastNode = -1;
    for (int u = 1; u < n; ++u) {
        Edge<int>* edge = graph.findEdge(u, 0);
        if (edge) {
            double currentCost = dp[(1 << n) - 1][u] + edge->getWeight();
            if (currentCost < bestCost) {
                bestCost = currentCost;
                lastNode = u;
            }
        }
    }

    // Reconstruct the path using the parent table.
    if (lastNode != -1) {
        bestPath.clear();
        int mask = (1 << n) - 1;
        int currentNode = lastNode;
        while (currentNode != 0) {
            bestPath.push_back(currentNode);
            int temp = parent[mask][currentNode];
            mask ^= (1 << currentNode);
            currentNode = temp;
        }
        bestPath.push_back(0); // Add the starting node to complete the cycle.
        std::reverse(bestPath.begin(), bestPath.end()); // Reverse the path to start from the starting node.
        bestPath.push_back(0); // Add the starting node at the end to complete the cycle.
    }
}



void display_OHmenu(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Other Heuristics Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the approach you want:\n";
        cout << "1. NearestNeighbour\n";
        cout << "2. K-means Clustering NearestNeighbour\n";
        cout << "3. Simulated Annealing\n";
        cout << "4. LinKernighan\n";
        cout << "5. HeldKarp (Optimal solution but feasible only on toy graphs)";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display_NNmenu();
                break;
            case '2': {
                string clusterChoice;
                int clusterOption;
                bool validInput = false;

                // Keep asking for input until valid input is provided
                while (!validInput) {
                    cout << "Enter the number of clusters: ";
                    cin >> clusterChoice;

                    try {
                        clusterOption = stoi(clusterChoice);
                        validInput = true;
                    } catch (const invalid_argument& e) {
                        cout << "Invalid input. Please enter a valid integer.\n";
                    } catch (const out_of_range& e) {
                        cout << "Invalid input. The number is out of range.\n";
                    }
                }

                display_CLUSTERmenu(clusterOption);
                break;
            }
            case '3':
                display_SANmenu();
                break;
            case '4':
                display_LINmenu();
                break;
            case '5':
                display_HeldKarp_menu();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }

}

void display_NNmenu() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Nearest Neighbour Algorithm Heuristic       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Small Graph\n";
        cout << "2. Medium Graph \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display_NNmenuSmallGraph();
                break;
            case '2':
                display_NNmenuMediumGraph();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}
void display_NNmenuSmallGraph(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Nearest Neighbour Heuristic SmallGraph Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the graph you want:\n";
        cout << "1. Stadium Graph\n";
        cout << "2. Shipping Graph \n";
        cout <<"3. Tourism Graph \n";
        cout << "e. Back to the other Heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                getValue_NNmenuSmallGraphStadium();
                break;
            case '2':
                getValue_NNmenuSmallGraphShipping();
                break;
            case '3':
                getValue_NNmenuSmallGraphTourism();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void getValue_NNmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc) {
    Reader reader;
    Graph<int> graph;
    graph = readAndParseFunc(reader);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.nearestNeighbour(graph);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double totalDistance = 0;
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i]->getInfo() << (i < tour.size() - 1 ? " -> " : "");
        if (i > 0) {
            Edge<int>* edge = graph.findEdge(tour[i - 1]->getInfo(), tour[i]->getInfo());
            if (edge)
                totalDistance += edge->getWeight();
        }
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";

}



void getValue_NNmenuSmallGraphStadium(){
    auto readAndParseStadium = [](Reader& reader) { return reader.readAndParseStadium(); };
    getValue_NNmenuSmallGraph(readAndParseStadium);
}

void getValue_NNmenuSmallGraphShipping() {
    auto readAndParseShipping = [](Reader& reader) { return reader.readAndParseShipping(); };
    getValue_NNmenuSmallGraph(readAndParseShipping);
}

void getValue_NNmenuSmallGraphTourism() {
    auto readAndParseTourism = [](Reader& reader) { return reader.readAndParseTourism(); };
    getValue_NNmenuSmallGraph(readAndParseTourism);
}

void display_NNmenuMediumGraph() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Nearest Neighbour Medium Graph Menu       \n";
        cout << "-----------------------------\n";
        cout << "How many edges does the graph have:\n";
        cout << "1. 25\n";
        cout << "2. 50\n";
        cout << "3. 75\n";
        cout << "4. 100\n";
        cout << "5. 200\n";
        cout << "6. 300\n";
        cout << "7. 400\n";
        cout << "8. 500\n";
        cout << "9. 600\n";
        cout << "10. 700\n";
        cout << "11. 800\n";
        cout << "12. 900\n";
        cout << "e. Back to the other heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice == "e") {
            cout << "Exiting menu system...\n";
            exitMenu = true;
        } else {
            // Convert choice to an integer for comparison
            int choiceNum = stoi(choice);
            switch (choiceNum) {
                case 1:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
                    break;
                case 2:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_50.csv");
                    break;
                case 3:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_75.csv");
                    break;
                case 4:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_100.csv");
                    break;
                case 5:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_200.csv");
                    break;
                case 6:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_300.csv");
                    break;
                case 7:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_400.csv");
                    break;
                case 8:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_500.csv");
                    break;
                case 9:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_600.csv");
                    break;
                case 10:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_700.csv");
                    break;
                case 11:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_800.csv");
                    break;
                case 12:
                    getValue_NNmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_900.csv");
                    break;
                default:
                    cout << "Invalid input. Please choose a valid option.\n";
            }
        }
    }
}

void getValue_NNmenuMediumGraph(const std::string &filename){
    Reader reader;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.nearestNeighbour(graph);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double totalDistance = 0;
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i]->getInfo() << (i < tour.size() - 1 ? " -> " : "");
        if (i > 0) {
            Edge<int>* edge = graph.findEdge(tour[i - 1]->getInfo(), tour[i]->getInfo());
            if (edge)
                totalDistance += edge->getWeight();
        }
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
}

void display_LINmenu() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to LinKernighan Algorithm Heuristic       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Small Graph\n";
        cout << "2. Medium Graph \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display_LINmenuSmallGraph();
                break;
            case '2':
                display_LINmenuMediumGraph();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}
void display_LINmenuSmallGraph(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to LinKernighan Heuristic SmallGraph Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the graph you want:\n";
        cout << "1. Stadium Graph\n";
        cout << "2. Shipping Graph \n";
        cout <<"3. Tourism Graph \n";
        cout << "e. Back to the Other Heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                getValue_LINmenuSmallGraphStadium();
                break;
            case '2':
                getValue_LINmenuSmallGraphShipping();
                break;
            case '3':
                getValue_LINmenuSmallGraphTourism();
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void getValue_LINmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc) {
    Reader reader;
    Graph<int> graph;
    graph = readAndParseFunc(reader);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.linKernighan(graph);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double totalDistance = 0;
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i]->getInfo() << (i < tour.size() - 1 ? " -> " : "");
        if (i > 0) {
            Edge<int>* edge = graph.findEdge(tour[i - 1]->getInfo(), tour[i]->getInfo());
            if (edge)
                totalDistance += edge->getWeight();
        }
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";

}



void getValue_LINmenuSmallGraphStadium(){
    auto readAndParseStadium = [](Reader& reader) { return reader.readAndParseStadium(); };
    getValue_LINmenuSmallGraph(readAndParseStadium);
}

void getValue_LINmenuSmallGraphShipping() {
    auto readAndParseShipping = [](Reader& reader) { return reader.readAndParseShipping(); };
    getValue_LINmenuSmallGraph(readAndParseShipping);
}

void getValue_LINmenuSmallGraphTourism() {
    auto readAndParseTourism = [](Reader& reader) { return reader.readAndParseTourism(); };
    getValue_LINmenuSmallGraph(readAndParseTourism);
}

void display_LINmenuMediumGraph() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to LinKernighan Medium Graph Menu       \n";
        cout << "-----------------------------\n";
        cout << "How many edges does the graph have:\n";
        cout << "1. 25\n";
        cout << "2. 50\n";
        cout << "3. 75\n";
        cout << "4. 100\n";
        cout << "5. 200\n";
        cout << "6. 300\n";
        cout << "7. 400\n";
        cout << "8. 500\n";
        cout << "9. 600\n";
        cout << "10. 700\n";
        cout << "11. 800\n";
        cout << "12. 900\n";
        cout << "e. Back to the other Heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice == "e") {
            cout << "Exiting menu system...\n";
            exitMenu = true;
        } else {
            // Convert choice to an integer for comparison
            int choiceNum = stoi(choice);
            switch (choiceNum) {
                case 1:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
                    break;
                case 2:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_50.csv");
                    break;
                case 3:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_75.csv");
                    break;
                case 4:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_100.csv");
                    break;
                case 5:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_200.csv");
                    break;
                case 6:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_300.csv");
                    break;
                case 7:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_400.csv");
                    break;
                case 8:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_500.csv");
                    break;
                case 9:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_600.csv");
                    break;
                case 10:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_700.csv");
                    break;
                case 11:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_800.csv");
                    break;
                case 12:
                    getValue_LINmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_900.csv");
                    break;
                default:
                    cout << "Invalid input. Please choose a valid option.\n";
            }
        }
    }
}

void getValue_LINmenuMediumGraph(const std::string &filename){
    Reader reader;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.linKernighan(graph);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double totalDistance = 0;
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i]->getInfo() << (i < tour.size() - 1 ? " -> " : "");
        if (i > 0) {
            Edge<int>* edge = graph.findEdge(tour[i - 1]->getInfo(), tour[i]->getInfo());
            if (edge)
                totalDistance += edge->getWeight();
        }
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
}

void display_SANmenu() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Simulated Annealing Algorithm Heuristic       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Small Graph\n";
        cout << "2. Medium Graph \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display_SANmenuSmallGraph();
                break;
            case '2':
                display_SANmenuMediumGraph();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}
void display_SANmenuSmallGraph(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Simulated Annealing Heuristic SmallGraph Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the graph you want:\n";
        cout << "1. Stadium Graph\n";
        cout << "2. Shipping Graph \n";
        cout <<"3. Tourism Graph \n";
        cout << "e. Back to the other Heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                getValue_SANmenuSmallGraphStadium();
                break;
            case '2':
                getValue_SANmenuSmallGraphShipping();
                break;
            case '3':
                getValue_SANmenuSmallGraphTourism();
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void getValue_SANmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc) {
    Reader reader;
    Graph<int> graph;
    graph = readAndParseFunc(reader);
    auto start = std::chrono::high_resolution_clock::now();

    // Generate initial tour using Nearest Neighbor heuristic
    std::vector<Vertex<int>*> initialTour = graph.nearestNeighbour(graph);
    std::vector<int> currentTour;
    for (const auto& vertex : initialTour) {
        currentTour.push_back(vertex->getInfo());
    }

    // Set parameters for simulated annealing
    double initialTemperature = 1000.0;
    double coolingRate = 0.99;
    int iterations = 10000;

    // Perform Simulated Annealing
    std::vector<int> bestTour = simulatedAnnealing(graph, initialTemperature, coolingRate, iterations);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // Output the best tour and its cost
    double totalDistance = graph.computeTourLength(bestTour, graph);
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < bestTour.size(); ++i) {
        std::cout << bestTour[i] << (i < bestTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
}




void getValue_SANmenuSmallGraphStadium(){
    auto readAndParseStadium = [](Reader& reader) { return reader.readAndParseStadium(); };
    getValue_SANmenuSmallGraph(readAndParseStadium);
}

void getValue_SANmenuSmallGraphShipping() {
    auto readAndParseShipping = [](Reader& reader) { return reader.readAndParseShipping(); };
    getValue_SANmenuSmallGraph(readAndParseShipping);
}

void getValue_SANmenuSmallGraphTourism() {
    auto readAndParseTourism = [](Reader& reader) { return reader.readAndParseTourism(); };
    getValue_SANmenuSmallGraph(readAndParseTourism);
}

void display_SANmenuMediumGraph() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Simulated Annealing Medium Graph Menu       \n";
        cout << "-----------------------------\n";
        cout << "How many edges does the graph have:\n";
        cout << "1. 25\n";
        cout << "2. 50\n";
        cout << "3. 75\n";
        cout << "4. 100\n";
        cout << "5. 200\n";
        cout << "6. 300\n";
        cout << "7. 400\n";
        cout << "8. 500\n";
        cout << "9. 600\n";
        cout << "10. 700\n";
        cout << "11. 800\n";
        cout << "12. 900\n";
        cout << "e. Back to the other heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice == "e") {
            cout << "Exiting menu system...\n";
            exitMenu = true;
        } else {
            // Convert choice to an integer for comparison
            int choiceNum = stoi(choice);
            switch (choiceNum) {
                case 1:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
                    break;
                case 2:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_50.csv");
                    break;
                case 3:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_75.csv");
                    break;
                case 4:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_100.csv");
                    break;
                case 5:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_200.csv");
                    break;
                case 6:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_300.csv");
                    break;
                case 7:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_400.csv");
                    break;
                case 8:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_500.csv");
                    break;
                case 9:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_600.csv");
                    break;
                case 10:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_700.csv");
                    break;
                case 11:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_800.csv");
                    break;
                case 12:
                    getValue_SANmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_900.csv");
                    break;
                default:
                    cout << "Invalid input. Please choose a valid option.\n";
            }
        }
    }
}

std::vector<int> generateInitialTour(Graph<int> graph) {
    std::vector<int> tour;
    std::vector<Vertex<int>*> nearestNeighborTour = graph.nearestNeighbour(graph);

    // Convert Vertex pointers to their respective IDs
    for (const auto& vertex : nearestNeighborTour) {
        tour.push_back(vertex->getInfo());
    }

    return tour;
}

std::vector<int> generateNeighbor(const std::vector<int>& tour) {
    std::vector<int> neighbor = tour;

    // Select two random indices to swap
    int index1 = rand() % (neighbor.size() - 1); // Avoid swapping the start/end city
    int index2 = rand() % (neighbor.size() - 1);

    // Ensure the two indices are different
    while (index1 == index2) {
        index2 = rand() % (neighbor.size() - 1);
    }

    // Swap the cities at the selected indices
    std::swap(neighbor[index1], neighbor[index2]);

    return neighbor;
}


bool shouldAccept(double delta, double temperature) {
    // Decide whether to accept a new solution based on the change in cost and temperature
    return (delta < 0) || (exp(-delta / temperature) > static_cast<double>(rand()) / RAND_MAX);
}

std::vector<int> simulatedAnnealing(Graph<int>& graph, double initialTemperature, double coolingRate, int iterations) {
    std::vector<int> currentTour = generateInitialTour(graph);
    std::vector<int> bestTour = currentTour;
    double currentCost = graph.computeTourLength(currentTour, graph);
    double bestCost = currentCost;
    double temperature = initialTemperature;

    for (int i = 0; i < iterations; ++i) {
        std::vector<int> neighbor = generateNeighbor(currentTour);
        double neighborCost = graph.computeTourLength(neighbor, graph);
        double delta = neighborCost - currentCost;

        if (shouldAccept(delta, temperature)) {
            currentTour = neighbor;
            currentCost = neighborCost;
            if (currentCost < bestCost) {
                bestTour = currentTour;
                bestCost = currentCost;
            }
        }

        temperature *= coolingRate;
    }

    return bestTour;
}

void getValue_SANmenuMediumGraph(const std::string& filename) {
    Reader reader;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename);

    auto start = std::chrono::high_resolution_clock::now();

    // Generate initial tour using Nearest Neighbor heuristic
    std::vector<Vertex<int>*> initialTour = graph.nearestNeighbour(graph);
    std::vector<int> currentTour;
    for (const auto& vertex : initialTour) {
        currentTour.push_back(vertex->getInfo());
    }

    // Set parameters for simulated annealing
    double initialTemperature = 1000.0;
    double coolingRate = 0.99;
    int iterations = 10000;

    // Perform Simulated Annealing
    std::vector<int> bestTour = simulatedAnnealing(graph, initialTemperature, coolingRate, iterations);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // Output the best tour and its cost
    double totalDistance = graph.computeTourLength(bestTour, graph);
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < bestTour.size(); ++i) {
        std::cout << bestTour[i] << (i < bestTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
}

void display_HeldKarp_menu() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to HeldKarp Algorithm for Small Graphs optimal solution       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Small Graph\n";
        cout << "2. Medium Graph \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display_HeldKarp_menuSmallGraph();
                break;
            case '2':
                display_HeldKarp_menuMediumGraph();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}
void display_HeldKarp_menuSmallGraph(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to HeldKarp Algorithm for Small Graphs optimal solution SmallGraph Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the graph you want:\n";
        cout << "1. Stadium Graph\n";
        cout << "2. Shipping Graph \n";
        cout <<"3. Tourism Graph \n";
        cout << "e. Back to the other Heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                getValue_HeldKarp_menuSmallGraphStadium();
                break;
            case '2':
                getValue_HeldKarp_menuSmallGraphShipping();
                break;
            case '3':
                getValue_HeldKarp_menuSmallGraphTourism();
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void getValue_HeldKarp_menuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc) {
    Reader reader;
    Graph<int> graph;
    graph = readAndParseFunc(reader);
    std::vector<int> bestPath;
    double bestCost = std::numeric_limits<double>::infinity();
    std::vector<int> path(1, 0);
    graph.findVertex(0)->setVisited(true);

    auto startTime = std::chrono::high_resolution_clock::now();
    heldKarp(graph,bestPath,bestCost);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTime = endTime - startTime;

    std::cout << "Execution Time: " << execTime.count() << " seconds\n";
    std::cout << "Best Path Cost: " << bestCost << "\n";
    std::cout << "Best Path: ";
    for (size_t i = 0; i < bestPath.size(); ++i) {
        std::cout << bestPath[i];
        if (i < bestPath.size() - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << "\n";
}



void getValue_HeldKarp_menuSmallGraphStadium(){
    auto readAndParseStadium = [](Reader& reader) { return reader.readAndParseStadium(); };
    getValue_HeldKarp_menuSmallGraph(readAndParseStadium);
}

void getValue_HeldKarp_menuSmallGraphShipping() {
    auto readAndParseShipping = [](Reader& reader) { return reader.readAndParseShipping(); };
    getValue_HeldKarp_menuSmallGraph(readAndParseShipping);
}

void getValue_HeldKarp_menuSmallGraphTourism() {
    auto readAndParseTourism = [](Reader& reader) { return reader.readAndParseTourism(); };
    getValue_HeldKarp_menuSmallGraph(readAndParseTourism);
}

void getValue_HeldKarp_menuMediumGraph(const std::string &filename){
    Reader reader;
    Graph<int> graph;
    graph = reader.readAndParseExtra_Fully_Connected_Graphs(filename);
    std::vector<int> bestPath;
    double bestCost = std::numeric_limits<double>::infinity();
    std::vector<int> path(1, 0);
    graph.findVertex(0)->setVisited(true);

    auto startTime = std::chrono::high_resolution_clock::now();
    tsp(0, path, 0.0, 1, graph, bestPath, bestCost);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTime = endTime - startTime;

    std::cout << "Execution Time: " << execTime.count() << " seconds\n";
    std::cout << "Best Path Cost: " << bestCost << "\n";
    std::cout << "Best Path: ";
    for (int node : bestPath) {
        std::cout << node << (node == 0 ? "\n" : " -> ");
    }

}

void display_HeldKarp_menuMediumGraph() {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to HeldKarp Medium Graph Menu       \n";
        cout << "-----------------------------\n";
        cout << "How many edges does the graph have:\n";
        cout << "1. 25\n";
        cout << "2. 50\n";
        cout << "3. 75\n";
        cout << "4. 100\n";
        cout << "5. 200\n";
        cout << "6. 300\n";
        cout << "7. 400\n";
        cout << "8. 500\n";
        cout << "9. 600\n";
        cout << "10. 700\n";
        cout << "11. 800\n";
        cout << "12. 900\n";
        cout << "e. Back to the other heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice == "e") {
            cout << "Exiting menu system...\n";
            exitMenu = true;
        } else {
            // Convert choice to an integer for comparison
            int choiceNum = stoi(choice);
            switch (choiceNum) {
                case 1:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
                    break;
                case 2:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_50.csv");
                    break;
                case 3:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_75.csv");
                    break;
                case 4:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_100.csv");
                    break;
                case 5:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_200.csv");
                    break;
                case 6:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_300.csv");
                    break;
                case 7:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_400.csv");
                    break;
                case 8:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_500.csv");
                    break;
                case 9:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_600.csv");
                    break;
                case 10:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_700.csv");
                    break;
                case 11:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_800.csv");
                    break;
                case 12:
                    getValue_HeldKarp_menuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_900.csv");
                    break;
                default:
                    cout << "Invalid input. Please choose a valid option.\n";
            }
        }
    }
}

void display_CLUSTERmenu(int clusterOption) {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to K-means Clustering NearestNeighbour Algorithm Heuristic ( Mixes clustering and NN )      \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Small Graph\n";
        cout << "2. Medium Graph \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':\
                display_CLUSTERmenuSmallGraph(clusterOption);
                break;
            case '2':
                display_CLUSTERmenuMediumGraph(clusterOption);
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}


void display_CLUSTERmenuSmallGraph(int clusterOption) {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to K-means Clustering NearestNeighbour Heuristic SmallGraph Menu       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the graph you want:\n";
        cout << "1. Stadium Graph\n";
        cout << "2. Shipping Graph \n";
        cout << "3. Tourism Graph \n";
        cout << "e. Back to the other Heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                getValue_CLUSTERmenuSmallGraphStadium(clusterOption);
                break;
            case '2':
                getValue_CLUSTERmenuSmallGraphShipping(clusterOption);
                break;
            case '3':
                getValue_CLUSTERmenuSmallGraphTourism(clusterOption);
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}


void getValue_CLUSTERmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc, int clusterOption)  {
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    Reader reader;
    Graph<int> graph = readAndParseFunc(reader);
    auto coordinates = reader.readCoordinates();

    // Step 1: Clustering the graph
    int numClusters = clusterOption; // Number of clusters based on the user's choice
    auto clusters = reader.kMeansClustering(graph, numClusters, coordinates);

    // Step 2: Solve TSP for each cluster using nearest neighbor
    std::vector<int> combinedTour;
    std::unordered_set<int> visitedNodes;

    // Initially start from node 0
    int currentNode = 0;
    combinedTour.push_back(currentNode);
    visitedNodes.insert(currentNode);

    for (const auto& cluster : clusters) {
        if (cluster.empty()) continue;

        // Ensure the current node is included in the cluster
        std::vector<int> extendedCluster = cluster;
        if (std::find(cluster.begin(), cluster.end(), currentNode) == cluster.end()) {
            extendedCluster.push_back(currentNode);
        }

        auto clusterTour = graph.nearestNeighborForCluster(graph, extendedCluster);

        // Skip the starting node since it's already added
        for (size_t i = 1; i < clusterTour.size(); ++i) {
            if (visitedNodes.insert(clusterTour[i]).second) {
                combinedTour.push_back(clusterTour[i]);
                currentNode = clusterTour[i];
            }
        }
    }

    combinedTour.push_back(0); // End at the zero-identifier node

    // Calculate total distance of the combined tour
    double totalDistance = 0;
    int previous = combinedTour.front();
    for (size_t i = 1; i < combinedTour.size(); ++i) {
        Edge<int>* edge = graph.findEdge(previous, combinedTour[i]);
        if (edge) {
            totalDistance += edge->getWeight();
        }
        previous = combinedTour[i];
    }

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> execTime = end - start; // Calculate execution time

    // Print the TSP tour and total distance
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < combinedTour.size(); ++i) {
        std::cout << combinedTour[i] << (i < combinedTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Execution Time: " << execTime.count() << " seconds\n"; // Print execution time
}

void getValue_CLUSTERmenuSmallGraphStadium(int clusterOption) {
    auto readAndParseStadium = [](Reader& reader) { return reader.readAndParseStadium(); };
    getValue_CLUSTERmenuSmallGraph(readAndParseStadium, clusterOption);
}

void getValue_CLUSTERmenuSmallGraphShipping(int clusterOption) {
    auto readAndParseShipping = [](Reader& reader) { return reader.readAndParseShipping(); };
    getValue_CLUSTERmenuSmallGraph(readAndParseShipping, clusterOption);
}

void getValue_CLUSTERmenuSmallGraphTourism(int clusterOption) {
    auto readAndParseTourism = [](Reader& reader) { return reader.readAndParseTourism(); };
    getValue_CLUSTERmenuSmallGraph(readAndParseTourism, clusterOption);
}


void getValue_CLUSTERmenuMediumGraph(const std::string &filename, int clusterOption) {
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    Reader reader;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename);
    auto coordinates = reader.readCoordinates();

    // Step 1: Clustering the graph
    int numClusters = clusterOption; // Number of clusters based on the user's choice
    auto clusters = reader.kMeansClustering(graph, numClusters, coordinates);

    // Step 2: Solve TSP for each cluster using nearest neighbor
    std::vector<int> combinedTour;
    std::unordered_set<int> visitedNodes;

    // Initially start from node 0
    int currentNode = 0;
    combinedTour.push_back(currentNode);
    visitedNodes.insert(currentNode);

    for (const auto& cluster : clusters) {
        if (cluster.empty()) continue;

        // Ensure the current node is included in the cluster
        std::vector<int> extendedCluster = cluster;
        if (std::find(cluster.begin(), cluster.end(), currentNode) == cluster.end()) {
            extendedCluster.push_back(currentNode);
        }

        auto clusterTour = graph.nearestNeighborForCluster(graph, extendedCluster);

        // Skip the starting node since it's already added
        for (size_t i = 1; i < clusterTour.size(); ++i) {
            if (visitedNodes.insert(clusterTour[i]).second) {
                combinedTour.push_back(clusterTour[i]);
                currentNode = clusterTour[i];
            }
        }
    }

    combinedTour.push_back(0); // End at the zero-identifier node

    // Calculate total distance of the combined tour
    double totalDistance = 0;
    int previous = combinedTour.front();
    for (size_t i = 1; i < combinedTour.size(); ++i) {
        Edge<int>* edge = graph.findEdge(previous, combinedTour[i]);
        if (edge) {
            totalDistance += edge->getWeight();
        }
        previous = combinedTour[i];
    }

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> execTime = end - start; // Calculate execution time

    // Print the TSP tour and total distance
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < combinedTour.size(); ++i) {
        std::cout << combinedTour[i] << (i < combinedTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Execution Time: " << execTime.count() << " seconds\n"; // Print execution time

}

void display_CLUSTERmenuMediumGraph(int clusterOption) {
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Nearest Neighbour Medium Graph Menu       \n";
        cout << "-----------------------------\n";
        cout << "How many edges does the graph have:\n";
        cout << "1. 25\n";
        cout << "2. 50\n";
        cout << "3. 75\n";
        cout << "4. 100\n";
        cout << "5. 200\n";
        cout << "6. 300\n";
        cout << "7. 400\n";
        cout << "8. 500\n";
        cout << "9. 600\n";
        cout << "10. 700\n";
        cout << "11. 800\n";
        cout << "12. 900\n";
        cout << "e. Back to the other heuristics Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice == "e") {
            cout << "Exiting menu system...\n";
            exitMenu = true;
        } else {
            // Convert choice to an integer for comparison
            int choiceNum = stoi(choice);
            switch (choiceNum) {
                case 1:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_25.csv", clusterOption);
                    break;
                case 2:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_50.csv", clusterOption);
                    break;
                case 3:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_75.csv", clusterOption);
                    break;
                case 4:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_100.csv", clusterOption);
                    break;
                case 5:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_200.csv", clusterOption);
                    break;
                case 6:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_300.csv", clusterOption);
                    break;
                case 7:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_400.csv", clusterOption);
                    break;
                case 8:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_500.csv", clusterOption);
                    break;
                case 9:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_600.csv", clusterOption);
                    break;
                case 10:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_700.csv", clusterOption);
                    break;
                case 11:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_800.csv", clusterOption);
                    break;
                case 12:
                    getValue_CLUSTERmenuMediumGraph("../Data/Extra_Fully_Connected_Graphs/edges_900.csv", clusterOption);
                    break;
                default:
                    cout << "Invalid input. Please choose a valid option.\n";
            }
        }
    }
}

void display_RWmenu(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Real World TSP       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Graph1('Small')\n";
        cout << "2. Graph2(Medium) \n";
        cout << "3. Graph3(Large)  \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;

        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                // getValue_RWsmallGraph();
                break;
            case '2':
                // getValue_RWmediumGraph();
                break;
            case '3':
                // getValue_RWlargeGraph();
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
}

void getValue_RWsmallGraph(){
    Reader reader;
    Graph<int> graph = reader.readAndParseRealWorld_Graphs(1);
    Vertex<int>* startVertexPtr = graph.findVertex(0);

    vector<Edge<int>*> mst = graph.primMST(startVertexPtr->getInfo());

    vector<int> preorderList = preOrderTraversal(startVertexPtr, mst);

    vector<int> tspTour;
    set<int> visited;
    for (int vertex : preorderList) {
        if (visited.insert(vertex).second) {
            tspTour.push_back(vertex);
        }
    }
    tspTour.push_back(tspTour.front());

    double totalDistance = 0;
    int previous = tspTour.front();
    for (size_t i = 1; i < tspTour.size(); i++) {
        Edge<int> *edge = graph.findEdge(previous, tspTour[i]);
        if (edge)
            totalDistance += edge->getWeight();
        previous = tspTour[i];
    }

    cout << "TSP Tour: ";
    for (size_t i = 0; i < tspTour.size(); ++i) {
        cout << tspTour[i] << (i < tspTour.size() - 1 ? " -> " : "");
    }
    cout << endl;
    cout << "Total Approximation Distance: " << totalDistance << "\n";
}

void App::run() {
    mainMenu();
}