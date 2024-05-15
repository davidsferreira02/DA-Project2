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
#include "App.h"

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
std::vector<int> preOrderTraversal(Vertex<int>* root, const std::vector<Edge<int>*>& mst);
void tsp(int currentNode, std::vector<int>& path, double currentCost, int level, Graph<int>& graph, std::vector<int>& bestPath, double& bestCost);
void heldKarp(const Graph<int>& graph, std::vector<int>& bestPath, double& bestCost);
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

void display_LINmenu();
void display_LINmenuSmallGraph();
void display_LINmenuMediumGraph();
void getValue_LINmenuSmallGraph(const std::function<Graph<int>(Reader&)>& readAndParseFunc);
void getValue_LINmenuSmallGraphStadium();
void getValue_LINmenuSmallGraphShipping();
void getValue_LINmenuSmallGraphTourism();
void getValue_LINmenuMediumGraph(const std::string &filename);





int mainMenu(){
    cout << "Loading...";

    Reader reader;


    // Graph<int*>  graph = reader.getGraph();

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
        cout <<"3. Large Graph \n";
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
    return;

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
        cout <<"3. Large Graph \n";
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
    return;

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
    return;

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
    for (int node : bestPath) {
        std::cout << node << (node == 0 ? "\n" : " -> ");
    }
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

void display4_1menuMediumGraph(){
    string choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Backtracing Medium Graph Menu       \n";
        cout << "-----------------------------\n";
        cout << "how many edges the graph has:\n";
        cout << "1. 25\n";
        cout << "2. 50 \n";
        cout <<"3. 75 \n";
        cout <<"4. 100 \n";
        cout <<"5. 200 \n";
        cout <<"6. 300 \n";
        cout <<"7. 400 \n";
        cout <<"8. 500\n";
        cout <<"9. 600 \n";
        cout <<"10. 700 \n";
        cout <<"11. 800 \n";
        cout <<"12. 900 \n";
        cout << "e. Back to the backtracing Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice;
        if (choice.length() != 1) {
            choice = "0";
        }

        switch (choice[0]) {
            case '1':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
                break;
            case '2':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_50.csv");
                break;
            case '3':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_75.csv");
                break;
            case '4':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_100.csv");
                break;
            case '5':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_200.csv");
                break;
            case '6':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_300.csv");
                break;
            case '7':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_400.csv");
                break;
            case '8':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_500.csv");
                break;
            case '9':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_600.csv");
                break;
            case '10':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_700.csv");
                break;
            case '11':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_800.csv");
                break;
            case '12':
                display4_1menuMedium("../Data/Extra_Fully_Connected_Graphs/edges_900.csv");
                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
    return;

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
    return;

}


void display4_2menuSmall(const std::function<Graph<int>(Reader&)>& readAndParseFunc) {
    Reader reader;
    Graph<int> graph = readAndParseFunc(reader);
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
    dp[1][0] = 0; // Starting at node 0.

    // Iterate over all subsets of nodes.
    for (int mask = 1; mask < (1 << n); ++mask) {
        for (int u = 0; u < n; ++u) {
            if (mask & (1 << u)) { // If u is in the subset represented by mask.
                for (int v = 0; v < n; ++v) {
                    if (!(mask & (1 << v))) { // If v is not in the subset.
                        Edge<int>* edge = graph.findEdge(u, v);
                        if (edge) {
                            dp[mask | (1 << v)][v] = std::min(dp[mask | (1 << v)][v], dp[mask][u] + edge->getWeight());
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

    // Reconstruct the path using the last node.
    if (lastNode != -1) {
        bestPath.clear();
        int mask = (1 << n) - 1;
        int currentNode = lastNode;
        while (currentNode != 0) {
            bestPath.push_back(currentNode);
            for (int v = 0; v < n; ++v) {
                if ((mask & (1 << v)) && graph.findEdge(v, currentNode)) {
                    mask ^= (1 << currentNode);
                    currentNode = v;
                    break;
                }
            }
        }
        bestPath.push_back(0); // Add the starting node to complete the cycle.
        std::reverse(bestPath.begin(), bestPath.end()); // Reverse the path to start from the starting node.
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
        cout << "2. LinKernighan\n";
        cout <<"3.  \n";
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
            case '2':
                display_LINmenu();
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
    return;

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
        cout <<"3. Large Graph \n";
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
            case '3':

                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
    return;

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
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
    return;

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
        cout <<"3. Large Graph \n";
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
            case '3':

                break;
            case 'e':
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
    return;

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
    return;

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
void App::run() {
    mainMenu();
}