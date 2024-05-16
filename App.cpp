#include "App.h"


#include <iostream>
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
void display4_2menuLargeGraph();
void display4_2menuSmallGraphStadium();
void display4_2menuSmallGraphShipping();
void display4_2menuSmallGraphTourism();
void display4_2menuMedium(const std::string &filename);
void display4_2menuLarge1(int nodeID);
void display4_2menuLarge2(int nodeID);
void display4_2menuLarge3(int nodeID);

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

void getValue_RWsmallGraph(int nodeID);
void getValue_RWmediumGraph(int nodeID);
void getValue_RWlargeGraph(int nodeID);

std::vector<int> preOrderTraversal(Vertex<int>* root, const std::vector<Edge<int>*>& mst);
void tsp(int currentNode, std::vector<int>& path, double currentCost, int level, Graph<int>& graph, std::vector<int>& bestPath, double& bestCost);
void heldKarp(const Graph<int>& graph, std::vector<int>& bestPath, double& bestCost);


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
            case '4':
                display_RWmenu();
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
        cout << "2. Medium Graph\n";
        cout << "3  Large Graph\n";
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
            case '3':
                display4_2menuLargeGraph();
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
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename,vertexMap,edgeMap);
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

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTime = end - start;

    cout << "TSP Tour: ";
    for (size_t i = 0; i < tspTour.size(); ++i) {
        cout << tspTour[i] << (i < tspTour.size() - 1 ? " -> " : "");
    }
    cout << endl;
    cout << "Total Approximation Distance: " << totalDistance << "\n";
    cout << "Execution Time: " << execTime.count() << " seconds\n";
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
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename,vertexMap,edgeMap);
    Vertex<int>* startVertexPtr = nullptr;
    auto vertexIter = vertexMap.find(0);
    if (vertexIter != vertexMap.end()) {
        startVertexPtr = vertexIter->second;
    }

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

void display4_2menuLargeGraph() {
    string choice_str;
    int choice;
    bool exitMenu = false;
    while (!exitMenu) {
        cout << "\n-----------------------------\n";
        cout << "     Welcome to Triangular Approximation for Large Graphs       \n";
        cout << "-----------------------------\n";
        cout << "Enter the number of the option of the size of the graph you want:\n";
        cout << "1. Graph1('Small')\n";
        cout << "2. Graph2(Medium) \n";
        cout << "3. Graph3(Large)  \n";
        cout << "e. Back to the main Menu\n";
        cout << "-----------------------------\n";
        cout << "Your choice: ";
        cin >> choice_str;

        try {
            choice = stoi(choice_str);
        } catch (const invalid_argument&) {
            choice = 0;
        }

        int nodeID = 0;
        string choiceNode;

        switch (choice) {
            case 1:
                while (true) {
                    cout << "Input the starting node for the algorithm: ";
                    cin >> choiceNode;
                    try {
                        nodeID = stoi(choiceNode);
                        display4_2menuLarge1(nodeID);
                        break;
                    } catch (const invalid_argument&) {
                        cout << "Invalid input. Please enter a valid integer.\n";
                    } catch (const out_of_range&) {
                        cout << "Input out of range. Please enter a valid integer within the range of int.\n";
                    }
                }
                break;
            case 2:
                while (true) {
                    cout << "Input the starting node for the algorithm: ";
                    cin >> choiceNode;
                    try {
                        nodeID = stoi(choiceNode);
                        display4_2menuLarge2(nodeID);
                        break;
                    } catch (const invalid_argument&) {
                        cout << "Invalid input. Please enter a valid integer.\n";
                    } catch (const out_of_range&) {
                        cout << "Input out of range. Please enter a valid integer within the range of int.\n";
                    }
                }
                break;
            case 3:
                while (true) {
                    cout << "Input the starting node for the algorithm: ";
                    cin >> choiceNode;
                    try {
                        nodeID = stoi(choiceNode);
                        display4_2menuLarge3(nodeID);
                        break;
                    } catch (const invalid_argument&) {
                        cout << "Invalid input. Please enter a valid integer.\n";
                    } catch (const out_of_range&) {
                        cout << "Input out of range. Please enter a valid integer within the range of int.\n";
                    }
                }
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

void display4_2menuLarge1(int nodeID) {
    Reader reader;
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    Graph<int> graph = reader.readAndParseRealWorld_Graphs(1,vertexMap,edgeMap);
    Vertex<int>* startVertexPtr = nullptr;
    auto vertexIter = vertexMap.find(nodeID);
    if (vertexIter != vertexMap.end()) {
        startVertexPtr = vertexIter->second;
    } else {
        std::cerr << "Node ID " << nodeID << " not found in the graph.\n";
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();

    vector<Edge<int>*> mst = graph.primMSTMaps(startVertexPtr->getInfo(),vertexMap,edgeMap);

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

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    totalDistance = graph.computeTourLength(tspTour, graph);
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < tspTour.size(); ++i) {
        std::cout << tspTour[i] << (i < tspTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
}

void display4_2menuLarge2(int nodeID) {
    Reader reader;
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    Graph<int> graph = reader.readAndParseRealWorld_Graphs(2,vertexMap,edgeMap);
    Vertex<int>* startVertexPtr = nullptr;
    auto vertexIter = vertexMap.find(nodeID);
    if (vertexIter != vertexMap.end()) {
        startVertexPtr = vertexIter->second;
    } else {
        std::cerr << "Node ID " << nodeID << " not found in the graph.\n";
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();

    vector<Edge<int>*> mst = graph.primMSTMaps(startVertexPtr->getInfo(),vertexMap,edgeMap);

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

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    totalDistance = graph.computeTourLength(tspTour, graph);
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < tspTour.size(); ++i) {
        std::cout << tspTour[i] << (i < tspTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
}

void display4_2menuLarge3(int nodeID) {
    Reader reader;
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    Graph<int> graph = reader.readAndParseRealWorld_Graphs(3,vertexMap,edgeMap);
    Vertex<int>* startVertexPtr = nullptr;
    auto vertexIter = vertexMap.find(nodeID);
    if (vertexIter != vertexMap.end()) {
        startVertexPtr = vertexIter->second;
    } else {
        std::cerr << "Node ID " << nodeID << " not found in the graph.\n";
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();

    vector<Edge<int>*> mst = graph.primMSTMaps(startVertexPtr->getInfo(),vertexMap,edgeMap);

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

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    totalDistance = graph.computeTourLength(tspTour, graph);
    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < tspTour.size(); ++i) {
        std::cout << tspTour[i] << (i < tspTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
}

/**
 * @brief Performs a pre-order traversal of a tree represented by an MST.
 * @param root The root vertex of the tree.
 * @param mst The minimum spanning tree (MST) represented by a vector of edges.
 * @return A vector containing the vertices visited in pre-order traversal.
 * @timecomplexity The time complexity of this function depends on the size of the MST and the efficiency of the pre-order traversal algorithm.
 * Typically, it has a time complexity of O(V + E), where V is the number of vertices and E is the number of edges in the MST.
 */
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

/**
 * @brief Solves the Traveling Salesman Problem (TSP) using a backtracking algorithm.
 * @param currentNode The current node being visited.
 * @param path The current path being explored.
 * @param currentCost The current cost of the path.
 * @param level The level of recursion.
 * @param graph The graph representing the TSP problem.
 * @param bestPath The best path found so far.
 * @param bestCost The cost of the best path found so far.
 * @details This function recursively explores possible paths using a backtracking algorithm and updates the best path and cost found.
 * @timecomplexity The time complexity of this function depends on the size of the graph and the efficiency of the backtracking algorithm used.
 * In the worst case, it has an exponential time complexity.
 */
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

/**
 * @brief Solves the Traveling Salesman Problem (TSP) using the Held-Karp algorithm.
 * @param graph The graph representing the TSP problem.
 * @param bestPath The best path found by the algorithm.
 * @param bestCost The cost of the best path found.
 * @details This function applies the Held-Karp algorithm to find the optimal solution for the TSP.
 * @timecomplexity The time complexity of this function is O(2^n * n^2), where n is the number of vertices in the graph.
 * @remarks Pros: Guarantees the optimal solution for the TSP. Suitable for small to medium-sized graphs.
 * Cons: High time complexity, impractical for large graphs due to exponential growth.
 */
void heldKarp(const Graph<int>& graph, std::vector<int>& bestPath, double& bestCost) {
    int n = graph.getNumVertex();

    std::vector<std::vector<double>> dp(1 << n, std::vector<double>(n, INF));
    std::vector<std::vector<int>> parent(1 << n, std::vector<int>(n, -1));
    dp[1][0] = 0;

    for (int mask = 1; mask < (1 << n); ++mask) {
        for (int u = 0; u < n; ++u) {
            if (mask & (1 << u)) {
                for (int v = 0; v < n; ++v) {
                    if (!(mask & (1 << v))) {
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
        bestPath.push_back(0);
        std::reverse(bestPath.begin(), bestPath.end());
        bestPath.push_back(0);
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
        cout << "3. LinKernighan\n";
        cout << "4. HeldKarp (Optimal solution but feasible only on toy graphs)";
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
                display_LINmenu();
                break;
            case '4':
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

/**
 * Displays a menu for selecting small graph instances and applying the Nearest Neighbor algorithm.
 */
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

/**
 * Applies the Nearest Neighbor algorithm to the stadium graph instance.
 */
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

/**
 * Applies the Nearest Neighbor algorithm to the shipping graph instance.
 */
void getValue_NNmenuSmallGraphStadium(){
    auto readAndParseStadium = [](Reader& reader) { return reader.readAndParseStadium(); };
    getValue_NNmenuSmallGraph(readAndParseStadium);
}

/**
 * Applies the Nearest Neighbor algorithm to the tourism graph instance.
 */
void getValue_NNmenuSmallGraphShipping() {
    auto readAndParseShipping = [](Reader& reader) { return reader.readAndParseShipping(); };
    getValue_NNmenuSmallGraph(readAndParseShipping);
}

/**
 * Applies the Nearest Neighbor algorithm to the tourism graph instance.
 */
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
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename,vertexMap,edgeMap);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.nearestNeighbourMedium(graph,vertexMap,edgeMap);
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
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.linKernighan(graph,vertexMap,edgeMap);
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
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename,vertexMap,edgeMap);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.linKernighan(graph,vertexMap,edgeMap);
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
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename,vertexMap,edgeMap);
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
    auto start = std::chrono::high_resolution_clock::now();

    Reader reader;
    Graph<int> graph = readAndParseFunc(reader);
    auto coordinates = reader.readCoordinates();

    int numClusters = clusterOption;
    auto clusters = reader.kMeansClustering(graph, numClusters, coordinates);

    std::vector<int> combinedTour;
    std::unordered_set<int> visitedNodes;

    int currentNode = 0;
    combinedTour.push_back(currentNode);
    visitedNodes.insert(currentNode);

    for (const auto& cluster : clusters) {
        if (cluster.empty()) continue;

        std::vector<int> extendedCluster = cluster;
        if (std::find(cluster.begin(), cluster.end(), currentNode) == cluster.end()) {
            extendedCluster.push_back(currentNode);
        }

        auto clusterTour = graph.nearestNeighborForCluster(graph, extendedCluster);

        for (size_t i = 1; i < clusterTour.size(); ++i) {
            if (visitedNodes.insert(clusterTour[i]).second) {
                combinedTour.push_back(clusterTour[i]);
                currentNode = clusterTour[i];
            }
        }
    }

    combinedTour.push_back(0);

    double totalDistance = 0;
    int previous = combinedTour.front();
    for (size_t i = 1; i < combinedTour.size(); ++i) {
        Edge<int>* edge = graph.findEdge(previous, combinedTour[i]);
        if (edge) {
            totalDistance += edge->getWeight();
        }
        previous = combinedTour[i];
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTime = end - start;

    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < combinedTour.size(); ++i) {
        std::cout << combinedTour[i] << (i < combinedTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Execution Time: " << execTime.count() << " seconds\n";
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
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;
    Graph<int> graph = reader.readAndParse4_2Extra_Fully_Connected_Graphs(filename,vertexMap,edgeMap);
    auto coordinates = reader.readCoordinates();

    int numClusters = clusterOption;
    auto clusters = reader.kMeansClustering(graph, numClusters, coordinates);

    std::vector<int> combinedTour;
    std::unordered_set<int> visitedNodes;

    int currentNode = 0;
    combinedTour.push_back(currentNode);
    visitedNodes.insert(currentNode);

    for (const auto& cluster : clusters) {
        if (cluster.empty()) continue;

        std::vector<int> extendedCluster = cluster;
        if (std::find(cluster.begin(), cluster.end(), currentNode) == cluster.end()) {
            extendedCluster.push_back(currentNode);
        }

        auto clusterTour = graph.nearestNeighborForCluster(graph, extendedCluster);

        for (size_t i = 1; i < clusterTour.size(); ++i) {
            if (visitedNodes.insert(clusterTour[i]).second) {
                combinedTour.push_back(clusterTour[i]);
                currentNode = clusterTour[i];
            }
        }
    }

    combinedTour.push_back(0);

    double totalDistance = 0;
    int previous = combinedTour.front();
    for (size_t i = 1; i < combinedTour.size(); ++i) {
        Edge<int>* edge = graph.findEdge(previous, combinedTour[i]);
        if (edge) {
            totalDistance += edge->getWeight();
        }
        previous = combinedTour[i];
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTime = end - start;

    std::cout << "TSP Tour: ";
    for (size_t i = 0; i < combinedTour.size(); ++i) {
        std::cout << combinedTour[i] << (i < combinedTour.size() - 1 ? " -> " : "");
    }
    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Execution Time: " << execTime.count() << " seconds\n";

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

void display_RWmenu() {
    string choice_str;
    int choice;
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
        cin >> choice_str;

        try {
            choice = stoi(choice_str);
        } catch (const invalid_argument&) {
            choice = 0;
        }

        int nodeID = 0;
        string choiceNode;

        switch (choice) {
            case 1:
                while (true) {
                    cout << "Input the starting node for the algorithm: ";
                    cin >> choiceNode;
                    try {
                        nodeID = stoi(choiceNode);
                        getValue_RWsmallGraph(nodeID);
                        break;
                    } catch (const invalid_argument&) {
                        cout << "Invalid input. Please enter a valid integer.\n";
                    } catch (const out_of_range&) {
                        cout << "Input out of range. Please enter a valid integer within the range of int.\n";
                    }
                }
                break;
            case 2:
                while (true) {
                    cout << "Input the starting node for the algorithm: ";
                    cin >> choiceNode;
                    try {
                        nodeID = stoi(choiceNode);
                        getValue_RWmediumGraph(nodeID);
                        break;
                    } catch (const invalid_argument&) {
                        cout << "Invalid input. Please enter a valid integer.\n";
                    } catch (const out_of_range&) {
                        cout << "Input out of range. Please enter a valid integer within the range of int.\n";
                    }
                }
                break;
            case 3:
                while (true) {
                    cout << "Input the starting node for the algorithm: ";
                    cin >> choiceNode;
                    try {
                        nodeID = stoi(choiceNode);
                        getValue_RWlargeGraph(nodeID);
                        break;
                    } catch (const invalid_argument&) {
                        cout << "Invalid input. Please enter a valid integer.\n";
                    } catch (const out_of_range&) {
                        cout << "Input out of range. Please enter a valid integer within the range of int.\n";
                    }
                }
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

void getValue_RWsmallGraph(int nodeID){

    Reader reader;
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;

    auto start1 = std::chrono::high_resolution_clock::now();
    Graph<int> graph = reader.readAndParseRealWorld_Graphs(1, vertexMap, edgeMap);
    auto end1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration1 = end1 - start1;

    cout << "Parsing durou: " << duration1.count() << "\n";

    Vertex<int>* startingVertex = vertexMap[nodeID];
    if(!startingVertex){
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.nearestNeighbourNode(graph,startingVertex, vertexMap, edgeMap);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;


    double totalDistance = 0;
    cout << "TSP Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        cout << tour[i]->getInfo() << (i < tour.size() - 1 ? " -> " : "");
        if (i > 0) {
            Edge<int>* edge = graph.findEdge(tour[i - 1]->getInfo(), tour[i]->getInfo());
            if (edge)
                totalDistance += edge->getWeight();
        }
    }

    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
    std::cout << "Tour Size: " << tour.size() << "\n";

}

void getValue_RWmediumGraph(int nodeID) {
    Reader reader;
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;

    auto start1 = std::chrono::high_resolution_clock::now();
    Graph<int> graph = reader.readAndParseRealWorld_Graphs(2, vertexMap, edgeMap);
    auto end1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration1 = end1 - start1;

    cout << "Parsing durou: " << duration1.count() << "\n";

    Vertex<int>* startingVertex = vertexMap[nodeID];
    if(!startingVertex){
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.nearestNeighbourNode(graph,startingVertex, vertexMap, edgeMap);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;


    double totalDistance = 0;
    cout << "TSP Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        cout << tour[i]->getInfo() << (i < tour.size() - 1 ? " -> " : "");
        if (i > 0) {
            Edge<int>* edge = graph.findEdge(tour[i - 1]->getInfo(), tour[i]->getInfo());
            if (edge)
                totalDistance += edge->getWeight();
        }
    }

    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
    std::cout << "Tour Size: " << tour.size() << "\n";
}

void getValue_RWlargeGraph(int nodeID) {
    Reader reader;
    unordered_map<int, Vertex<int>*> vertexMap;
    unordered_map<std::string, Edge<int>*> edgeMap;

    auto start1 = std::chrono::high_resolution_clock::now();
    Graph<int> graph = reader.readAndParseRealWorld_Graphs(3, vertexMap, edgeMap);
    auto end1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration1 = end1 - start1;

    cout << "Parsing durou: " << duration1.count() << "\n";
    Vertex<int>* startingVertex = vertexMap[nodeID];
    if(!startingVertex){
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vertex<int>*> tour = graph.nearestNeighbourNode(graph,startingVertex, vertexMap, edgeMap);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;


    double totalDistance = 0;
    cout << "TSP Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        cout << tour[i]->getInfo() << (i < tour.size() - 1 ? " -> " : "");
        if (i > 0) {
            Edge<int>* edge = graph.findEdge(tour[i - 1]->getInfo(), tour[i]->getInfo());
            if (edge)
                totalDistance += edge->getWeight();
        }
    }

    std::cout << std::endl;
    std::cout << "Total Approximation Distance: " << totalDistance << "\n";
    std::cout << "Time: " << duration.count() << "\n";
    std::cout << "Tour Size: " << tour.size() << "\n";
}

void App::run() {
    mainMenu();
}