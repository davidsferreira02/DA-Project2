//
// Created by David Ferreira on 11/05/2024.
//

#include "App.h"


#include <iostream>
#include <cmath>
#include <utility>
#include "reader.h"
#include "App.h"

using namespace std;



using namespace std::chrono;


//Backtracking

int n;
vector<vector<double>> graph;
vector<double> best_path;
double best_cost = numeric_limits<double>::infinity();

bool is_valid(vector<double>& path, int pos, int city) {
    // Check if the city has not been visited in the current path
    for (int i = 0; i < pos; ++i) {
        if (path[i] == city)
            return false;
    }
    return true;
}

void recursive_backtracking(vector<double>& path, int pos, double cost) {
    if (pos == n) {

        // Check if the path is a complete solution
        if (cost + graph[path[pos - 1]][0] < best_cost) {
            best_cost = cost + graph[path[pos - 1]][0];
            best_path = path;

        }
        return;
    }
    for (int i = 1; i < n; ++i) {
        if (is_valid(path, pos, i)) {
            path[pos] = i;
            recursive_backtracking(path, pos + 1, cost + graph[path[pos - 1]][i]);
        }
    }
}

vector<double> tsp(vector<vector<double>>& distance_matrix) {
    n = distance_matrix.size();
    graph = distance_matrix;
    vector<double> path(n);
    path[0] = 0; // Starting from city 0
    recursive_backtracking(path, 1, 0);
    return best_path;

}
void display4_1menuSmallGraphStadium(){
    Reader reader;
    std::vector<std::vector<double>> matrix=reader.readAndParseStadium();
    auto start = high_resolution_clock::now(); // Start timing
    vector<double> solution = tsp(matrix);
    auto stop = high_resolution_clock::now(); // Stop timing
    auto duration = duration_cast<milliseconds>(stop - start); // Calculate duration

    cout << "Melhor caminho encontrado: ";
    for (int city : solution) {
        cout << city << " ";
    }
    cout << endl;

    cout << "Custo do melhor caminho: " << best_cost << endl;
    cout << "Tempo de execução: " << duration.count() << " milliseconds" << endl;
}

void display4_1menuSmallGraphShipping(){
    Reader reader;
    std::vector<std::vector<double>> matrix=reader.readAndParseShipping();
    auto start = high_resolution_clock::now(); // Start timing
    vector<double> solution = tsp(matrix);
    auto stop = high_resolution_clock::now(); // Stop timing
    auto duration = duration_cast<milliseconds>(stop - start); // Calculate duration

    cout << "Melhor caminho encontrado: ";
    for (int city : solution) {
        cout << city << " ";
    }
    cout << endl;

    cout << "Custo do melhor caminho: " << best_cost << endl;
    cout << "Tempo de execução: " << duration.count() << " milliseconds" << endl;
}

void display4_1menuSmallGraphTourism()
{
    Reader reader;
    std::vector<std::vector<double>> matrix=reader.readAndParseTourism();
    auto start = high_resolution_clock::now(); // Start timing
    vector<double> solution = tsp(matrix);
    auto stop = high_resolution_clock::now(); // Stop timing
    auto duration = duration_cast<milliseconds>(stop - start); // Calculate duration

    cout << "Melhor caminho encontrado: ";
    for (int city : solution) {
        cout << city << " ";
    }
    cout << endl;

    cout << "Custo do melhor caminho: " << best_cost << endl;
    cout << "Tempo de execução: " << duration.count() << " milliseconds" << endl;
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
                cout << "Exiting menu system...\n";
                exitMenu = true;
                break;
            default:
                cout << "Invalid input. Please choose a valid option.\n";
        }
    }
    return;

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
    return 0;
}


void App::run() {
    mainMenu();
}

