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

int bestCost;
vector<int> bestPath;

bool is_valid(vector<int>& path, int pos, int city) {
    for (int i = 0; i < pos; ++i) {
        if (path[i] == city)
            return false;
    }
    return true;
}

 void recursive_backtracking(vector<int>& path, int pos, int cost,vector<vector<int>>& distance_matrix) {
    if (pos == distance_matrix.size()) {

        if (cost + distance_matrix[path[pos - 1]][0] < bestCost) {
            bestCost = cost + distance_matrix[path[pos - 1]][0];
            bestPath = path;
        }
       return;
    }

    for (int i = 1; i < distance_matrix.size(); ++i) {
        if (is_valid(path, pos, i)) {
            path[pos] = i;
            recursive_backtracking(path, pos + 1, cost + distance_matrix[path[pos - 1]][i],distance_matrix);
        }
    }
}

vector<int> tsp(vector<vector<int>>& distance_matrix) {
    int n = distance_matrix.size();
    vector<int> path(n);
    path[0] = 0; // Começando da cidade 0
    recursive_backtracking(path, 1, 0,distance_matrix);
    return bestPath;
}

void display4_1menuSmallGraphStadium(){
    Reader reader;
    reader.readAndParseStadium();
    int best_cost=99999999;
    Graph<int*> graph = reader.getGraph();
    vector<vector<int>> distance_matrix;
    int i=0;
    int j=0;
    for(auto v:graph.getVertexSet()){
        for(auto e:v->getAdj()){
            distance_matrix[i][j]=e->getWeight();
            j++;
        }
        i++;
    }
    vector<int> solution = tsp(distance_matrix);


    cout << "Melhor caminho encontrado: ";
    for (int node : solution) {
        cout << node << " ";
    }
    cout << endl;

    cout << "Custo do melhor caminho: " << best_cost << endl;


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
        cout <<"3. tourism Graph \n";
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
void display4_1menuSmallGraph();



int mainMenu(){
    cout << "Loading...";

    Reader reader;


    Graph<int*>  graph = reader.getGraph();

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

