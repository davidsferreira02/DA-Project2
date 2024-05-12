#include "reader.h"
#include "Classes/Shipping.h"


/**
 * @brief Read and parse data from the Stations.csv file to create Station objects.
 *
 * This function reads data from the Stations.csv file, creates Station objects, and adds them to the graph.
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>


std::vector<std::vector<double>> Reader::readAndParseStadium() {
    std::ifstream file("../Data/Toy-Graphs/stadiums.csv");
    std::string line;

    std::unordered_set<int> nodes;
    std::vector<Stadium*> stadiums;
    getline(file, line); // Skip the header line
    while (getline(file, line)) {
        std::stringstream ss(line);
        int source, destination;
        double dist;
        char comma;
        ss >> source >> comma >> destination >> comma >> dist;
        nodes.insert(source);
        nodes.insert(destination);
        stadiums.push_back(new Stadium(source, destination, dist));
    }

    int num_nodes = nodes.size();
    std::vector<std::vector<double>> matrix(num_nodes, std::vector<double>(num_nodes, 0.0));

    for (auto s : stadiums) {
        matrix[s->getSrc()][s->getDst()] = s->getDist();
    }

    // Clean up dynamically allocated memory
    for (auto s : stadiums) {
        delete s;
    }

    return matrix;
}

std::vector<std::vector<double>> Reader::readAndParseShipping() {
    std::ifstream file("../Data/Toy-Graphs/shipping.csv");
    std::string line;

    std::unordered_set<int> nodes;
    std::vector<Shipping*> shipping;
    getline(file, line); // Skip the header line
    while (getline(file, line)) {
        std::stringstream ss(line);
        int source, destination;
        double dist;
        char comma;
        ss >> source >> comma >> destination >> comma >> dist;
        nodes.insert(source);
        nodes.insert(destination);
        shipping.push_back(new Shipping(source, destination, dist));
    }

    int num_nodes = nodes.size();
    std::vector<std::vector<double>> matrix(num_nodes, std::vector<double>(num_nodes, 0.0));

    for (auto s : shipping) {
        matrix[s->getSrc()][s->getDst()] = s->getDist();
    }

    // Clean up dynamically allocated memory
    for (auto s : shipping) {
        delete s;
    }

    return matrix;
}



std::vector<std::vector<double>> Reader::readAndParseTourism() {
    std::ifstream file("../Data/Toy-Graphs/tourism.csv");
    std::string line;

    std::unordered_set<int> nodes;
    std::vector<Tourism*> tourism;
   // getline(file, line); // Skip the header line
    while (getline(file, line)) {
        std::stringstream ss(line);
        int source, destination;
        double dist;
        std::string label_org;
        std::string label_dst;
        char comma;
        ss >> source >> comma >> destination >> comma >> dist>> comma>>label_org >>comma >> label_dst;
        nodes.insert(source);
        nodes.insert(destination);
        tourism.push_back(new Tourism(source, destination, dist,label_org,label_dst));
    }

    int num_nodes = nodes.size();
    std::vector<std::vector<double>> matrix(num_nodes, std::vector<double>(num_nodes, 0.0));

    for (auto s : tourism) {
        matrix[s->getSrc()][s->getDst()] = s->getDist();
    }

    // Clean up dynamically allocated memory
    for (auto s : tourism) {
        delete s;
    }

    return matrix;
}












