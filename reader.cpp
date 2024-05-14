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
#include <complex>


struct Coordinates {
    double latitude;
    double longitude;
};

std::unordered_map<int, Coordinates> readCoordinates() {
    std::ifstream file("../Data/Extra_Fully_Connected_Graphs/nodes.csv");
    std::unordered_map<int, Coordinates> nodeCoordinates;
    std::string line;
    if (file.is_open()) {
        std::getline(file, line); // Skip header
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            std::getline(iss, token, ',');
            int id = std::stoi(token);
            std::getline(iss, token, ',');
            double longitude = std::stod(token);
            std::getline(iss, token, ',');
            double latitude = std::stod(token);
            nodeCoordinates[id] = {latitude, longitude};
        }
        file.close();
    }
    return nodeCoordinates;
}

// Function to convert degrees to radians
double convert_to_radians(double coord) {
    return coord * 3.14 / 180.0;
}

// Haversine function to calculate the distance between two points on Earth
double Haversine(double lat1, double lon1, double lat2, double lon2) {

    double rad_lat1 = convert_to_radians(lat1);
    double rad_lon1 = convert_to_radians(lon1);
    double rad_lat2 = convert_to_radians(lat2);
    double rad_lon2 = convert_to_radians(lon2);


    double delta_lat = rad_lat2 - rad_lat1;
    double delta_lon = rad_lon2 - rad_lon1;

    double a = std::pow(std::sin(delta_lat / 2), 2) +
               std::cos(rad_lat1) * std::cos(rad_lat2) * std::pow(std::sin(delta_lon / 2), 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    // Earth's radius in meters
    const double earthRadius = 6371000;

    // Calculate and return the distance
    return earthRadius * c;
}





Graph<int> Reader::readAndParse4_2Extra_Fully_Connected_Graphs(const std::string filename) {
    std::ifstream file("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
    std::string line;

    Graph<int> graph;
    std::unordered_map<int, Coordinates> coordinates = readCoordinates();

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' '); // Replace commas with spaces for easier parsing

        std::istringstream iss(line);
        int source, dest;
        double dist;

        if (!(iss >> source >> dest >> dist)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }
        graph.addVertex(source);
        graph.addVertex(dest);
        graph.addEdge(source, dest, dist);
        graph.addEdge(dest, source, dist);
    }

    auto vertices = graph.getVertexSet();
    for (auto& v : vertices) {
        for (auto& w : vertices) {
            if (v->getInfo() != w->getInfo() && !graph.findEdge(v->getInfo(), w->getInfo())) {
                double dist = Haversine(
                        coordinates[v->getInfo()].latitude, coordinates[v->getInfo()].longitude,
                        coordinates[w->getInfo()].latitude, coordinates[w->getInfo()].longitude
                );
                graph.addEdge(v->getInfo(), w->getInfo(), dist);
                graph.addEdge(w->getInfo(), v->getInfo(), dist);
            }
        }
    }

    file.close();
    return graph;
}



Graph<int> Reader::readAndParseStadium() {
    std::ifstream file("../Data/Toy-Graphs/stadiums.csv");
    std::string line;

    // Assuming T is int and the Graph is for an undirected graph with distances as weights
    Graph<int> graph;

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' '); // Replace commas with spaces for easier parsing

        std::istringstream iss(line);
        int source, dest;
        double dist;

        if (!(iss >> source >> dest >> dist)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }

        // Add vertices and edge
        graph.addVertex(source);
        graph.addVertex(dest);
        graph.addEdge(source, dest, dist);
        graph.addEdge(dest, source, dist);
    }

    file.close();
    return graph;



}


Graph<int> Reader::readAndParseExtra_Fully_Connected_Graphs(const std::string filename) {
    std::ifstream file("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
    std::string line;

    // Assuming T is int and the Graph is for an undirected graph with distances as weights
    Graph<int> graph;

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' '); // Replace commas with spaces for easier parsing

        std::istringstream iss(line);
        int source, dest;
        double dist;

        if (!(iss >> source >> dest >> dist)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }

        // Add vertices and edge
        graph.addVertex(source);
        graph.addVertex(dest);
        graph.addEdge(source, dest, dist);
        graph.addEdge(dest, source, dist);
    }

    file.close();
    return graph;



}

Graph<int> Reader::readAndParseShipping() {
    std::ifstream file("../Data/Toy-Graphs/shipping.csv");
    std::string line;

    // Assuming T is int and the Graph is for an undirected graph with distances as weights
    Graph<int> graph;

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' '); // Replace commas with spaces for easier parsing

        std::istringstream iss(line);
        int source, dest;
        double dist;

        if (!(iss >> source >> dest >> dist)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }

        // Add vertices and edge
        graph.addVertex(source);
        graph.addVertex(dest);
        graph.addEdge(source, dest, dist);
        graph.addEdge(dest, source, dist);
    }

    file.close();
    return graph;
}


Graph<int> Reader::readAndParseTourism() {
    std::ifstream file("../Data/Toy-Graphs/tourism.csv");
    std::string line;

    // Assuming T is int and the Graph is for an undirected graph with distances as weights
    Graph<int> graph;

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' '); // Replace commas with spaces for easier parsing

        std::istringstream iss(line);
        int source, dest, label_source, label_dest;
        double dist;
        std::string label_source_str, label_dest_str;

        if (!(iss >> source >> dest >> dist >> label_source_str >> label_dest_str)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }

        graph.addVertex(source);
        graph.addVertex(dest);
        graph.addEdge(source, dest, dist);
        graph.addEdge(dest, source, dist);
    }

    file.close();
    return graph;
}


















