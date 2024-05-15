#include "reader.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <complex>

std::unordered_map<int, Reader::Coordinates> Reader::readCoordinates() {
    std::ifstream file("../Data/Extra_Fully_Connected_Graphs/nodes.csv");
    std::unordered_map<int, Coordinates> nodeCoordinates;
    std::string line;
    if (file.is_open()) {
        std::getline(file, line);
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

double convert_to_radians(double coord) {
    return coord * 3.14 / 180.0;
}

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

    const double earthRadius = 6371000;

    return earthRadius * c;
}

std::vector<std::vector<int>> Reader::kMeansClustering(const Graph<int>& graph, int k, const std::unordered_map<int, Coordinates>& coordinates) {
    std::vector<std::vector<int>> clusters(k);
    std::vector<int> vertexIds;

    for (const auto& vertex : graph.getVertexSet()) {
        vertexIds.push_back(vertex->getInfo());
    }

    std::vector<int> centroids;
    centroids.push_back(vertexIds[0]);
    std::unordered_set<int> visitedCentroids;
    visitedCentroids.insert(vertexIds[0]);

    while (centroids.size() < k) {
        double maxDist = -1.0;
        int farthestVertex = -1;
        for (int vertex : vertexIds) {
            if (visitedCentroids.find(vertex) == visitedCentroids.end()) {
                double minDist = INF;
                for (int centroid : centroids) {
                    double dist = Haversine(
                            coordinates.at(vertex).latitude, coordinates.at(vertex).longitude,
                            coordinates.at(centroid).latitude, coordinates.at(centroid).longitude
                    );
                    if (dist < minDist) {
                        minDist = dist;
                    }
                }
                if (minDist > maxDist) {
                    maxDist = minDist;
                    farthestVertex = vertex;
                }
            }
        }
        if (farthestVertex != -1) {
            centroids.push_back(farthestVertex);
            visitedCentroids.insert(farthestVertex);
        } else {
            break;
        }
    }

    bool changed = true;
    while (changed) {
        for (auto& cluster : clusters) {
            cluster.clear();
        }

        for (const auto& vertex : graph.getVertexSet()) {
            int nearestCentroid = -1;
            double minDistance = INF;
            for (int i = 0; i < k; ++i) {
                double dist = Haversine(
                        coordinates.at(vertex->getInfo()).latitude, coordinates.at(vertex->getInfo()).longitude,
                        coordinates.at(centroids[i]).latitude, coordinates.at(centroids[i]).longitude
                );
                if (dist < minDistance) {
                    minDistance = dist;
                    nearestCentroid = i;
                }
            }
            clusters[nearestCentroid].push_back(vertex->getInfo());
        }

        changed = false;
        for (int i = 0; i < k; ++i) {
            double avgLat = 0.0, avgLon = 0.0;
            for (int vertex : clusters[i]) {
                avgLat += coordinates.at(vertex).latitude;
                avgLon += coordinates.at(vertex).longitude;
            }
            avgLat /= clusters[i].size();
            avgLon /= clusters[i].size();

            double minDistance = INF;
            int newCentroid = centroids[i];
            for (int vertex : clusters[i]) {
                double dist = Haversine(avgLat, avgLon, coordinates.at(vertex).latitude, coordinates.at(vertex).longitude);
                if (dist < minDistance) {
                    minDistance = dist;
                    newCentroid = vertex;
                }
            }

            if (newCentroid != centroids[i]) {
                centroids[i] = newCentroid;
                changed = true;
            }
        }
    }

    return clusters;
}

Graph<int> Reader::readAndParse4_2Extra_Fully_Connected_Graphs(const std::string filename) {
    std::ifstream file(filename);
    std::string line;

    Graph<int> graph;
    std::unordered_map<int, Coordinates> coordinates = readCoordinates();

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
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

    Graph<int> graph;

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');

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

    file.close();
    return graph;



}


Graph<int> Reader::readAndParseExtra_Fully_Connected_Graphs(const std::string filename) {
    std::ifstream file("../Data/Extra_Fully_Connected_Graphs/edges_25.csv");
    std::string line;

    Graph<int> graph;

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');

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

    file.close();
    return graph;



}

Graph<int> Reader::readAndParseShipping() {
    std::ifstream file("../Data/Toy-Graphs/shipping.csv");
    std::string line;

    Graph<int> graph;

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');

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

    file.close();
    return graph;
}


Graph<int> Reader::readAndParseTourism() {
    std::ifstream file("../Data/Toy-Graphs/tourism.csv");
    std::string line;

    Graph<int> graph;

    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return graph;
    }

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');

        std::istringstream iss(line);
        int source, dest;
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

Graph<int> Reader::readAndParseRealWorld_Graphs(int graphNumber, std::unordered_map<int, Vertex<int>*> &vertexMap)
{
    Graph<int> graph;
    std::string line;

    std::string filePath;
    switch (graphNumber) {
        case 1:
            filePath = "../Data/Real-world Graphs/graph1/edges.csv";
            break;
        case 2:
            filePath = "../Data/Real-world Graphs/graph2/edges.csv";
            break;
        case 3:
            filePath = "../Data/Real-world Graphs/graph3/edges.csv";
            break;
        default:
            std::cerr << "Invalid graph number!" << std::endl;
            return graph;
    }

    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return graph;
    }

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');

        std::istringstream iss(line);
        int source, dest;
        double dist;

        if (!(iss >> source >> dest >> dist)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }
        graph.addEdge(source, dest, dist);
        graph.addEdge(dest, source, dist);
        if (vertexMap.find(source) == vertexMap.end()) {
            vertexMap[source] = new Vertex<int>(static_cast<int>(source));
            graph.addVertex(source);

        }
        if (vertexMap.find(dest) == vertexMap.end()) {
            vertexMap[dest] = new Vertex<int>(static_cast<int>(dest));
            graph.addVertex(dest);
        }
    }

    return graph;
}













