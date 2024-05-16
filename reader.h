#ifndef READER_H
#define READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Classes/Graph.h"

class Reader {
private:


public:
    Graph<int> readAndParseStadium();
    Graph<int>readAndParseShipping();
    Graph<int> readAndParseTourism();

    Graph<int> readAndParse4_2Extra_Fully_Connected_Graphs(const std::string filename,std::unordered_map<int, Vertex<int>*> &vertexMap, std::unordered_map<std::string, Edge<int>*> &edgeMap);
    Graph<int> readAndParseRealWorld_Graphs(int graphNumber, std::unordered_map<int, Vertex<int>*> &vertexMap, std::unordered_map<std::string, Edge<int>*> &edgeMap);

    struct Coordinates {
        double latitude;
        double longitude;
    };
    std::unordered_map<int, Coordinates> readCoordinates();
    std::vector<std::vector<int>> kMeansClustering(const Graph<int>& graph, int k, const std::unordered_map<int, Coordinates>& coordinates);
};

#endif /* READER_H */