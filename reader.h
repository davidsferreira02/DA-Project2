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

    Graph<int> readAndParseExtra_Fully_Connected_Graphs(const std::string filename);
    Graph<int> readAndParse4_2Extra_Fully_Connected_Graphs(const std::string filename);
    Graph<int> readAndParseRealWorld_Graphs(int graphNumber,  std::unordered_map<int, Vertex<int>*> &vertexMap);

    struct Coordinates {
        double latitude;
        double longitude;
    };
    std::unordered_map<int, Coordinates> readCoordinates();
    std::vector<std::vector<int>> kMeansClustering(const Graph<int>& graph, int k, const std::unordered_map<int, Coordinates>& coordinates);
};

#endif /* READER_H */