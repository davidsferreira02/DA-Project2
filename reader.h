#ifndef READER_H
#define READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Classes/Graph.h"
#include "Classes/Stadium.h"
#include "Classes/Tourism.h"
#include "Classes/Shipping.h"


/**
 * @brief The Reader class is responsible for reading and parsing input data
 *        to construct the graph representing the water distribution network.
 */
class Reader {
private:
    Graph<int*> graph;  ///< The graph representing the water distribution network.
    std::unordered_map<int, Stadium*> srcMapStadium;
    std::unordered_map<int, Stadium*> dstMapStadium;
    std::unordered_map<int, Stadium*> distMapStadium;
    std::unordered_map<int, Shipping*> srcMapShipping;
    std::unordered_map<int, Shipping*> dstMapShipping;
    std::unordered_map<int, Shipping*> distMapShipping;
    std::unordered_map<int, Tourism*> srcMapTourism;
    std::unordered_map<int, Tourism*> dstMapTourism;
    std::unordered_map<int, Tourism*> distMapTourism;
    std::unordered_map<std::string,Tourism*> labelOrgTourism;
    std::unordered_map<std::string, Tourism*> labelDstTourism;

public:
    Reader();

    void readAndParseStadium();

    void readAndParseShipping();

    void readAndParseTourism();


    int* getNode(int src);

    Graph<int*> getGraph(){
        return graph;
    }
};

#endif /* READER_H */