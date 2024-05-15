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
/*class Reader {
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
    std::vector<double> distance;
    int num_nodes;


public:
    Reader();cpp:443: undefined reference to `Reader::kMeansClustering(Graph<int> const&, int, std::unordered_map<int, Reader::Coordinates, std::hash<int>, std::equal_to<int>, std::allocator<std::pair<int const, Reader::Coordinates> > > const&)'
/usr/bin/ld: CMakeFiles/YourExecutable.dir/reader.cpp.o: in function `Reader::readAndParse4_2Extra_Fully_Connected_Graphs(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)':
/home/rebelojoao/Documents/DA/prj2/reader.cpp:147: undefined

    void readAndParseStadium();

    void readAndParseShipping();

    void readAndParseTourism();


    int* getNode(int src);

    Graph<int*> getGraph(){
        return graph;
    }
    int getNumNodes() const { return num_nodes; }

};
*/
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