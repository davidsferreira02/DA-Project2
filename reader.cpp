#include "reader.h"
#include "Classes/Shipping.h"


/**
 * @brief Read and parse data from the Stations.csv file to create Station objects.
 *
 * This function reads data from the Stations.csv file, creates Station objects, and adds them to the graph.
 */
void Reader::readAndParseStadium() {
    std::ifstream file("../Data/Toy-Graphs/stadiums.csv");
    std::string line;

    getline(file, line); // Skip the header line

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::string src_str, dst_str, dist_str;

        getline(ss, src_str, ',');
        int src = std::stoi(src_str);

        getline(ss, dst_str, ',');
        int dst = std::stoi(dst_str);

        getline(ss, dist_str);
        double dist = std::stod(dist_str);

        graph.addVertex(&src);
        graph.addVertex(&dst);

        graph.addEdge(&src, &dst, dist);
        auto *DS= new Stadium(src,dst,dist);
        srcMapStadium[src]=DS;
        dstMapStadium[dst]=DS;
    }

    file.close();
}



void Reader::readAndParseTourism() {
    std::ifstream file("../Data/Toy-Graphs/stadiums.csv");
    std::string line;

    getline(file, line); // Skip the header line

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::string src_str, dst_str, dist_str,label_org,label_dst;

        getline(ss, src_str, ',');
        int src = std::stoi(src_str);

        getline(ss, dst_str, ',');
        int dst = std::stoi(dst_str);

        getline(ss, dist_str);
        double dist = std::stod(dist_str);

        getline(ss, label_org);
        getline(ss,label_dst);

        graph.addVertex(&src);
        graph.addVertex(&dst);

        graph.addEdge(&src, &dst, dist);
        auto *DS= new Tourism(src,dst,dist,label_org,label_dst);
        labelOrgTourism[label_org]=DS;
        labelDstTourism[label_dst]=DS;
    }

    file.close();
}


void Reader::readAndParseShipping() {
    std::ifstream file("../Data/Toy-Graphs/shipping.csv");
    std::string line;

    getline(file, line); // Skip the header line

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::string src_str, dst_str, dist_str;

        getline(ss, src_str, ',');
        int src = std::stoi(src_str);

        getline(ss, dst_str, ',');
        int dst = std::stoi(dst_str);

        getline(ss, dist_str);
        double dist = std::stod(dist_str);

        graph.addVertex(&src);
        graph.addVertex(&dst);

        graph.addEdge(&src, &dst, dist);
        auto *DS= new Shipping(src,dst,dist);
        srcMapShipping[src]=DS;
        dstMapShipping[dst]=DS;
    }

    file.close();
}









Reader::Reader() = default;