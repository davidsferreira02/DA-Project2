//
// Created by David Ferreira on 04/05/2024.
//

#ifndef PROJ1_TOURISM_H
#define PROJ1_TOURISM_H

#include <string>

class Tourism {


protected:
    int src;
    int dst;
    int dist;
    std::string label_org;
    std::string label_dst;


public:
    Tourism(int src, int dst , int dist,std::string label_org,std::string label_dst);
    virtual ~Tourism();
    int getSrc();
    int getDst() ;
    int getDist() ;
    std::string getLabel_org();
    std::string getLabel_dst();


};


#endif //PROJ1_TOURISM_H
