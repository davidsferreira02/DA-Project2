//
// Created by David Ferreira on 04/05/2024.
//

#ifndef PROJ1_STADIUM_H
#define PROJ1_STADIUM_H


class Stadium {

protected:
    int src;
    int dst;
    double dist;

public:
    Stadium(int src, int dst , double dist);
    virtual ~Stadium();
    int getSrc();
    int getDst() ;
    double getDist() ;
};


#endif //PROJ1_STADIUM_H
