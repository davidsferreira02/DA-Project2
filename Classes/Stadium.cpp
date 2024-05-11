//
// Created by David Ferreira on 04/05/2024.
//

#include "Stadium.h"



Stadium::Stadium(int src,int dst,double dist) : src(src),dst(dst),dist(dist){}


Stadium::~Stadium() = default;


int Stadium::getSrc() {
    return src;
}

double Stadium::getDist() {
    return dist;
}

int Stadium::getDst() {
    return dst;
}