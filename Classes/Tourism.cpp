//
// Created by David Ferreira on 04/05/2024.
//

#include "Tourism.h"


Tourism::Tourism(int src, int dst, int dist, std::string label_org, std::string label_dst) :src(src),dst(dst),dist(dist),label_org(label_org),label_dst(label_dst){}


Tourism::~Tourism() = default;


int Tourism::getDst() {
    return dst;
}

int Tourism::getSrc() {
    return src;
}

int Tourism::getDist() {
    return dist;
}

std::string Tourism::getLabel_dst() {
    return label_dst;
}

std::string Tourism::getLabel_org() {
    return label_org;
}

