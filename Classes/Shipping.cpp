//
// Created by David Ferreira on 04/05/2024.
//

#include "Shipping.h"



Shipping::Shipping(int src,int dst,double dist) : src(src),dst(dst),dist(dist){}

Shipping::~Shipping() = default;

int Shipping::getSrc() {
    return src;
}

double Shipping::getDist() {
    return dist;
}

int Shipping::getDst() {
    return dst;
}