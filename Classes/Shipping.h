//
// Created by David Ferreira on 04/05/2024.
//

#ifndef PROJ1_SHIPPING_H
#define PROJ1_SHIPPING_H


class Shipping {



    protected:
        int src;
        int dst;
        double dist;

    public:
        Shipping(int src, int dst , double dist);
        virtual ~Shipping();
        int getSrc();
        int getDst() ;
        double getDist() ;
    };



#endif //PROJ1_SHIPPING_H
