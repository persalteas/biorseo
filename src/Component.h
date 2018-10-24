#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>

class Component 
{
    public:
        Component();
        Component(std::string& cons_seq, uint k);
    private:
        std::string     _cons_seq;
        uint            k;
};


istream& operator>>(istream& is, Component& m);
ostream& operator<<(ostream& os, const Component& m);

#endif