#ifndef MOTIF_H
#define MOTIF_H

#include <vector>
#include <string>
#include "Component.h"

class Motif
{
    public:
        Motif();
        Motif(std::string filename);
        ~Motif();
        std::vector<Component>& getComponents() const;
        Component& getComponent(uint k) const;
        std::string name;

    private:
        std::vector<Component>  _comps;

};

istream& operator>>(istream& is, Motif& m);
ostream& operator<<(ostream& os, const Motif& m);

#endif