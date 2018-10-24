#ifndef RNA_H
#define RNA_H

#include <string>

class RNA {
    public:
        RNA();
        RNA(std::string);
        ~RNA();
        uint            length() const;
        std::string     str() const;
        // char            seq(uint k) const;
        // bool            isConsensus() const;
        // bool            containsPseudoBases() const;
    private:
        std::string     _seq;
        // bool            _coding;
        // bool            _consensus;
        // bool            _containsPseudoNTs;
};

inline uint RNA::length() const { return _seq.length(); }
inline std::string RNA::str() const { return _seq; }

#endif