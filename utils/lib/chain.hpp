#ifndef CHAIN_H
#define CHAIN_H

#include "alignment.hpp"
#include <vector>
#include <map>

class chain{
    typedef long int LI;

private:

    LI matchingBases;
    LI tLen         ;
    LI qLen         ;
    LI tMin         ;
    LI tMax         ;
    LI qMin         ;
    LI qMax         ;

    std::string tName;
    std::string qName;

    bool twoStrands;

    void _printBed(void);

public:

    std::vector<alignment *> alignments;

    chain(void);
    ~chain(void);

    bool addAlignment(alignment * );
    bool addAlignment(std::vector<alignment *> &);
    bool chainToBed(void);

    std::string getTName(void);
    std::string getQName(void);

    LI getMatchingBases( void     );
    int getNAlignments( void    );

    LI getQMax(void);
    LI getQMin(void);
    LI getTMax(void);
    LI getTMin(void);

    void printBed(void);
    void cleanUp(void);

    inline friend std::ostream & operator<< (std::ostream & os, chain & ch);

};




std::ostream & operator<< (std::ostream & os, chain & ch)
{
    os << "chain score: " << ch.matchingBases  << std::endl;
    os << "n chain:  " << ch.alignments.size() << std::endl;
    for(std::vector<alignment *>::iterator it = ch.alignments.begin();
        it != ch.alignments.end(); it++){

        if((*it)->flipped){
            (*it)->revCompT();
        }

        os << (**it) << std::endl;
    }

    return os;
}

#endif /*  */
