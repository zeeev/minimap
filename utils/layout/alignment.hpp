#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>
#include <iostream>

class alignment{

    typedef long int LI ;

public:

    std::string qName  ;
    LI   qLen   ;
    LI   qStart ;
    LI   qEnd   ;

    std::string tName  ;
    LI   tLen   ;
    LI   tStart ;
    LI   tEnd   ;

    char strand ;

    std::string line   ;

    LI match;

    bool flipped;

    alignment( std::string a,
               LI b, LI c, LI d,
               std::string e,
               LI f, LI g, LI h,
               char i,
               LI j  ,
               bool k,
               std::string l) :
        qName(a),
        qStart(b),
        qEnd(c),
        qLen(d),
        tName(e),
        tStart(f),
        tEnd(g),
        tLen(h),
        strand(i),
        match(j),
        flipped(k),
        line(l)
    {}

    ~alignment(void);

    void revComp(void);

    char getStrand(void);

    inline friend std::ostream & operator<< (std::ostream & os, const alignment & al);

};


std::ostream & operator<< (std::ostream & os, const alignment & al)
{
    os << al.qName
       << "\t" << al.qStart
       << "\t" << al.qEnd
       << "\t" << al.tName
       << "\t" << al.tStart
       << "\t" << al.tEnd
       << "\t" << al.strand
       << "\t" <<  al.line;

    return os;
}

#endif /*  */
