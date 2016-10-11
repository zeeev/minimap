#include "alignment.hpp"

bool alignment::parseMMLine(std::string & line){

    std::vector<std::string> ld = split(line, "\t");

    if(ld.size() < 10){
        return false;
    }

    this->qName   = ld[0];
    this->qStart  = atol(ld[2].c_str());
    this->qEnd    = atol(ld[3].c_str());
    this->qLen    = atol(ld[1].c_str());
    this->tName   = ld[5];
    this->tStart  = atol(ld[7].c_str());
    this->tEnd    = atol(ld[8].c_str());
    this->tLen    = atol(ld[6].c_str());
    this->strand  = ld[4][0];
    this->match   = atol(ld[9].c_str());
    this->flipped = false;
    this->line    = line;

    return true;

}



alignment::alignment(void){
}

alignment::~alignment(void){
}

void alignment::revCompQ(void){

   LI tmp = qStart       ;
   qStart = qLen - qEnd  ;
   qEnd   = qLen - tmp   ;

   if(strand == 45){
       strand = 43;
   }
   else{
       strand = 45;
   }
   if(flipped == true){
       flipped = false;
   }
   else{
       flipped = true;
   }
 }


void alignment::revCompT(void){

    LI tmp = tStart       ;
    tStart = tLen - tEnd  ;
    tEnd   = tLen - tmp   ;

    if(flipped == true){
        flipped = false;
    }
    else{
        flipped = true;
    }
}


char alignment::getStrand(void){
    return strand;
}

