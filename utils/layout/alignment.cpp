#include "alignment.hpp"

alignment::~alignment(void){}

void alignment::revComp(void){

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

char alignment::getStrand(void){
    return strand;
}

