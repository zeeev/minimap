#include "chain.hpp"

void chain::_printBed(void){

    char strand = this->alignments.front()->strand;



    std::cout << this->tName
              << "\t"
              << this->tMin
              << "\t"
              << this->tMax
              << "\t"
              << this->qName
              << "\t"
              << this->qMin
              << "\t"
              << this->qMax
              << "\t"
              << strand
              << "\t"
              << this->matchingBases
              << "\t"
              << this->alignments.size()
              << "\t"
              << this->tLen
              << "\t"
              << this->qLen
              << std::endl;
}


long int chain::getTMin(void){
    return this->tMin;
}
long int chain::getTMax(void){
    return this->tMax;
}
long int chain::getQMin(void){
    return this->qMin;
}
long int chain::getQMax(void){
    return this->qMax;
}

std::string chain::getTName(void){
    return this->tName;
}

std::string chain::getQName(void){
    return this->qName;
}

int chain::getNAlignments(void){
    return alignments.size();
}

chain::chain(void){
    matchingBases = 0;
    qLen  = 0;
    tLen  = 0;
    tMin  = 0;
    tMax  = 0;
    qMin  = 0;
    qMax  = 0;

    twoStrands = false;
}

void chain::cleanUp(void){
    for(std::vector<alignment *>::iterator it
            = alignments.begin(); it != alignments.end(); it++){
        delete *it;
    }
}

chain::~chain(void){
}

bool chain::addAlignment(std::vector<alignment *> & als){
    for(std::vector<alignment *>::iterator it = als.begin();
        it != als.end(); it++){

        if(qName.empty()){
            alignments.push_back(*it);
            qName = (*it)->qName;
            tName = (*it)->tName;
            tLen  = (*it)->tLen;
            qLen  = (*it)->qLen;
        }
        else{
            if((*it)->qName != qName){
                std::cerr << "WANRING: name not the same" << std::endl;
                return false;
            }
            alignments.push_back(*it);
        }
        this->matchingBases += (*it)->match;
    }
    return true;
}

bool chain::addAlignment(alignment * al){

    if(qName.empty()){
        qName = al->qName;
        tName = al->tName;
    }
    else{
        if(al->qName != qName){
            std::cerr << "WANRING: name not the same" << std::endl;
            return false;
        }
    }

    this->tLen = al->tLen;
    this->qLen = al->qLen;

    this->matchingBases += al->match;

    alignments.push_back(al);

    return true;
}

long int chain::getMatchingBases(void){
    return this->matchingBases;
}

void chain::printBed(void){

    bool flipped = this->alignments.front()->flipped;

    if(flipped == true){
        flipped = false;
    }
    else{
        flipped = true;
    }

    std::vector<chain *> chains;

    for(std::vector<alignment *>::iterator it = this->alignments.begin();
        it != this->alignments.end(); it++){

        if((*it)->flipped != flipped){
             chain * nc = new chain;
            nc->addAlignment(*it);
            chains.push_back(nc);
            flipped = (*it)->flipped;
        }
        else{
            chains.back()->addAlignment(*it);
        }
    }

    for(std::vector<chain *>::iterator it = chains.begin();
        it != chains.end(); it++){

        (*it)->chainToBed();
        (*it)->_printBed();

        delete *it;
    }
}

bool chain::chainToBed(void){

    this->tMin = this->alignments.front()->tStart;
    this->tMax = this->alignments.front()->tEnd;
    this->qMin = this->alignments.front()->qStart;
    this->qMax = this->alignments.front()->qEnd;

    for(std::vector<alignment *>::iterator it = this->alignments.begin();
        it != this->alignments.end(); it++){

        if((*it)->tStart < tMin ){
            tMin = (*it)->tStart;
        }
        if((*it)->qStart < qMin ){
            qMin = (*it)->qStart;
        }
        if((*it)->tEnd > tMax){
            tMax = (*it)->tEnd;
        }
        if((*it)->qEnd > qMax){
            qMax = (*it)->qEnd;
        }
    }
    return true;
}
