#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <bitset>
#include <algorithm>
#include "split.hpp"
#include "alignment.hpp"
#include "chain.hpp"
#include "interval_tree.h"

typedef std::vector<chain *> chains;

bool revSortByScore(chain * L, chain * R){
    return L->getMatchingBases() > R->getMatchingBases();
}

bool processSeqid(std::vector<std::string> & lines,
                  chains & mc){

    if(lines.size() == 0){
        return true;
    }

    long int lastQend  =  -1 ;
    long int lastTend  =  -1 ;
    char lastStrand    = '-' ;

    std::vector<alignment *> als;

    std::vector<std::string>::iterator it = lines.begin();

    for( ; it != lines.end(); it++){

        alignment * al = new alignment;

        if(!al->parseMMLine(*it)){
            std::cerr << "FATAL: did not parse" << std::endl;
            return false;
        }

        if(al->strand == '-'){
            al->revCompT();
        }

        long int qGap = al->qStart - lastQend;
        long int tGap = al->tStart - lastTend;

        if(( al->tStart < lastTend  ||
            al->qStart < lastQend )  && als.size() > 0){
                chain * lc = new chain;
                lc->addAlignment(als);
                mc.push_back(lc);
                als.clear();
                als.push_back(al);
                lastTend = -1;
                lastQend = -1;
                lastStrand = al->strand;
        }
        else{
            //            std::cerr << qGap << "\t" << tGap << std::endl;
            als.push_back(al);
            lastTend = al->tEnd;
            lastQend = al->qEnd;
            lastStrand = al->strand;
        }
    }
    if(als.size() > 0){
        chain * lc = new chain;
        lc->addAlignment(als);
        mc.push_back(lc);
    }
    return true;
}


bool parseFile(chains & myChains){
    std::string line ;
    std::string tName;
    std::string qName;
    std::vector<std::string> lines;

    long int nAlignments = 0;

    while(getline(std::cin, line)){

        nAlignments += 1;
        std::vector<std::string> lineDat = split(line, "\t");

        if(lineDat[0] != tName || lineDat[5] != qName ){
            if(!processSeqid(lines, myChains)){
                return false;
            }
            lines.clear();
            lines.push_back(line);
            tName = lineDat[0];
            qName = lineDat[5];
        }
        else{
            lines.push_back(line);
        }
    }
    if(!processSeqid(lines, myChains)){
        return false;
    }
    std::cerr << "INFO: there are " << nAlignments << " alignments" << std::endl;
    return true;
}

int main(int argc, char ** argv)
{

    chains myChains;

    if(!parseFile(myChains)){
        std::cerr << "FATAL: could not parse file.";
        return 1;
    }

    std::cerr << "INFO: parsed all aligments." << std::endl;
    std::cerr << "INFO: there are " << myChains.size() << " chains" << std::endl;

    sort(myChains.begin(), myChains.end(), revSortByScore);

    IntervalTree< chain * >  intervals;

    for(chains::iterator it = myChains.begin();
        it != myChains.end(); it++){

        (*it)->chainToBed();

        std::vector< chain *> results = intervals.fetch(
                                                        (*it)->getTMin(),
                                                        (*it)->getTMax());
        bool hit = false;

        for(chains::iterator iz = results.begin();
            iz != results.end(); iz++){
            if((*iz)->getTName() == (*it)->getTName()){
                hit = true;
            }
        }
        if(hit){
        }
        else{
            (*it)->printBed();
            intervals.insert((*it), (*it)->getTMin(),
                             (*it)->getTMax());
        }
        std::cerr << **it << std::endl;
    }
    return 0;
}
