#include <stdlib.h>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "alignment.hpp"
#include "split.hpp"

struct matchInfo{
    long int posStrand;
    long int negStrand;
};

void addMatch(matchInfo * mi, char strand, long int match)
{
    if( strand == '-' ){
        mi->negStrand += match ;
    }
    else{
        mi->posStrand += match ;
    }
}

typedef std::map<std::string, std::map< std::string, matchInfo * > > qtCount;

bool sortT(const alignment * al1, const alignment * al2 ){
    if(al1->tName == al2->tName){
        return al1->tStart < al2->tEnd;
    }
    return al1->tName < al2->tName;
}

bool sortQ(const alignment * al1, const alignment * al2 ){
    if(al1->qName == al2->qName){
        return al1->qStart < al2->qEnd;
    }
    return al1->qName < al2->qName;
}

int main(int argc, char ** argv)
{
    std::vector<std::string> tOrder ;
    tOrder.push_back("chr1");
    tOrder.push_back("chr2");
    tOrder.push_back("chr3");
    tOrder.push_back("chr4");
    tOrder.push_back("chr5");
    tOrder.push_back("chr6");
    tOrder.push_back("chr7");
    tOrder.push_back("chr8");
    tOrder.push_back("chr9");
    tOrder.push_back("chr10");
    tOrder.push_back("chr11");
    tOrder.push_back("chr12");
    tOrder.push_back("chr13");
    tOrder.push_back("chr14");
    tOrder.push_back("chr15");
    tOrder.push_back("chr16");
    tOrder.push_back("chr17");
    tOrder.push_back("chr18");
    tOrder.push_back("chr19");
    tOrder.push_back("chr20");
    tOrder.push_back("chr21");
    tOrder.push_back("chr22");
    tOrder.push_back("chrX");
    tOrder.push_back("chrY");

    std::vector<alignment *> records_tSorted;
    std::vector<alignment *> records_qSorted;

    qtCount qtMatch;

    std::map<std::string, std::string> bestTarget;
    std::map<std::string, long int>         tLens;

    std::string line;

    while(getline(std::cin, line)){

        std::vector<std::string> ld = split(line, "\t");

        alignment * al = new alignment (ld[0],
                                        atol(ld[2].c_str()),
                                        atol(ld[3].c_str()),
                                        atol(ld[1].c_str()),
                                        ld[5],
                                        atol(ld[7].c_str()),
                                        atol(ld[8].c_str()),
                                        atol(ld[6].c_str()),
                                        ld[4][0],
                                        atol(ld[9].c_str()),
                                        false,
                                        line);

        tLens[ld[5]] = atol(ld[6].c_str());

        if(qtMatch.find(al->qName)
           == qtMatch.end()){

            matchInfo * mi = new matchInfo;
            mi->negStrand = 0;
            mi->posStrand = 0;

            addMatch(mi, al->strand, al->match);
            qtMatch[al->qName][al->tName] = mi;
        }

        else if(qtMatch[al->qName].find(al->tName)
                == qtMatch[al->qName].end() ){

            matchInfo * mi = new matchInfo;
            mi->negStrand = 0;
            mi->posStrand = 0;

            addMatch(mi, al->strand, al->match);
            qtMatch[al->qName][al->tName] = mi;

        }
        else{
            matchInfo * mi = qtMatch[al->qName][al->tName];


            addMatch(mi, al->strand, al->match);
        }

        if(al->strand == '-'){
            al->revComp();
        }

        records_tSorted.push_back(al);
        records_qSorted.push_back(al);

    }

    std::cerr << "INFO: loaded "  << records_tSorted.size()
              << " alignments " << std::endl;
    std::cerr << "INFO: Sorting alignment blocks." << std::endl;

    std::sort(records_tSorted.begin(), records_tSorted.end(), sortT);
    std::sort(records_qSorted.begin(), records_qSorted.end(), sortQ);

    for(qtCount::iterator it = qtMatch.begin();
        it != qtMatch.end(); it++){

        std::string best    ;
        long int highest = 0;

        for(std::map<std::string, matchInfo * >::iterator iz
                = it->second.begin();
            iz != it->second.end(); iz++){

            if(iz->second->negStrand +
               iz->second->posStrand  > highest){
                best   = iz->first;
                highest = iz->second->negStrand + iz->second->posStrand;
            }
        }
        bestTarget[it->first] = best;
    }

    std::map<std::string, long int > qOffset;
    std::map<std::string, long int > tOffset;

    std::map<std::string, long int > highQmatch;

    long int qSum = 0;


    for(std::vector<std::string>::iterator it = tOrder.begin();
        it != tOrder.end(); it++){
        for(std::vector< alignment * >::iterator i = records_tSorted.begin();
            i != records_tSorted.end(); i++){

            if((*it) != (*i)->tName){
                continue;
            }

            if(bestTarget[(*i)->qName] == (*i)->tName){
                if(highQmatch.find((*i)->qName) == highQmatch.end()){
                    qOffset[(*i)->qName] = qSum;
                    highQmatch[(*i)->qName] = (*i)->match;
                }
                else{
                    if(highQmatch[(*i)->qName] < (*i)->match){
                        qOffset[(*i)->qName] = qSum;
                        highQmatch[(*i)->qName] = (*i)->match;
                    }
                }
            }
            qSum += abs((*i)->qEnd - (*i)->qStart);
        }
    }
    long int tSum = 0;

    for(std::vector<std::string>::iterator it = tOrder.begin();
        it != tOrder.end(); it++){
        if(tLens.find(*it) == tLens.end()){
            std::cerr << "FATAL: target length not found " << std::endl;
            return 1;
        }
        tOffset[*it] = tSum       ;
        tSum         += tLens[*it];
        tLens.erase(*it)          ;
    }

    for(std::map<std::string, long int>::iterator i = tLens.begin();
        i != tLens.end(); i++){
        std::cerr << i->first << "\t" << tSum << std::endl;
        tOffset[i->first] = tSum;
        tSum += i->second;
    }

    for(std::vector< alignment * >::iterator i = records_tSorted.begin();
        i != records_tSorted.end(); i++){


        if(qtMatch[(*i)->qName][(*i)->tName]->negStrand
           > qtMatch[(*i)->qName][(*i)->tName]->posStrand){
            if(!(*i)->flipped){
                (*i)->revComp();
            }
        }
        else{
            if((*i)->flipped){
                (*i)->revComp();
            }
        }

        (*i)->qStart = (*i)->qStart + qOffset[(*i)->qName];
        (*i)->qEnd   = (*i)->qEnd   + qOffset[(*i)->qName];
        (*i)->tStart = (*i)->tStart + tOffset[(*i)->tName];
        (*i)->tEnd   = (*i)->tEnd   + tOffset[(*i)->tName];

        std::cout << **i << std::endl;
    }

    return 0;
}
