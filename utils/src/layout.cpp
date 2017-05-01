#include <unistd.h>
#include <stdlib.h>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "alignment.hpp"
#include "split.hpp"

struct opts{
    bool include ;
    bool ordered ;
    bool qi ;
    std::vector<std::string> tOrder;
    std::map<std::string, bool> toInclude;

}globalOpts;

static const char *optString = "i:hq";

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
        return al1->tStart <= al2->tStart;
    }
    return al1->tName < al2->tName;
}

bool sortQ(const alignment * al1, const alignment * al2 ){
    if(al1->qName == al2->qName){
        return al1->qStart <= al2->qStart;
    }
    return al1->qName < al2->qName;
}

bool loadData(qtCount & qtMatch,
              std::vector<alignment *> & records_tSorted,
              std::map<std::string, long int> & tLens){

    std::string line;

    while(getline(std::cin, line)){

        std::vector<std::string> ld = split(line, "\t");

        if(globalOpts.include){
            if(globalOpts.toInclude.find(ld[5])
               == globalOpts.toInclude.end() ){
                continue;
            }
        }

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
            al->revCompQ();
        }
        records_tSorted.push_back(al);
    }
    return true;
}

void printHelp(void){
    std::cerr << "Layout provides a simple algorithm to layout minimap" << std::endl;
    std::cerr << "alignments for plotting. It works on both chained" << std::endl ;
    std::cerr << "and unchained alignments." << std::endl << std::endl;
    std::cerr << "Usage:" << std::endl;
    std::cerr << " layout -h" << std::endl;
    std::cerr << " cat minimap.txt | layout -i [target_seqid_order] " << std::endl;
    std::cerr << " cat minimap.txt | layout " << std::endl << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << " -h  Show this message." << std::endl;
    std::cerr << " -i  A comma seperated list of targets.     " << std::endl;
    std::cerr << " -q  flag no query offset.                  " << std::endl;
    std::cerr << "     This also dictates target order.       " << std::endl;
}

int parseOpts(int argc, char** argv)
{
    int opt = 0;
    opt = getopt(argc, argv, optString);
    while(opt != -1){
        switch(opt){
        case 'i':
            {

                globalOpts.include = true;
                globalOpts.ordered = true;

                std::vector<std::string> include
                    = split((std::string)optarg, ",");

                for(std::vector<std::string>::iterator it = include.begin();
                    it != include.end(); it++){
                    globalOpts.toInclude[*it] = true;
                    globalOpts.tOrder.push_back(*it);
                }

                if(globalOpts.toInclude.empty()){
                    std::cerr
                        << "FATAL: check -i, it should be a comma sep list."
                        << std::endl;
                    return 0;
                }
                break;
            }
        case 'q':
            {
                globalOpts.qi = false;
                break;
            }
        case 'h':

            {
                return 0;
            }
        default:
            {
                std::cerr << "FATAL: unable to parse command line correctly."
                          << std::endl;
                return 0;
            }
        }
        opt = getopt( argc, argv, optString );
    }
    return 1;
}

int main(int argc, char ** argv)
{
    globalOpts.include  = false;
    globalOpts.ordered  = false;
    globalOpts.qi       = true ;

    int parse = parseOpts(argc, argv);
    if(parse != 1){
        printHelp();
        exit(1);
    }

    std::vector<alignment *> records_tSorted;

    qtCount qtMatch;

    std::map<std::string, std::string> bestTarget;
    std::map<std::string, long int>         tLens;

    loadData(qtMatch, records_tSorted, tLens);

    if(globalOpts.ordered == false){
        for(std::map<std::string, long int>::iterator it = tLens.begin();
            it != tLens.end() ; it++){
            globalOpts.tOrder.push_back(it->first);
        }
    }

    std::cerr << "INFO: Loaded "  << records_tSorted.size()
              << " alignments " << std::endl;
    std::cerr << "INFO: Sorting alignment blocks." << std::endl;

    std::sort(records_tSorted.begin(), records_tSorted.end(), sortT);

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

    for(std::vector<std::string>::iterator it = globalOpts.tOrder.begin();
        it != globalOpts.tOrder.end(); it++){
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

    for(std::vector<std::string>::iterator it = globalOpts.tOrder.begin();
        it != globalOpts.tOrder.end(); it++){
        if(tLens.find(*it) == tLens.end()){
            std::cerr << "FATAL: " << *it << " not found " << std::endl;
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
                (*i)->revCompQ();
            }
        }
        else{
            if((*i)->flipped){
                (*i)->revCompQ();
            }
        }

        if(globalOpts.qi){
            (*i)->qStart = (*i)->qStart + qOffset[(*i)->qName];
            (*i)->qEnd   = (*i)->qEnd   + qOffset[(*i)->qName];
        }
        else{
            (*i)->qStart = (*i)->qStart ;
            (*i)->qEnd   = (*i)->qEnd   ;
        }
        if(globalOpts.tOrder.size() > 1){
            (*i)->tStart = (*i)->tStart + tOffset[(*i)->tName];
            (*i)->tEnd   = (*i)->tEnd   + tOffset[(*i)->tName];
        }
        if((*i)->strand == '-'){
            long int tmp = (*i)->qStart;
            (*i)->qStart = (*i)->qEnd;
            (*i)->qEnd   = tmp;
        }

        std::cout << **i << std::endl;
    }

    for(std::vector<alignment *>::iterator it = records_tSorted.begin();
        it != records_tSorted.end(); it++){
        delete (*it);
    }

    for(qtCount::iterator it = qtMatch.begin();
        it != qtMatch.end(); it++){

        for(std::map<std::string, matchInfo * >::iterator iz
                = it->second.begin();
            iz != it->second.end(); iz++){

            delete iz->second;

        }
    }
    return 0;
}
