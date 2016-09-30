#include <map>
#include <string>
#include <vector>
#include <iostream>
#include "chain.hpp"
#include "split.hpp"

void processSeqid(std::vector<std::string> & lines,
                  std::map<std::string, bool> & sentinel){


    chain qChain;

    std::string seqid;

    for(std::vector<std::string>::iterator it = lines.begin();
        it != lines.end(); it++){

        std::vector<std::string> lineDat = split(*it, "\t");

        seqid = lineDat[0];

        if(sentinel.find(seqid) != sentinel.end()){
            std::cerr << "FATAL: file must be sorted by tName" << std::endl;
            exit(1);
        }

        long int qLen   = atol(lineDat[1].c_str());
        long int qStart = atol(lineDat[2].c_str());
        long int qEnd   = atol(lineDat[3].c_str());
        long int matchB = atol(lineDat[9].c_str());

        qChain.addAlignment(qStart, qEnd, matchB);
    }

    std::vector<int> indiciesOfAlignments;
    qChain.buildLinks();
    qChain.traceback(indiciesOfAlignments);

    for(std::vector<int>::iterator it = indiciesOfAlignments.begin();
        it != indiciesOfAlignments.end() ; it++){
        std::cout << lines[*it] << std::endl;
    }

    sentinel[seqid] = true;
}

int main(int argc, char ** argv)
{
    std::string line ;
    std::string seqid;
    std::vector<std::string> lines;
    std::map<std::string, bool> sanity;

    while(getline(std::cin, line)){

        std::vector<std::string> lineDat = split(line, "\t");

        if(lineDat[0] != seqid){
            processSeqid(lines, sanity);
            lines.clear();
            lines.push_back(line);
            seqid = lineDat[0];
        }
        else{
            lines.push_back(line);
        }
    }
    processSeqid(lines, sanity);

    return 0;
}
