#include <map>
#include <string>
#include <vector>
#include <iostream>
#include "chain.hpp"
#include "split.hpp"


int main(int argc, char ** argv)
{
    chain qChain;

    std::vector<std::string> lines;

    std::string line;

    std::string seqid ;

    while(getline(std::cin, line)){

        lines.push_back(line);

        std::vector<std::string> lineDat = split(line, "\t");

        seqid = lineDat[0];

        long int qLen   = atol(lineDat[1].c_str());
        long int qStart = atol(lineDat[2].c_str());
        long int qEnd   = atol(lineDat[3].c_str());
        long int matchB = atol(lineDat[9].c_str());

        if(lineDat[4] == "-"){
            long int tmp = qStart ;
            qStart = qLen - qEnd  ;
            qEnd   = qLen - tmp   ;
        }

        qChain.addAlignment(qStart, qEnd, matchB);
    }
    std::vector<int> indiciesOfAlignments;
    qChain.buildLinks();
    qChain.traceback(indiciesOfAlignments);


    std::cerr << "BEFORE CHAIN: " << lines.size() << " AFTER CHAIN: " << indiciesOfAlignments.size() << " "  << seqid << std::endl;

    for(std::vector<int>::iterator it = indiciesOfAlignments.begin();
        it != indiciesOfAlignments.end() ; it++){
        std::cout << lines[*it] << std::endl;
    }

    return 0;
}
