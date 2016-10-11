#ifndef CHAINING_H
#define CHAINING_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

class node{
public:
    double matches       ;
    double overallScore  ;
    int start            ;
    int end              ;
    int index            ;

    std::vector<node *> children;

};

class chaining{
private:
    int current_index        ;
    node               * last;
    std::vector<node *> nodes;

public:
    chaining(void);
    ~chaining(void);
    bool buildLinks(  void    );
    bool addAlignment(int, int, int         );
    bool traceback(std::vector<int> &       );
};

    bool _endCmp(const node *, const node * );

#endif
