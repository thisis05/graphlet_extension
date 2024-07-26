#pragma once

#include "graph.h"
#include "algorithm"

using Count = int64_t;

struct Size3Info {
    Count total;
    Count* perVertex;
    Count* perEdge;
};

Size3Info getTriangle4(Graph* graph);
void countThree(Graph* dag, double (&mcounts)[6]);