#pragma once

#include "graph.h"

using Count = int64_t;

struct TriangleInfo {
    Count total;
    Count* perVertex;
    Count* perEdge;
};

TriangleInfo getTriangle(Graph* graph);
void countThree(Graph* dag, double (&mcounts)[6]);