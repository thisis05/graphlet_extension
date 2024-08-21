#pragma once

#include "digraph.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>

#define DEBUG_PRINT(fmt, ...) printf(fmt, __VA_ARGS__)

struct ThreeSizeInfo {
    Count tri1;
    Count tri2;
    Count tri3;
    Count tri4;

    ThreeSizeInfo(VertexIdx nVertices, EdgeIdx nEdges, EdgeIdx nEdges2) 
        : tri1(0), tri2(0), tri3(0), tri4(0) {}

};

struct TriangleInfo {
    VertexIdx* triend;
    Count count;
    Count maxTri;

    TriangleInfo(Count degree) {
        maxTri = degree;
        count = 0;
        triend = new VertexIdx[maxTri]();
    }

    ~TriangleInfo() {
        delete[] triend;
    }
};

struct FourSizeInfo {
    Count n;
    Count e1;
    Count e2;
    Count star1;
    Count star2;
    Count tri1;
    Count tri2;
    Count tri3;
    Count tri4;
    Count clique1, clique2, clique3, clique4, clique5, clique6, clique7, clique8, clique9, clique10, clique11;
    Count chord1, chord2, chord3, chord4, chord5, chord6, chord7, chord8; 
    Count tailed1, tailed2, tailed3, tailed4, tailed5, tailed6, tailed7, tailed8;
    Count cycle1, cycle2, cycle3;
    Count path1, path2, path3, path4;

    FourSizeInfo(VertexIdx nVertices, EdgeIdx nEdges, EdgeIdx nEdges2) 
        : star1(0), star2(0), tri1(0), tri2(0), tri3(0), tri4(0),
        clique1(0), clique2(0), clique3(0), clique4(0), clique5(0), clique6(0), clique7(0), clique8(0), clique9(0), clique10(0), clique11(0),
        chord1(0), chord2(0), chord3(0), chord4(0), chord5(0), chord6(0), chord7(0), chord8(0),
        tailed1(0), tailed2(0), tailed3(0), tailed4(0), tailed5(0), tailed6(0), tailed7(0), tailed8(0),
        cycle1(0), cycle2(0), cycle3(0),
        path1(0), path2(0), path3(0), path4(0) {
        n = nVertices;
        e1 = nEdges;
        e2 = nEdges2;
    }
};

struct CycleInfo {
    Count cycle1;
    Count cycle2; 
    Count cycle3; 
    Count cycle4;
    CycleInfo() 
    : cycle1(0), cycle2(0), cycle3(0), cycle4(0) {}
};


ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2);
FourSizeInfo get4size(CGraph *gout, CGraph *gin, CGraph *gin_2, CGraph *gout_2);

void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[6]);
void countFour(CDAG *dag, CDAG *dag_2, double (&mcounts)[36]);
void mEquation3(double (&mcounts)[6]);
void mEquation4(double (&mcounts)[36]);
