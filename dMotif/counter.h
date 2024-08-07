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
    EdgeIdx** tri_k;        // 각 엣지에 대한 k 값 배열
    EdgeIdx* tri_k_counts;  // 각 엣지에 대한 k 값 개수
    EdgeIdx* tri_idx; //idx mapping
    Count maxEdges;
    Count maxTri;
    Count idx;

    TriangleInfo(Count edges, Count degree) {
        maxEdges = edges;
        maxTri = degree;
        idx = 0;
        tri_k = new EdgeIdx*[maxEdges];
        tri_idx = new EdgeIdx[maxEdges]();
        tri_k_counts = new EdgeIdx[maxEdges]();

        for (EdgeIdx e = 0; e < maxEdges; ++e) {
            tri_k[e] = new EdgeIdx[maxTri]();
        }
    }

    ~TriangleInfo() {
        for (EdgeIdx e = 0; e < maxEdges; ++e) {
            delete[] tri_k[e];
        }
        delete[] tri_k;
        delete[] tri_idx;
        delete[] tri_k_counts;
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

    FourSizeInfo(VertexIdx nVertices, EdgeIdx nEdges, EdgeIdx nEdges2) 
        : star1(0), star2(0), tri1(0), tri2(0), tri3(0), tri4(0),
        clique1(0), clique2(0), clique3(0), clique4(0), clique5(0), clique6(0), clique7(0), clique8(0), clique9(0), clique10(0), clique11(0),
        chord1(0), chord2(0), chord3(0), chord4(0), chord5(0), chord6(0), chord7(0), chord8(0){
        
        n = nVertices;
        e1 = nEdges;
        e2 = nEdges2;
    }
};


ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2);
FourSizeInfo get4size(CGraph *gout, CGraph *gin, CGraph *gin_2, CGraph *gout_2);

void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[6]);
void countFour(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[40]);
void mEquation3(double (&mcounts)[6]);
void mEquation4(double (&mcounts)[40]);
