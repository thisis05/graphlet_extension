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
    
    Count *perVertex1;
    Count *perEdge1;
    Count *perVertex2;
    Count *perEdge2_1;
    Count *perEdge2_2;
    Count *perVertex4;
    Count *perEdge4;

    ThreeSizeInfo(VertexIdx nVertices, EdgeIdx nEdges, EdgeIdx nEdges2) 
        : tri1(0), tri2(0), tri3(0), tri4(0) {
        perVertex1 = new Count[nVertices+1]();
        perEdge1 = new Count[nEdges2+1]();
        
        perVertex2 = new Count[nVertices+1]();
        perEdge2_1 = new Count[nEdges+1]();
        perEdge2_2 = new Count[nEdges2+1]();
        
        perVertex4 = new Count[nVertices+1]();
        perEdge4 = new Count[nEdges+1]();
    }

    ~ThreeSizeInfo() {
        delete[] perVertex1;
        delete[] perEdge1;
        delete[] perVertex2;
        delete[] perEdge2_1;
        delete[] perEdge2_2;
        delete[] perVertex4;
        delete[] perEdge4;
    }
};

struct FourSizeInfo {
    Count tri1;
    Count tri2;
    Count tri4;
    Count clique1, clique2, clique3, clique4, clique5, clique6, clique7, clique8, clique9, clique10, clique11;
    Count chord1, chord2, chord3, chord4, chord5, chord6, chord7, chord8; 
    
    Count *perVertex1;
    Count *perEdge1;
    Count *perVertex2;
    Count *perEdge2_1;
    Count *perEdge2_2;
    Count *perVertex4;
    Count *perEdge4;

    FourSizeInfo(VertexIdx nVertices, EdgeIdx nEdges, EdgeIdx nEdges2) 
        : tri1(0), tri2(0), tri4(0),
        clique1(0), clique2(0), clique3(0), clique4(0), clique5(0), clique6(0), clique7(0), clique8(0), clique9(0), clique10(0), clique11(0),
        chord1(0), chord2(0), chord3(0), chord4(0), chord5(0), chord6(0), chord7(0), chord8(0){
        perVertex1 = new Count[nVertices+1]();
        perEdge1 = new Count[nEdges2+1]();
        
        perVertex2 = new Count[nVertices+1]();
        perEdge2_1 = new Count[nEdges+1]();
        perEdge2_2 = new Count[nEdges2+1]();
        
        perVertex4 = new Count[nVertices+1]();
        perEdge4 = new Count[nEdges+1]();
    }

    ~FourSizeInfo() {
        delete[] perVertex1;
        delete[] perEdge1;
        delete[] perVertex2;
        delete[] perEdge2_1;
        delete[] perEdge2_2;
        delete[] perVertex4;
        delete[] perEdge4;
    }
};

ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2);
FourSizeInfo get4size(CGraph *gout, CGraph *gout_2);
void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[6]);
void countFour(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[40]);
void mEquation3(double (&mcounts)[6]);
void mEquation4(double (&mcounts)[40]);
