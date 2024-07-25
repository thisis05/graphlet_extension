#include "counter.h"
#include <cstdio>


TriangleInfo getTriangle(Graph* graph) {
    TriangleInfo ret;
    ret.total = 0;
    ret.perVertex = new Count[graph->nVertices]();
    ret.perEdge = new Count[graph->nEdges]();

    for (VertexIdx i = 0; i < graph->nVertices; ++i) {
        EdgeIdx start = graph->out_offsets[i];
        EdgeIdx end = graph->out_offsets[i+1];
        
        for (EdgeIdx j = start; j < end; ++j) {
            for (EdgeIdx k = j + 1; k < end; ++k) {
                VertexIdx end1 = graph->out_dsts[j];
                VertexIdx end2 = graph->out_dsts[k];

                EdgeIdx loc = graph->getEdgeBinary(end1, end2);
                
                if (loc != -1) {
                    ret.total++;

                    ret.perVertex[i]++;
                    ret.perVertex[end1]++;
                    ret.perVertex[end2]++;

                    ret.perEdge[j]++;
                    ret.perEdge[k]++;
                    ret.perEdge[loc]++;
                }
            }
        }
    }

    return ret;
}
void countThree(Graph* dag, double (&mcounts)[6]){

    double w1 = 0, w2 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0; // # of 3-size d-Motifs : wedge 1, 2 & Triangle 1, 2, 3, 4 
    VertexIdx n = dag->nVertices;

    for (VertexIdx i = 0; i < n; i++){
        VertexIdx deg = dag->verticesDegree[i];
        t3 += deg * (deg-1) / 2; 
    }

    mcounts[4] = t3;
    TriangleInfo t4count = getTriangle(dag);
    mcounts[5] = t4count.total;

    printf("Tri3 : %.1f\n", mcounts[4] - 3*mcounts[5]);
    printf("Tri4 : %.1f\n", mcounts[5]);
}
