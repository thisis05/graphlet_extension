#include "counter.h"
#include <cstdio>


Size3Info getTriangle4(Graph* graph) {
    Size3Info ret;
    ret.total = 0;
    ret.perVertex = new Count[graph->nVertices]();
    ret.perEdge = new Count[graph->nEdges]();
    graph->adjVector2 = new std::map<VertexIdx, std::vector<VertexIdx>>();
    std::map<VertexIdx, std::map<VertexIdx, bool>> wedgeCheck;
    EdgeIdx total_E2 = 0;

    for (VertexIdx i = 0; i < graph->nVertices; ++i) { // i : start node
        EdgeIdx start = graph->out_offsets[i]; // start : start edge index of i 
        EdgeIdx end = graph->out_offsets[i+1]; // end : end edge index of i  
        for (EdgeIdx j = start; j < end; ++j) { 
            VertexIdx end1 = graph->out_dsts[j]; // end1 : first neighbor of i
            VertexIdx end1_nsize_out = graph->out_offsets[end1+1] - graph->out_offsets[end1];
            VertexIdx* found = new VertexIdx[end1_nsize_out]();
            for (EdgeIdx k = j + 1; k < end; ++k) {
                VertexIdx end2 = graph->out_dsts[k]; // end2 : second neighbor of i
                EdgeIdx loc = graph->getEdgeBinary(end1, end2, found);

                if (loc != -1) {
                    ret.total++;

                    ret.perVertex[i]++;
                    ret.perVertex[end1]++;
                    ret.perVertex[end2]++;

                    ret.perEdge[j]++;
                    ret.perEdge[k]++;
                    ret.perEdge[loc]++;
                }
                else{
                    if (wedgeCheck[end1][end2] == false) {
                        wedgeCheck[end1][end2] = true;
                        (*graph->adjVector2)[end1].push_back(end2);
                        (*graph->adjVector2)[end2].push_back(end1);
                        total_E2++;
                    }
                }
            }
            for (VertexIdx n = 0; n < end1_nsize_out; ++n){
                VertexIdx offset = graph->out_offsets[end1];
                if(found[n] == 0 && wedgeCheck[i][end1] == false){
                    wedgeCheck[i][end1] = true;
                    (*graph->adjVector2)[i].push_back(n+offset);
                    (*graph->adjVector2)[n+offset].push_back(i);
                    total_E2++;
                } 
            }
            delete[] found;           
        }
        // EdgeIdx start_in = graph->in_offsets[i]; 
        // EdgeIdx end_in = graph->in_offsets[i+1]; 
        // for (EdgeIdx j = start_in; j < end_in; ++j){   
        //     VertexIdx end1= graph->in_dsts[j];
        //     VertexIdx end1_nsize_out = graph->out_offsets[end1+1] - graph->out_offsets[end1];
        //     VertexIdx* found = new VertexIdx[end1_nsize_out]();
        //     for (EdgeIdx k = j + 1; k < end; ++k) {
        //         VertexIdx end2 = graph->in_dsts[k]; // end2 : second neighbor of i
        //         EdgeIdx loc = graph->getEdgeBinary(end1, end2, found);
        //         if (loc == -1) {
        //             if (wedgeCheck[end1][end2] == false) {
        //                 wedgeCheck[end1][end2] = true;
        //                 (*graph->adjVector2)[end1].push_back(end2);
        //                 (*graph->adjVector2)[end2].push_back(end1);
        //                 total_E2++;
        //             }
        //         }
        //     }
        //     delete[] found;    
        // }       
    }  
    graph->nEdges2 = total_E2;
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
    Size3Info t4count = getTriangle4(dag);
    mcounts[5] = t4count.total;
    
    (*dag).createE2attr();

    printf("Tri3 : %.1f\n", mcounts[4] - 3*mcounts[5]);
    printf("Tri4 : %.1f\n", mcounts[5]);
}
