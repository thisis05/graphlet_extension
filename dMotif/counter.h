#pragma once

#include "graph.h"
#include "digraph.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>

#define DEBUG_PRINT(fmt, ...) printf(fmt, __VA_ARGS__)

struct TraingleInfo {
    Count tri1;
    Count tri2;
    Count tri4;
    
    std::vector<Count> perVertex1;
    std::vector<Count> perEdge1;
    std::vector<Count> perVertex2;
    std::vector<Count> perEdge2_1;
    std::vector<Count> perEdge2_2;
    std::vector<Count> perVertex4;
    std::vector<Count> perEdge4;

    TraingleInfo(VertexIdx nVertices, EdgeIdx nEdges, EdgeIdx nEdges2) 
        : tri1(0), tri2(0), tri4(0),
          perVertex1(nVertices), perEdge1(nEdges2),
          perVertex2(nVertices), perEdge2_1(nEdges), perEdge2_2(nEdges2),
          perVertex4(nVertices), perEdge4(nEdges) {}
};

TraingleInfo getTriangle(CGraph *gout, CGraph *gout_2) {
    TraingleInfo ret(gout->nVertices, gout->nEdges, gout_2->nEdges);

    #pragma omp parallel
    {
        TraingleInfo local_ret(gout->nVertices, gout->nEdges, gout_2->nEdges);

        #pragma omp for schedule(dynamic)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout->nbors[j];
                
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout->nbors[k];
                    
                    EdgeIdx loc_111 = gout->getEdgeBinary(end1, end2);
                    if (loc_111 != -1) {
                        local_ret.tri4++;
                        local_ret.perVertex4[i]++;
                        local_ret.perVertex4[end1]++;
                        local_ret.perVertex4[end2]++;
                        local_ret.perEdge4[j]++;
                        local_ret.perEdge4[k]++;
                        local_ret.perEdge4[loc_111]++;
                    }
                }
                
                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx n2 = gout_2->nbors[k];
                    EdgeIdx loc_122 = (end1 > n2) ? gout_2->getEdgeBinary(n2, end1) : gout_2->getEdgeBinary(end1, n2);
                    
                    if (loc_122 != -1) {
                        local_ret.tri2++;
                        local_ret.perVertex2[i]++;
                        local_ret.perVertex2[end1]++;
                        local_ret.perVertex2[n2]++;
                        local_ret.perEdge2_1[j]++;
                        local_ret.perEdge2_2[k]++;
                        local_ret.perEdge2_2[loc_122]++;
                    }
                }
            }
        }

        #pragma omp for schedule(dynamic)
        for (VertexIdx i = 0; i < gout_2->nVertices; ++i) {
            const EdgeIdx start = gout_2->offsets[i];
            const EdgeIdx end = gout_2->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    
                    EdgeIdx loc_222 = gout_2->getEdgeBinary(end1, end2);
                    if (loc_222 != -1) {
                        local_ret.tri1++;
                        local_ret.perVertex1[i]++;
                        local_ret.perVertex1[end1]++;
                        local_ret.perVertex1[end2]++;
                        local_ret.perEdge1[j]++;
                        local_ret.perEdge1[k]++;
                        local_ret.perEdge1[loc_222]++;
                    } else {
                        EdgeIdx loc_221 = gout->getEdgeBinary(end1, end2);
                        if (loc_221 != -1) {
                            local_ret.tri2++;
                            local_ret.perVertex2[i]++;
                            local_ret.perVertex2[end1]++;
                            local_ret.perVertex2[end2]++;
                            local_ret.perEdge2_2[j]++;
                            local_ret.perEdge2_2[k]++;
                            local_ret.perEdge2_1[loc_221]++;
                        }
                    }
                }
            }
        }

        #pragma omp critical
        {
            ret.tri1 += local_ret.tri1;
            ret.tri2 += local_ret.tri2;
            ret.tri4 += local_ret.tri4;
            
            for (VertexIdx i = 0; i < gout->nVertices; ++i) {
                ret.perVertex1[i] += local_ret.perVertex1[i];
                ret.perVertex2[i] += local_ret.perVertex2[i];
                ret.perVertex4[i] += local_ret.perVertex4[i];
            }
            
            for (EdgeIdx j = 0; j < gout->nEdges; ++j) {
                ret.perEdge2_1[j] += local_ret.perEdge2_1[j];
                ret.perEdge4[j] += local_ret.perEdge4[j];
            }
            
            for (EdgeIdx j = 0; j < gout_2->nEdges; ++j) {
                ret.perEdge1[j] += local_ret.perEdge1[j];
                ret.perEdge2_2[j] += local_ret.perEdge2_2[j];
            }
        }
    }

    return ret;
}

void countAll(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[40]) {
    // 1. Count Three
    double w1 = 0, w2 = 0, t3 = 0;
    VertexIdx n = cg->nVertices;

    #pragma omp parallel for reduction(+:w1,w2,t3)
    for (VertexIdx i = 0; i < n; i++) {
        VertexIdx deg = cg->degree(i);
        VertexIdx deg_2 = cg_2->degree(i);

        w1 += (deg_2 * (deg_2 - 1)) / 2.0;
        w2 += deg * deg_2;
        t3 += (deg * (deg - 1)) / 2.0;
    }

    mcounts[0] = w1;
    mcounts[1] = w2;
    mcounts[4] = t3;

    TraingleInfo tricount = getTriangle(&(dag->outlist), &(dag_2->outlist));
    
    mcounts[2] = tricount.tri1;
    mcounts[3] = tricount.tri2;
    mcounts[5] = tricount.tri4;

    // 2. Count Four 
    // ... (이 부분은 원래 코드에 구현되어 있지 않아 생략했습니다)
}