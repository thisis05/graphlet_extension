#include "counter.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>

#define DEBUG_PRINT(fmt, ...) printf(fmt, __VA_ARGS__)

ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2) {

    ThreeSizeInfo ret (gout->nVertices, gout->nEdges, gout_2->nEdges);

    #pragma omp parallel
    {
        ThreeSizeInfo local_ret (gout->nVertices, gout->nEdges, gout_2->nEdges);
        
        #pragma omp for schedule(dynamic, 64)
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
                    }
                }
                
                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx n2 = gout_2->nbors[k];
                    EdgeIdx loc_122 = (end1 > n2) ? gout_2->getEdgeBinary(n2, end1) : gout_2->getEdgeBinary(end1, n2);
                    
                    if (loc_122 != -1) {
                        local_ret.tri2++;
                    }
                }
            }
        }

        #pragma omp for schedule(dynamic, 64)
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
                    } else {
                        EdgeIdx loc_221 = gout->getEdgeBinary(end1, end2);
                        if (loc_221 != -1) {
                            local_ret.tri2++;
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
        }
    }

    return ret;
}

FourSizeInfo get4size(CGraph *gout, CGraph *gout_2) {

    FourSizeInfo ret (gout->nVertices, gout->nEdges, gout_2->nEdges);

    VertexIdx step = 0;
    #pragma omp parallel
    {
        FourSizeInfo local_ret (gout->nVertices, gout->nEdges, gout_2->nEdges);
        
        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout->nbors[j];
                
                Count count = 0;
                Count *tri4ends = new Count[end-start]();
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
                        tri4ends[count] = end2;
                        count++;
                    }
                }
                
                for (VertexIdx posk = 0; posk < count; ++posk) // loop over all pairs of triangles formed by (i,j)
                {
                    VertexIdx k = tri4ends[posk]; // k is vertex as index posk in triends
                    VertexIdx degk = gout->offsets[k+1] - gout->offsets[k]; // gettting degree of k in gout
                    VertexIdx remaining = count-posk; // number of vertices that k needs to be checked with
                    if (degk >= remaining){   
                        // We will search all other vertices in triends in k's adj list
                        for (VertexIdx posell = posk+1; posell < count; ++posell){
                            VertexIdx ell = tri4ends[posell]; 
                            if (gout->isEdgeBinary(k,ell)){ // (k,ell) is an end, thus (i,j,k,ell) form a 4-clique
                               local_ret.clique11++;
                            }
                            else local_ret.clique10++;
                        }
                    }
                    else{
                        // We will search all vertices in k's adj list in the remaining portion of triends
                        for (EdgeIdx posell = gout->offsets[k]; posell < gout->offsets[k+1]; posell++){
                            VertexIdx ell = gout->nbors[posell];
                            if (binarySearch(tri4ends+posk+1,count-posk-1,ell) != -1){
                                local_ret.clique11++;
                            }
                            else local_ret.clique10++;
                        }
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

        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx i = 0; i < gout_2->nVertices; ++i) {
            //if(step % 100 == 0) printf("Progress : %lld / %lld\n", step, gout->nVertices);
            //step++;
            const EdgeIdx start = gout_2->offsets[i];
            const EdgeIdx end = gout_2->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                
                Count count = 0;
                Count *tri1ends = new Count[end-start]();

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
                        tri1ends[count] = end2;
                        count++;

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
                for (Count k1 = 0; k1 < count; ++k1){
                    VertexIdx e1 = tri1ends[k1];
                    for (Count k2 = k1+1; k2 < count; ++k2){
                        VertexIdx e2 = tri1ends[k2];
                        EdgeIdx clique1 = gout_2->getEdgeBinary(e1, e2);
                        if (clique1 != -1)
                            local_ret.clique1++;
                        else{
                            EdgeIdx clique2 = gout->getEdgeBinary(e1, e2);
                            if (clique2 != -1)
                                local_ret.clique2++;
                            else 
                                local_ret.chord1++;
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
            ret.clique1 += local_ret.clique1;
            ret.clique2 += local_ret.clique2;
            ret.clique3 += local_ret.clique3;
            ret.clique4 += local_ret.clique4;
            ret.clique5 += local_ret.clique5;
            ret.clique6 += local_ret.clique6;
            ret.clique7 += local_ret.clique7;
            ret.clique8 += local_ret.clique8;
            ret.clique9 += local_ret.clique9;
            ret.clique10 += local_ret.clique10;
            ret.clique11 += local_ret.clique11;
            ret.chord1 += local_ret.chord1;

            
            
            // for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            //     ret.perVertex1[i] += local_ret.perVertex1[i];
            //     ret.perVertex2[i] += local_ret.perVertex2[i];
            //     ret.perVertex4[i] += local_ret.perVertex4[i];
            // }
            
            // for (EdgeIdx j = 0; j < gout->nEdges; ++j) {
            //     ret.perEdge2_1[j] += local_ret.perEdge2_1[j];
            //     ret.perEdge4[j] += local_ret.perEdge4[j];
            // }
            
            // for (EdgeIdx j = 0; j < gout_2->nEdges; ++j) {
            //     ret.perEdge1[j] += local_ret.perEdge1[j];
            //     ret.perEdge2_2[j] += local_ret.perEdge2_2[j];
            // }
        }
    }
    return ret;
}


void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[6]){

    // 1. Count Three
    double w1 = 0, w2 = 0, t3 = 0; // # of 3-size d-Motifs : wedge 1, 2 & Triangle 1, 2, 3, 4 
    VertexIdx n = cg->nVertices;

    for (VertexIdx i = 0; i < n; i++){
        VertexIdx deg = cg->degree(i); // degree of i
        VertexIdx deg_2 = cg_2->degree(i); // degree2 of i

        w1 += deg_2 * (deg_2-1) / 2;
        w2 += deg * deg_2; 
        t3 += deg * (deg-1) / 2; 
    }

    mcounts[0] = w1;
    mcounts[1] = w2;
    mcounts[4] = t3;
    ThreeSizeInfo tricount = get3size(&(dag->outlist), &(dag_2->outlist));
    mcounts[2] = tricount.tri1;
    mcounts[3] = tricount.tri2;
    mcounts[5] = tricount.tri4;
}

void countFour(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[40]){

    // 1. Count Three
    double w1 = 0, w2 = 0, t3 = 0; // # of 3-size d-Motifs : wedge 1, 2 & Triangle 1, 2, 3, 4 
    VertexIdx n = cg->nVertices;

    for (VertexIdx i = 0; i < n; i++){
        VertexIdx deg = cg->degree(i); // degree of i
        VertexIdx deg_2 = cg_2->degree(i); // degree2 of i

        w1 += deg_2 * (deg_2-1) / 2;
        w2 += deg * deg_2; 
        t3 += deg * (deg-1) / 2; 
    }

    mcounts[0] = w1;
    mcounts[1] = w2;
    mcounts[4] = t3;
    FourSizeInfo motifcounts = get4size(&(dag->outlist), &(dag_2->outlist));
    mcounts[2] = motifcounts.tri1;
    mcounts[3] = motifcounts.tri2;
    mcounts[5] = motifcounts.tri4;
    
    
    mcounts[6] = motifcounts.clique1;
    mcounts[7] = motifcounts.clique2;
    mcounts[15] = motifcounts.clique10;
    mcounts[16] = motifcounts.clique11;
    mcounts[17] = motifcounts.chord1;
}

void mEquation3(double (&mcounts)[6]){
    double tri3 =  mcounts[4] - 3*mcounts[5];
    printf("Star1 : %.1f\n", mcounts[0] - 3*mcounts[2] - mcounts[3]);
    printf("Star2 : %.1f\n", mcounts[1] - 2*mcounts[3] - 2*tri3);
    printf("Tri1 : %.1f\n", mcounts[2]);
    printf("Tri2 : %.1f\n", mcounts[3]);
    printf("Tri3 : %.1f\n", tri3);
    printf("Tri4 : %.1f\n", mcounts[5]);
}

void mEquation4(double (&mcounts)[40]){
    double tri3 =  mcounts[4] - 3*mcounts[5];
    printf("Star1 : %.1f\n", mcounts[0] - 3*mcounts[2] - mcounts[3]);
    printf("Star2 : %.1f\n", mcounts[1] - 2*mcounts[3] - 2*tri3);
    printf("Tri1 : %.1f\n", mcounts[2]);
    printf("Tri2 : %.1f\n", mcounts[3]);
    printf("Tri3 : %.1f\n", tri3);
    printf("Tri4 : %.1f\n", mcounts[5]);
    printf("d1-1 : %.1f\n", mcounts[6]);
    printf("d1-2 : %.1f\n", mcounts[7]);
    printf("d1-10 : %.1f\n", mcounts[15]);
    printf("d1-11 : %.1f\n", mcounts[16]);
    printf("d2-1 : %.1f\n", mcounts[17]);
}