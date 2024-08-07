#include "counter.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <windows.h>
#include <psapi.h>
#include <iostream>

#define DEBUG_PRINT(fmt, ...) printf(fmt, __VA_ARGS__)



void PrintMemoryUsage() {
    // Get the current process handle
    HANDLE hProcess = GetCurrentProcess();

    // Structure to store memory information
    PROCESS_MEMORY_COUNTERS pmc;

    // Get the memory information for the current process
    if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) {
        std::cout << "Current Memory Usage: " << std::endl;
        std::cout << "Working Set Size: " << pmc.WorkingSetSize / 1024 << " KB" << std::endl;
        std::cout << "Peak Working Set Size: " << pmc.PeakWorkingSetSize / 1024 << " KB" << std::endl;
        std::cout << "Pagefile Usage: " << pmc.PagefileUsage / 1024 << " KB" << std::endl;
        std::cout << "Peak Pagefile Usage: " << pmc.PeakPagefileUsage / 1024 << " KB" << std::endl;

        // Get total physical memory
        MEMORYSTATUSEX memInfo;
        memInfo.dwLength = sizeof(MEMORYSTATUSEX);
        if (GlobalMemoryStatusEx(&memInfo)) {
            DWORDLONG totalPhysicalMemory = memInfo.ullTotalPhys;
            DWORDLONG usedPhysicalMemory = pmc.WorkingSetSize;
            
            std::cout << "Total Physical Memory: " << totalPhysicalMemory / 1024 << " KB" << std::endl;
            std::cout << "Used Physical Memory by Process: " << usedPhysicalMemory / 1024 << " KB" << std::endl;
            std::cout << "Percentage of Total Physical Memory Used by Process: " 
                      << (double(usedPhysicalMemory) / double(totalPhysicalMemory)) * 100.0 << " %" << std::endl;
        } else {
            std::cerr << "Failed to get total physical memory info." << std::endl;
        }
    } else {
        std::cerr << "Failed to get process memory info." << std::endl;
    }

    // Close the process handle
    CloseHandle(hProcess);
}

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
                
                // Tri4 : i < j < k (1, 1, 1)
                //        i->j : 1, i->k : 1, j->k : 1
                
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout->nbors[k];
                    
                    EdgeIdx loc_111 = gout->getEdgeBinary(end1, end2);
                    if (loc_111 != -1) {
                        local_ret.tri4++;
                    }
                }

                // Tri2 : i < j < k 
                //    (1) i->j : 2, i->k : 1, j->k : 2
                //    (2) i->j : 1, i->k : 2, j->k : 2  

                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    EdgeIdx loc221 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                    if (loc221 != -1) {
                        local_ret.tri2++;
                    }
                }                 
            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];

                    // Tri4 : i < j < k 
                    //    (3) i->j : 2, i->k : 2, j->k : 2

                    EdgeIdx loc_222 = gout_2->getEdgeBinary(end1, end2);
                    if (loc_222 != -1) {
                        local_ret.tri1++;
                    } else {
                        
                        // Tri2 : i < j < k 
                        //    (3) i->j : 2, i->k : 2, j->k : 1

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


FourSizeInfo get4size(CGraph *gout, CGraph *gin, CGraph *gout_2, CGraph *gin_2) {

    FourSizeInfo ret (gout->nVertices, gout->nEdges, gout_2->nEdges);
    Count count = 0;
    #pragma omp parallel
    {
        FourSizeInfo local_ret (gout->nVertices, gout->nEdges, gout_2->nEdges);
        
        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];

            TriangleInfo local_tri1(gout_2->degree(i), gout_2->degree(i)-1);
            TriangleInfo local_tri2_1(gout->degree(i), gout_2->degree(i));
            TriangleInfo local_tri2_2_1(gout_2->degree(i), gout->degree(i));
            TriangleInfo local_tri2_2_2(gout_2->degree(i), gout_2->degree(i)-1);
            TriangleInfo local_tri3_1_1(gout->degree(i), gout->degree(i)-1);
            TriangleInfo local_tri3_1_2(gout->degree(i), gout_2->degree(i));
            TriangleInfo local_tri3_2(gout_2->degree(i), gout->degree(i));
            TriangleInfo local_tri4(gout->degree(i), gout->degree(i)-1);
            
            for (EdgeIdx e_j = start; e_j < end; ++e_j) {
                const VertexIdx j = gout->nbors[e_j];
                EdgeIdx idx3_1_1 = local_tri3_1_1.idx;
                EdgeIdx idx3_1_2 = local_tri3_1_2.idx;
                EdgeIdx idx4 = local_tri4.idx;
                EdgeIdx idx2_1 = local_tri2_1.idx;
                Count isIncluded3_1_1 = 0;
                Count isIncluded3_1_2 = 0;
                Count isIncluded4 = 0;
                Count isIncluded2_1 = 0;

                for (EdgeIdx e_k = e_j+1; e_k < end; ++e_k) {
                    const VertexIdx k = gout->nbors[e_k];
                    EdgeIdx loc112 = gout_2->getEdgeBinary(j, k);
                    if (loc112 != -1) {
                        local_ret.tri3++;
                        if (isIncluded3_1_1 == 0){
                            local_tri3_1_1.tri_idx[idx3_1_1] = e_j;
                            isIncluded3_1_1 = 1;
                        }
                        //printf("Debug: i=%lld, idx3_1_1=%lld, local_tri3_1_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx3_1_1, local_tri3_1_1.tri_k_counts[idx3_1_1], local_tri3_1_1.maxEdges, local_tri3_1_1.maxTri);
                        local_tri3_1_1.tri_k[idx3_1_1][local_tri3_1_1.tri_k_counts[idx3_1_1]] = k;
                        local_tri3_1_1.tri_k_counts[idx3_1_1]++;
                    }
                    else{
                        EdgeIdx loc111 = gout->getEdgeBinary(j, k);
                        if (loc111 != -1) {
                            local_ret.tri4++;
                            if (isIncluded4 == 0){
                                local_tri4.tri_idx[idx4] = e_j;
                                isIncluded4 = 1;
                            }
                            //printf("Debug: i=%lld, idx4=%lld, local_tri4 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx4, local_tri4.tri_k_counts[idx4], local_tri4.maxEdges, local_tri4.maxTri);
                            local_tri4.tri_k[idx4][local_tri4.tri_k_counts[idx4]] = k;
                            local_tri4.tri_k_counts[idx4]++;
                        }
                    }
                }

                for (EdgeIdx e_k = start_2; e_k < end_2; ++e_k) {
                    const VertexIdx k = gout_2->nbors[e_k];
                    if (k <= j) {continue;}
                    EdgeIdx loc122 = gout_2->getEdgeBinary(j, k);
                    if (loc122 != -1) {
                        local_ret.tri2++;
                        if (isIncluded2_1 == 0){
                            local_tri2_1.tri_idx[idx2_1] = e_j;
                            isIncluded2_1 = 1;
                        }
                        //printf("Debug: i=%lld, idx2_1=%lld, local_tri2_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx2_1, local_tri2_1.tri_k_counts[idx2_1], local_tri2_1.maxEdges, local_tri2_1.maxTri);
                        local_tri2_1.tri_k[idx2_1][local_tri2_1.tri_k_counts[idx2_1]] = k;
                        local_tri2_1.tri_k_counts[idx2_1]++;
                    }
                    else{
                        EdgeIdx loc121 = gout->getEdgeBinary(j, k);
                        if (loc121 != -1) {
                            local_ret.tri3++;
                            if (isIncluded3_1_2 == 0){
                                local_tri3_1_2.tri_idx[idx3_1_2] = e_j;
                                isIncluded3_1_2 = 1;
                            }
                            //printf("Debug: i=%lld, idx3_1_2=%lld, local_tri3_1_2 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx3_1_2, local_tri3_1_2.tri_k_counts[idx3_1_2], local_tri3_1_2.maxEdges, local_tri3_1_2.maxTri);
                            local_tri3_1_2.tri_k[idx3_1_2][local_tri3_1_2.tri_k_counts[idx3_1_2]] = k;
                            local_tri3_1_2.tri_k_counts[idx3_1_2]++;
                        }
                    }
                }

                for (VertexIdx k1_idx = 0; k1_idx < local_tri2_1.tri_k_counts[idx2_1]; ++k1_idx){
                    VertexIdx k1 = local_tri2_1.tri_k[idx2_1][k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri2_1.tri_k_counts[idx2_1]; ++k2_idx){
                        VertexIdx k2 = local_tri2_1.tri_k[idx2_1][k2_idx];
                        EdgeIdx loc2 = gout_2->getEdgeBinary(k1, k2);
                        if(loc2 != -1){local_ret.clique2++;}
                    }
                }

                if (isIncluded3_1_1 == 1){local_tri3_1_1.idx++;}
                if (isIncluded3_1_2 == 1){local_tri3_1_2.idx++;}
                if (isIncluded4 == 1){local_tri4.idx++;}
                if (isIncluded2_1 == 1){local_tri2_1.idx++;}                
            }
            

            for (EdgeIdx e_j = start_2; e_j < end_2; ++e_j) {
                const VertexIdx j = gout_2->nbors[e_j];
                EdgeIdx idx3_2 = local_tri3_2.idx;
                EdgeIdx idx1 = local_tri1.idx;
                EdgeIdx idx2_2_1 = local_tri2_2_1.idx;
                EdgeIdx idx2_2_2 = local_tri2_2_2.idx;
                Count isIncluded2_2_1 = 0;
                Count isIncluded2_2_2 = 0;
                Count isIncluded1 = 0;
                Count isIncluded3_2 = 0;

                for (EdgeIdx e_k = e_j+1; e_k < end_2; ++e_k) {
                    const VertexIdx k = gout_2->nbors[e_k];
                    EdgeIdx loc_222 = gout_2->getEdgeBinary(j, k);
                    if (loc_222 != -1) {
                        local_ret.tri1++;
                        if (isIncluded1 == 0){
                            local_tri1.tri_idx[idx1] = e_j;
                            isIncluded1 = 1;
                        }
                        local_tri1.tri_k[idx1][local_tri1.tri_k_counts[idx1]] = k;
                        local_tri1.tri_k_counts[idx1]++;
                    }
                    else{
                        EdgeIdx loc221 = gout->getEdgeBinary(j, k);
                        if (loc221 != -1) {
                            local_ret.tri2++;
                            if (isIncluded2_2_2 == 0){
                                local_tri2_2_2.tri_idx[idx2_2_2] = e_j;
                                isIncluded2_2_2 = 1;
                            }
                            local_tri2_2_2.tri_k[idx2_2_2][local_tri2_2_2.tri_k_counts[idx2_2_2]] = k;
                            local_tri2_2_2.tri_k_counts[idx2_2_2]++;
                        }
                    }
                }

                for (EdgeIdx e_k = start; e_k < end; ++e_k) {
                    const VertexIdx k = gout->nbors[e_k];
                    if (k <= j) {continue;}
                    EdgeIdx loc212 = gout_2->getEdgeBinary(j, k);
                    if (loc212 != -1) {
                        local_ret.tri2++;
                        if (isIncluded2_2_1 == 0){
                            local_tri2_2_1.tri_idx[idx2_2_1] = e_j;
                            isIncluded2_2_1 = 1;
                        }
                        //printf("Debug: i=%lld, idx2_2_1=%lld, local_tri2_2_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx2_2_1, local_tri2_2_1.tri_k_counts[idx2_2_1], local_tri2_2_1.maxEdges, local_tri2_2_1.maxTri);
                        local_tri2_2_1.tri_k[idx2_2_1][local_tri2_2_1.tri_k_counts[idx2_2_1]] = k;
                        local_tri2_2_1.tri_k_counts[idx2_2_1]++;
                    }
                    else{
                        EdgeIdx loc211 = gout->getEdgeBinary(j, k);
                        if (loc211 != -1) {
                            local_ret.tri3++;
                            if (isIncluded3_2 == 0){
                                local_tri3_2.tri_idx[idx3_2] = e_j;
                                isIncluded3_2 = 1;
                            }
                            local_tri3_2.tri_k[idx3_2][local_tri3_2.tri_k_counts[idx3_2]] = k;
                            local_tri3_2.tri_k_counts[idx3_2]++;
                        }
                    }
                }
                //d1-1
                for (VertexIdx k1_idx = 0; k1_idx < local_tri1.tri_k_counts[idx1]; ++k1_idx){
                    VertexIdx k1 = local_tri1.tri_k[idx1][k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri1.tri_k_counts[idx1]; ++k2_idx){
                        VertexIdx k2 = local_tri1.tri_k[idx1][k2_idx];
                        EdgeIdx loc2 = gout_2->getEdgeBinary(k1, k2);
                        if(loc2 != -1){local_ret.clique1++;}
                        else{
                            EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                            if(loc1 != -1) {local_ret.clique2++;}
                        }
                    }
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri2_2_1.tri_k_counts[idx2_2_1]; ++k2_idx){
                        VertexIdx k2 = local_tri2_2_1.tri_k[idx2_2_1][k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique2++;}
                    }
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri2_2_2.tri_k_counts[idx2_2_2]; ++k2_idx){
                        VertexIdx k2 = local_tri2_2_2.tri_k[idx2_2_2][k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique2++;}
                    }
                } 
                
                if (isIncluded1 == 1){local_tri1.idx++;}
                if (isIncluded3_2 == 1){local_tri3_2.idx++;}
                if (isIncluded2_2_1 == 1){local_tri2_2_1.idx++;}
                if (isIncluded2_2_2 == 1){local_tri2_2_2.idx++;}         
            }      
        }
        

        #pragma omp critical
        {
            ret.tri1 += local_ret.tri1;
            ret.tri2 += local_ret.tri2;
            ret.tri3 += local_ret.tri3;
            ret.tri4 += local_ret.tri4;
            ret.clique1 += local_ret.clique1;
            ret.clique2 += local_ret.clique2;
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
    
    FourSizeInfo motifcounts = get4size(&(dag->outlist), &(dag->inlist), &(dag_2->outlist), &(dag_2->inlist));
    mcounts[2] = motifcounts.tri1;
    mcounts[3] = motifcounts.tri2;
    mcounts[4] = motifcounts.tri3;
    mcounts[5] = motifcounts.tri4;
    
    mcounts[6] = motifcounts.clique1;
    mcounts[7] = motifcounts.clique2;
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
    printf("Tri1 : %.1f\n", mcounts[2]);
    printf("Tri2 : %.1f\n", mcounts[3]);
    printf("Tri3 : %.1f\n", mcounts[4]);
    printf("Tri4 : %.1f\n", mcounts[5]);
    printf("d1-1 : %.1f\n", mcounts[6]);
    printf("d1-2 : %.1f\n", mcounts[7]);
}