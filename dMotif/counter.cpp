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
    #pragma omp parallel
    {
        FourSizeInfo local_ret (gout->nVertices, gout->nEdges, gout_2->nEdges);
        
        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];

            VertexIdx *local_star1_outout = new VertexIdx[gout->nVertices]();
            VertexIdx *local_star1_outin = new VertexIdx[gout->nVertices]();
            VertexIdx *local_star2_outout = new VertexIdx[gout->nVertices]();
            VertexIdx *local_star2_outin = new VertexIdx[gout->nVertices]();
            
            for (EdgeIdx e_j = start; e_j < end; ++e_j) {
                const VertexIdx j = gout->nbors[e_j];
                

                TriangleInfo local_tri2_1(gout_2->degree(i));
                TriangleInfo local_tri3_1_1(gout->degree(i)-1);
                TriangleInfo local_tri3_1_2(gout_2->degree(i));
                TriangleInfo local_tri4(gout->degree(i)-1);

                for (EdgeIdx e_k = e_j+1; e_k < end; ++e_k) {
                    const VertexIdx k = gout->nbors[e_k];
                    EdgeIdx loc112 = gout_2->getEdgeBinary(j, k);
                    if (loc112 != -1) {
                        local_ret.tri3++;
                        //printf("Debug: i=%lld, idx3_1_1=%lld, local_tri3_1_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx3_1_1, local_tri3_1_1.tri_k_counts[idx3_1_1], local_tri3_1_1.maxEdges, local_tri3_1_1.maxTri);
                        local_tri3_1_1.triend[local_tri3_1_1.count] = k;
                        local_tri3_1_1.count++;
                    }
                    else{
                        EdgeIdx loc111 = gout->getEdgeBinary(j, k);
                        if (loc111 != -1) {
                            local_ret.tri4++;
                            //printf("Debug: i=%lld, idx4=%lld, local_tri4 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx4, local_tri4.tri_k_counts[idx4], local_tri4.maxEdges, local_tri4.maxTri);
                            local_tri4.triend[local_tri4.count] = k;
                            local_tri4.count++;
                        }
                    }
                }

                for (EdgeIdx e_k = start_2; e_k < end_2; ++e_k) {
                    const VertexIdx k = gout_2->nbors[e_k];
                    if (k <= j) {continue;}
                    EdgeIdx loc122 = gout_2->getEdgeBinary(j, k);
                    if (loc122 != -1) {
                        local_ret.tri2++;
                        //printf("Debug: i=%lld, idx2_1=%lld, local_tri2_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx2_1, local_tri2_1.tri_k_counts[idx2_1], local_tri2_1.maxEdges, local_tri2_1.maxTri);
                        local_tri2_1.triend[local_tri2_1.count] = k;
                        local_tri2_1.count++;
                    }
                    else{
                        EdgeIdx loc121 = gout->getEdgeBinary(j, k);
                        if (loc121 != -1) {
                            local_ret.tri3++;
                            //printf("Debug: i=%lld, idx3_1_2=%lld, local_tri3_1_2 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx3_1_2, local_tri3_1_2.tri_k_counts[idx3_1_2], local_tri3_1_2.maxEdges, local_tri3_1_2.maxTri);
                            local_tri3_1_2.triend[local_tri3_1_2.count] = k;
                            local_tri3_1_2.count++;
                        }
                    }
                }

                // 1. Tri4-based
                for (VertexIdx k1_idx = 0; k1_idx < local_tri4.count; ++k1_idx){
                    VertexIdx k1 = local_tri4.triend[k1_idx];
                    // Tri4 - Tri4
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri4.count; ++k2_idx){
                        VertexIdx k2 = local_tri4.triend[k2_idx];
                        EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1){local_ret.clique11++;}
                        else{local_ret.clique10++;}
                    }
                    // Tri4 - Tri3_1_1
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri3_1_1.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_1_1.triend[k2_idx];
                        EdgeIdx loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                        if(loc1 != -1){local_ret.clique10++;}
                        else{local_ret.clique9++;}
                    }
                    //Tri4 - Tri3_1_2
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri3_1_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_1_2.triend[k2_idx];
                        EdgeIdx loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                        if(loc1 != -1) {local_ret.clique10++;}
                        else{local_ret.clique9++;}
                    }
                    //Tri4 - Tri2_1
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri2_1.count; ++k2_idx){
                        VertexIdx k2 = local_tri2_1.triend[k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique7++;}
                        else{
                            EdgeIdx loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique9++;}
                            else{local_ret.chord7++;}
                        }
                    }
                }

                // 2. Tri2_1-based
                for (VertexIdx k1_idx = 0; k1_idx < local_tri2_1.count; ++k1_idx){
                    VertexIdx k1 = local_tri2_1.triend[k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri2_1.count; ++k2_idx){
                        VertexIdx k2 = local_tri2_1.triend[k2_idx];
                        EdgeIdx loc2 = gout_2->getEdgeBinary(k1, k2);
                        if(loc2 != -1){local_ret.clique2++;}
                        else{
                            EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                            if(loc1 != -1) {local_ret.clique3++;}
                            else{local_ret.chord2++;}
                        }
                    }
                    // Tri2_1 - Tri3_1_1
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri3_1_1.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_1_1.triend[k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique4++;}
                        else{
                            EdgeIdx loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique5++;}
                            else{local_ret.chord5_2++;}
                        }
                    }
                    //Tri2_1 - Tri3_1_2
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri3_1_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_1_2.triend[k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique4++;}
                        else{
                            EdgeIdx loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique5++;}
                            else{local_ret.chord5_2++;}
                        }
                    }
                }

                // 3. Tri3_1_1 based
                for (VertexIdx k1_idx = 0; k1_idx < local_tri3_1_1.count; ++k1_idx){
                    VertexIdx k1 = local_tri3_1_1.triend[k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri3_1_1.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_1_1.triend[k2_idx];
                        EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique9++;}
                        else{local_ret.clique6++;}
                    }
                    // Tri3_1_1 - Tri3_1_2
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri3_1_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_1_2.triend[k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique5++;}
                        else{
                            EdgeIdx loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique8++;}
                            else{local_ret.chord8++;}
                        }
                    }
                }

                // 3. Tri3_1_2 based
                for (VertexIdx k1_idx = 0; k1_idx < local_tri3_1_2.count; ++k1_idx){
                    VertexIdx k1 = local_tri3_1_2.triend[k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri3_1_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_1_2.triend[k2_idx];
                        EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique9++;}
                        else{local_ret.clique6++;}
                    }
                }           
            }
            

            for (EdgeIdx e_j = start_2; e_j < end_2; ++e_j) {
                const VertexIdx j = gout_2->nbors[e_j];

                TriangleInfo local_tri1(gout_2->degree(i)-1);
                TriangleInfo local_tri3_2(gout->degree(i));
                TriangleInfo local_tri2_2_1(gout->degree(i));
                TriangleInfo local_tri2_2_2(gout_2->degree(i)-1);

                for (EdgeIdx e_k = e_j+1; e_k < end_2; ++e_k) {
                    const VertexIdx k = gout_2->nbors[e_k];
                    EdgeIdx loc_222 = gout_2->getEdgeBinary(j, k);
                    if (loc_222 != -1) {
                        local_ret.tri1++;
                        local_tri1.triend[local_tri1.count]= k;
                        local_tri1.count++;
                    }
                    else{
                        EdgeIdx loc221 = gout->getEdgeBinary(j, k);
                        if (loc221 != -1) {
                            local_ret.tri2++;
                            local_tri2_2_2.triend[local_tri2_2_2.count] = k;
                            local_tri2_2_2.count++;
                        }
                    }
                }

                for (EdgeIdx e_k = start; e_k < end; ++e_k) {
                    const VertexIdx k = gout->nbors[e_k];
                    if (k <= j) {continue;}
                    EdgeIdx loc212 = gout_2->getEdgeBinary(j, k);
                    if (loc212 != -1) {
                        local_ret.tri2++;
                        //printf("Debug: i=%lld, idx2_2_1=%lld, local_tri2_2_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx2_2_1, local_tri2_2_1.tri_k_counts[idx2_2_1], local_tri2_2_1.maxEdges, local_tri2_2_1.maxTri);
                        local_tri2_2_1.triend[local_tri2_2_1.count] = k;
                        local_tri2_2_1.count++;
                    }
                    else{
                        EdgeIdx loc211 = gout->getEdgeBinary(j, k);
                        if (loc211 != -1) {
                            local_ret.tri3++;
                            local_tri3_2.triend[local_tri3_2.count] = k;
                            local_tri3_2.count++;
                        }
                    }
                }

                // get Cycle
                // for (EdgeIdx e_k = gout_2->offsets[j]; e_k < gout_2->offsets[j+1]; ++e_k) {
                //     const VertexIdx k = gout_2->nbors[e_k];
                //     local_star1_outout[k]++;
                // }
                // for (EdgeIdx e_k = gin_2->offsets[j]; e_k < gin_2->offsets[j+1]; ++e_k) {
                //     const VertexIdx k = gin_2->nbors[e_k];
                //     local_star1_outin[k]++;
                // }

                // 1. Tri1-based
                for (VertexIdx k1_idx = 0; k1_idx < local_tri1.count; ++k1_idx){
                    VertexIdx k1 = local_tri1.triend[k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri1.count; ++k2_idx){
                        VertexIdx k2 = local_tri1.triend[k2_idx];
                        EdgeIdx loc2 = gout_2->getEdgeBinary(k1, k2);
                        if(loc2 != -1){local_ret.clique1++;}
                        else{
                            EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                            if(loc1 != -1) {local_ret.clique2++;}
                            else{local_ret.chord1++;}
                        }
                    }
                    // Tri1 - Tri2_2_1
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri2_2_1.count; ++k2_idx){
                        VertexIdx k2 = local_tri2_2_1.triend[k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique2++;}
                        else{
                            EdgeIdx loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique4++;}
                            else{local_ret.chord3++;}
                        }
                    }
                    // Tri1 - Tri2_2_2
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri2_2_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri2_2_2.triend[k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique2++;}
                        else{
                            EdgeIdx loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique4++;}
                            else{local_ret.chord3++;}
                        }
                    }
                    // Tri1 - Tri3_2
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri3_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_2.triend[k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique4++;}
                        else{
                            EdgeIdx loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique6++;}
                            else{local_ret.chord5_1++;}
                        }
                    }
                } 

                // 2. Tri2_2_1-based
                for (VertexIdx k1_idx = 0; k1_idx < local_tri2_2_1.count; ++k1_idx){
                    VertexIdx k1 = local_tri2_2_1.triend[k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri2_2_1.count; ++k2_idx){
                        VertexIdx k2 = local_tri2_2_1.triend[k2_idx];
                        EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique7++;}
                        else{local_ret.clique4++;}
                    }
                    // Tri2_2_1 - Tri2_2_2
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri2_2_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri2_2_2.triend[k2_idx];
                        EdgeIdx loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique3++;}
                        else{
                            EdgeIdx loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique5++;}
                            else{local_ret.chord4++;}
                        }
                    }
                    // Tri2_2_1 - Tri3_2
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri3_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_2.triend[k2_idx];
                        EdgeIdx loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                        if(loc1 != -1) {local_ret.clique9++;}
                        else{local_ret.clique5++;}
                    }
                }

                // 3. Tri2_2_2-based
                for (VertexIdx k1_idx = 0; k1_idx < local_tri2_2_2.count; ++k1_idx){
                    VertexIdx k1 = local_tri2_2_2.triend[k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri2_2_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri2_2_2.triend[k2_idx];
                        EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique7++;}
                        else{local_ret.clique4++;}
                    }
                    // Tri2_2_2 - Tri3_2
                    for (VertexIdx k2_idx = 0; k2_idx < local_tri3_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_2.triend[k2_idx];
                        EdgeIdx loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                        if(loc1 != -1) {local_ret.clique9++;}
                        else{local_ret.clique5++;}
                    }
                }

                // 4. Tri3_2-based
                for (VertexIdx k1_idx = 0; k1_idx < local_tri3_2.count; ++k1_idx){
                    VertexIdx k1 = local_tri3_2.triend[k1_idx];
                    for (VertexIdx k2_idx = k1_idx+1; k2_idx < local_tri3_2.count; ++k2_idx){
                        VertexIdx k2 = local_tri3_2.triend[k2_idx];
                        EdgeIdx loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique10++;}
                        else{local_ret.clique8++;}
                    }
                }  
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
            ret.clique3 += local_ret.clique3;
            ret.clique4 += local_ret.clique4;
            ret.clique5 += local_ret.clique5;
            ret.clique6 += local_ret.clique6;
            ret.clique7 += local_ret.clique7;
            ret.clique8 += local_ret.clique8;
            ret.clique9 += local_ret.clique9;
            ret.clique10 += local_ret.clique10;
            ret.clique11 += local_ret.clique11;
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
    mcounts[8] = motifcounts.clique3;
    mcounts[9] = motifcounts.clique4;
    mcounts[10] = motifcounts.clique5;
    mcounts[11] = motifcounts.clique6;
    mcounts[12] = motifcounts.clique7;
    mcounts[13] = motifcounts.clique8;
    mcounts[14] = motifcounts.clique9;
    mcounts[15] = motifcounts.clique10;
    mcounts[16] = motifcounts.clique11;
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
    printf("d1-3 : %.1f\n", mcounts[8]);
    printf("d1-4 : %.1f\n", mcounts[9]);
    printf("d1-5 : %.1f\n", mcounts[10]);
    printf("d1-6 : %.1f\n", mcounts[11]);
    printf("d1-7 : %.1f\n", mcounts[12]);
    printf("d1-8 : %.1f\n", mcounts[13]);
    printf("d1-9 : %.1f\n", mcounts[14]);
    printf("d1-10 : %.1f\n", mcounts[15]);
    printf("d1-11 : %.1f\n", mcounts[16]);
}