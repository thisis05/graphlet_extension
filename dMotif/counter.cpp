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
    
    const Count num_threads = 11;

    omp_set_num_threads(num_threads);
    EdgeIdx** all_local_tri1 = new EdgeIdx*[num_threads];
    EdgeIdx** all_local_tri2_1 = new EdgeIdx*[num_threads];
    EdgeIdx** all_local_tri2_2 = new EdgeIdx*[num_threads];
    EdgeIdx** all_local_tri3_1 = new EdgeIdx*[num_threads];
    EdgeIdx** all_local_tri3_2 = new EdgeIdx*[num_threads];
    EdgeIdx** all_local_tri4 = new EdgeIdx*[num_threads];

    for (int t = 0; t < num_threads; ++t) {
        all_local_tri1[t] = new EdgeIdx[gout_2->nEdges + 1]();
        all_local_tri2_1[t] = new EdgeIdx[gout->nEdges + 1]();
        all_local_tri2_2[t] = new EdgeIdx[gout_2->nEdges + 1]();
        all_local_tri3_1[t] = new EdgeIdx[gout->nEdges + 1]();
        all_local_tri3_2[t] = new EdgeIdx[gout_2->nEdges + 1]();
        all_local_tri4[t] = new EdgeIdx[gout->nEdges + 1]();
    }
    
    #pragma omp parallel
    {
        FourSizeInfo local_ret (gout->nVertices, gout->nEdges, gout_2->nEdges);
        Count thread_id = omp_get_thread_num();

        EdgeIdx e_j, e_k, loc111, loc112, loc121, loc122, loc211, loc212, loc221, loc_222, loc1, loc2;
        VertexIdx i, j, k, k1, k2, k1_idx, k2_idx;

        #pragma omp for schedule(guided)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];

            VertexIdx* local_star1 = new VertexIdx[gout->nVertices+1]();
            VertexIdx* local_star2_1 = new VertexIdx[gout->nVertices+1]();   
            VertexIdx* local_star2_2 = new VertexIdx[gout->nVertices+1]();   

            for (e_j = start; e_j < end; ++e_j) {
                j = gout->nbors[e_j];
                
                TriangleInfo local_tri2_1(gout_2->degree(i));
                TriangleInfo local_tri3_1_1(gout->degree(i)-1);
                TriangleInfo local_tri3_1_2(gout_2->degree(i));
                TriangleInfo local_tri4(gout->degree(i)-1);

                for (e_k = e_j+1; e_k < end; ++e_k) {
                    k = gout->nbors[e_k];
                    loc111 = gout->getEdgeBinary(j, k);
                    if (loc111 != -1) {
                        local_ret.tri4++;
                        //printf("Debug: i=%lld, idx4=%lld, local_tri4 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx4, local_tri4.tri_k_counts[idx4], local_tri4.maxEdges, local_tri4.maxTri);
                        local_tri4.triend[local_tri4.count] = k;
                        local_tri4.count++;
                        all_local_tri4[thread_id][e_j]++;
                        all_local_tri4[thread_id][e_k]++;
                        all_local_tri4[thread_id][loc111]++;
                    }
                    else{
                        loc112 = gout_2->getEdgeBinary(j, k);
                        local_ret.tri3++;
                        //printf("Debug: i=%lld, idx3_1_1=%lld, local_tri3_1_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx3_1_1, local_tri3_1_1.tri_k_counts[idx3_1_1], local_tri3_1_1.maxEdges, local_tri3_1_1.maxTri);
                        local_tri3_1_1.triend[local_tri3_1_1.count] = k;
                        local_tri3_1_1.count++;
                        all_local_tri3_1[thread_id][e_j]++;
                        all_local_tri3_1[thread_id][e_k]++;
                        all_local_tri3_2[thread_id][loc112]++;
                    }
                }

                for (e_k = start_2; e_k < end_2; ++e_k) {
                    k = gout_2->nbors[e_k];
                    if (k <= j) {continue;}
                    loc122 = gout_2->getEdgeBinary(j, k);
                    if (loc122 != -1) {
                        local_ret.tri2++;
                        //printf("Debug: i=%lld, idx2_1=%lld, local_tri2_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx2_1, local_tri2_1.tri_k_counts[idx2_1], local_tri2_1.maxEdges, local_tri2_1.maxTri);
                        local_tri2_1.triend[local_tri2_1.count] = k;
                        local_tri2_1.count++;
                        all_local_tri2_1[thread_id][e_j]++;
                        all_local_tri2_2[thread_id][e_k]++;
                        all_local_tri2_2[thread_id][loc122]++;
                    }
                    else{
                        loc121 = gout->getEdgeBinary(j, k);
                        if (loc121 != -1) {
                            local_ret.tri3++;
                            //printf("Debug: i=%lld, idx3_1_2=%lld, local_tri3_1_2 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx3_1_2, local_tri3_1_2.tri_k_counts[idx3_1_2], local_tri3_1_2.maxEdges, local_tri3_1_2.maxTri);
                            local_tri3_1_2.triend[local_tri3_1_2.count] = k;
                            local_tri3_1_2.count++;
                            all_local_tri3_1[thread_id][e_j]++;
                            all_local_tri3_2[thread_id][e_k]++;
                            all_local_tri3_1[thread_id][loc121]++;
                        }
                    }
                }

                //e1 out - e2 out
                for (e_k = gout_2->offsets[j]; e_k < gout_2->offsets[j+1]; ++e_k) {
                    k = gout_2->nbors[e_k];
                    local_star2_1[k]++;
                }

                //e1 out - e2 in
                for (e_k = gin_2->offsets[j]; e_k < gin_2->offsets[j+1]; ++e_k){
                    k = gin_2->nbors[e_k];
                    if (k > i){local_star2_1[k]++;}
                }

                // 1. Tri4-based
                for (k1_idx = 0; k1_idx < local_tri4.count; ++k1_idx){
                    k1 = local_tri4.triend[k1_idx];
                    // Tri4 - Tri4
                    for (k2_idx = k1_idx+1; k2_idx < local_tri4.count; ++k2_idx){
                        k2 = local_tri4.triend[k2_idx];
                        loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1){local_ret.clique11++;}
                        else{local_ret.clique10++;}
                    }
                    // Tri4 - Tri3_1_1
                    for (k2_idx = 0; k2_idx < local_tri3_1_1.count; ++k2_idx){
                        k2 = local_tri3_1_1.triend[k2_idx];
                        loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                        if(loc1 != -1){local_ret.clique10++;}
                        else{local_ret.clique9++;}
                    }
                    //Tri4 - Tri3_1_2
                    for (k2_idx = 0; k2_idx < local_tri3_1_2.count; ++k2_idx){
                        k2 = local_tri3_1_2.triend[k2_idx];
                        loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                        if(loc1 != -1) {local_ret.clique10++;}
                        else{local_ret.clique9++;}
                    }
                    //Tri4 - Tri2_1
                    for (k2_idx = 0; k2_idx < local_tri2_1.count; ++k2_idx){
                        k2 = local_tri2_1.triend[k2_idx];
                        loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique7++;}
                        else{
                            loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique9++;}
                        }
                    }
                }

                // 2. Tri2_1-based
                for (k1_idx = 0; k1_idx < local_tri2_1.count; ++k1_idx){
                    k1 = local_tri2_1.triend[k1_idx];
                    for (k2_idx = k1_idx+1; k2_idx < local_tri2_1.count; ++k2_idx){
                        k2 = local_tri2_1.triend[k2_idx];
                        loc2 = gout_2->getEdgeBinary(k1, k2);
                        if(loc2 != -1){local_ret.clique2++;}
                        else{
                            loc1 = gout->getEdgeBinary(k1, k2);
                            if(loc1 != -1) {local_ret.clique3++;}
                        }
                    }
                    // Tri2_1 - Tri3_1_1
                    for (k2_idx = 0; k2_idx < local_tri3_1_1.count; ++k2_idx){
                        k2 = local_tri3_1_1.triend[k2_idx];
                        loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique4++;}
                        else{
                            loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique5++;}
                        }
                    }
                    //Tri2_1 - Tri3_1_2
                    for (k2_idx = 0; k2_idx < local_tri3_1_2.count; ++k2_idx){
                        k2 = local_tri3_1_2.triend[k2_idx];
                        loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique4++;}
                        else{
                            loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique5++;}
                        }
                    }
                }

                // 3. Tri3_1_1 based
                for (k1_idx = 0; k1_idx < local_tri3_1_1.count; ++k1_idx){
                    k1 = local_tri3_1_1.triend[k1_idx];
                    for (k2_idx = k1_idx+1; k2_idx < local_tri3_1_1.count; ++k2_idx){
                        k2 = local_tri3_1_1.triend[k2_idx];
                        loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique9++;}
                        else{local_ret.clique6++;}
                    }
                    // Tri3_1_1 - Tri3_1_2
                    for (k2_idx = 0; k2_idx < local_tri3_1_2.count; ++k2_idx){
                        k2 = local_tri3_1_2.triend[k2_idx];
                        loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique5++;}
                        else{
                            loc1 = k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique8++;}
                        }
                    }
                }

                // 3. Tri3_1_2 based
                for (k1_idx = 0; k1_idx < local_tri3_1_2.count; ++k1_idx){
                    k1 = local_tri3_1_2.triend[k1_idx];
                    for (k2_idx = k1_idx+1; k2_idx < local_tri3_1_2.count; ++k2_idx){
                        k2 = local_tri3_1_2.triend[k2_idx];
                        loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique9++;}
                        else{local_ret.clique6++;}
                    }
                }           
            }
            

            for (e_j = start_2; e_j < end_2; ++e_j) {
                j = gout_2->nbors[e_j];

                TriangleInfo local_tri1(gout_2->degree(i)-1);
                TriangleInfo local_tri3_2(gout->degree(i));
                TriangleInfo local_tri2_2_1(gout->degree(i));
                TriangleInfo local_tri2_2_2(gout_2->degree(i)-1);

                for (EdgeIdx e_k = e_j+1; e_k < end_2; ++e_k) {
                    k = gout_2->nbors[e_k];
                    loc_222 = gout_2->getEdgeBinary(j, k);
                    if (loc_222 != -1) {
                        local_ret.tri1++;
                        local_tri1.triend[local_tri1.count]= k;
                        local_tri1.count++;
                        all_local_tri1[thread_id][e_j]++;
                        all_local_tri1[thread_id][e_k]++;
                        all_local_tri1[thread_id][loc_222]++;
                    }
                    else{
                        loc221 = gout->getEdgeBinary(j, k);
                        if (loc221 != -1) {
                            local_ret.tri2++;
                            local_tri2_2_2.triend[local_tri2_2_2.count] = k;
                            local_tri2_2_2.count++;
                            all_local_tri2_2[thread_id][e_j]++;
                            all_local_tri2_2[thread_id][e_k]++;
                            all_local_tri2_1[thread_id][loc221]++;
                        }
                    }
                }

                for (e_k = start; e_k < end; ++e_k) {
                    k = gout->nbors[e_k];
                    if (k <= j) {continue;}
                    loc212 = gout_2->getEdgeBinary(j, k);
                    if (loc212 != -1) {
                        local_ret.tri2++;
                        //printf("Debug: i=%lld, idx2_2_1=%lld, local_tri2_2_1 counts=%lld, maxEdges=%lld, maxTri=%lld\n", i, idx2_2_1, local_tri2_2_1.tri_k_counts[idx2_2_1], local_tri2_2_1.maxEdges, local_tri2_2_1.maxTri);
                        local_tri2_2_1.triend[local_tri2_2_1.count] = k;
                        local_tri2_2_1.count++;
                        all_local_tri2_2[thread_id][e_j]++;
                        all_local_tri2_1[thread_id][e_k]++;
                        all_local_tri2_2[thread_id][loc212]++;
                    }
                    else{
                        loc211 = gout->getEdgeBinary(j, k);
                        if (loc211 != -1) {
                            local_ret.tri3++;
                            local_tri3_2.triend[local_tri3_2.count] = k;
                            local_tri3_2.count++;
                            all_local_tri3_2[thread_id][e_j]++;
                            all_local_tri3_1[thread_id][e_k]++;
                            all_local_tri3_1[thread_id][loc211]++;
                        }
                    }
                }

                //e2 out - e2 out
                for (e_k = gout_2->offsets[j]; e_k < gout_2->offsets[j+1]; ++e_k) {
                    k = gout_2->nbors[e_k];
                    local_star1[k]++;
                }
                
                //e2 out - e2 in
                for (e_k = gin_2->offsets[j]; e_k < gin_2->offsets[j+1]; ++e_k) {
                    k = gin_2->nbors[e_k];
                    if (k > i){local_star1[k]++;}
                }

                //e2 out - e1 out
                for (e_k = gout->offsets[j]; e_k < gout->offsets[j+1]; ++e_k) {
                    k = gout->nbors[e_k];
                    local_star2_2[k]++;
                }

                //e2 out - e1 in
                for (e_k = gin->offsets[j]; e_k < gin->offsets[j+1]; ++e_k) {
                    k = gin->nbors[e_k];
                    if (k > i){local_star2_2[k]++;}
                }
                

                // 1. Tri1-based
                for (k1_idx = 0; k1_idx < local_tri1.count; ++k1_idx){
                    k1 = local_tri1.triend[k1_idx];
                    for (k2_idx = k1_idx+1; k2_idx < local_tri1.count; ++k2_idx){
                        k2 = local_tri1.triend[k2_idx];
                        local_ret.chord1++;
                        loc2 = gout_2->getEdgeBinary(k1, k2);
                        if(loc2 != -1){local_ret.clique1++;}
                        else{
                            loc1 = gout->getEdgeBinary(k1, k2);
                            if(loc1 != -1) {local_ret.clique2++;}
                        }
                    }
                    // Tri1 - Tri2_2_1
                    for (k2_idx = 0; k2_idx < local_tri2_2_1.count; ++k2_idx){
                        k2 = local_tri2_2_1.triend[k2_idx];
                        loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique2++;}
                        else{
                            loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique4++;}
                        }
                    }
                    // Tri1 - Tri2_2_2
                    for (k2_idx = 0; k2_idx < local_tri2_2_2.count; ++k2_idx){
                        k2 = local_tri2_2_2.triend[k2_idx];
                        loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique2++;}
                        else{
                            loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique4++;}
                        }
                    }
                    // Tri1 - Tri3_2
                    for (k2_idx = 0; k2_idx < local_tri3_2.count; ++k2_idx){
                        k2 = local_tri3_2.triend[k2_idx];
                        loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique4++;}
                        else{
                            loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique6++;}
                        }
                    }
                } 

                // 2. Tri2_2_1-based
                for (k1_idx = 0; k1_idx < local_tri2_2_1.count; ++k1_idx){
                    k1 = local_tri2_2_1.triend[k1_idx];
                    for (k2_idx = k1_idx+1; k2_idx < local_tri2_2_1.count; ++k2_idx){
                        k2 = local_tri2_2_1.triend[k2_idx];
                        loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique7++;}
                        else{local_ret.clique4++;}
                    }
                    // Tri2_2_1 - Tri2_2_2
                    for (k2_idx = 0; k2_idx < local_tri2_2_2.count; ++k2_idx){
                        k2 = local_tri2_2_2.triend[k2_idx];
                        loc2 = k1 < k2 ? gout_2->getEdgeBinary(k1, k2): gout_2->getEdgeBinary(k2, k1);
                        if(loc2 != -1){local_ret.clique3++;}
                        else{
                            loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                            if(loc1 != -1) {local_ret.clique5++;}
                        }
                    }
                    // Tri2_2_1 - Tri3_2
                    for (k2_idx = 0; k2_idx < local_tri3_2.count; ++k2_idx){
                        k2 = local_tri3_2.triend[k2_idx];
                        loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                        if(loc1 != -1) {local_ret.clique9++;}
                        else{local_ret.clique5++;}
                    }
                }

                // 3. Tri2_2_2-based
                for (k1_idx = 0; k1_idx < local_tri2_2_2.count; ++k1_idx){
                    k1 = local_tri2_2_2.triend[k1_idx];
                    for (k2_idx = k1_idx+1; k2_idx < local_tri2_2_2.count; ++k2_idx){
                        k2 = local_tri2_2_2.triend[k2_idx];
                        loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique7++;}
                        else{local_ret.clique4++;}
                    }
                    // Tri2_2_2 - Tri3_2
                    for (k2_idx = 0; k2_idx < local_tri3_2.count; ++k2_idx){
                        k2 = local_tri3_2.triend[k2_idx];
                        loc1 =  k1 < k2 ? gout->getEdgeBinary(k1, k2): gout->getEdgeBinary(k2, k1);
                        if(loc1 != -1) {local_ret.clique9++;}
                        else{local_ret.clique5++;}
                    }
                }

                // 4. Tri3_2-based
                for (k1_idx = 0; k1_idx < local_tri3_2.count; ++k1_idx){
                    k1 = local_tri3_2.triend[k1_idx];
                    for (k2_idx = k1_idx+1; k2_idx < local_tri3_2.count; ++k2_idx){
                        k2 = local_tri3_2.triend[k2_idx];
                        loc1 = gout->getEdgeBinary(k1, k2);
                        if(loc1 != -1) {local_ret.clique10++;}
                        else{local_ret.clique8++;}
                    }
                }  
            }

            for (e_j = gout_2->offsets[i]; e_j < gout_2->offsets[i+1]; ++e_j){
                j = gout_2->nbors[e_j];
                
                for (e_k = gout_2->offsets[j]; e_k < gout_2->offsets[j+1]; ++e_k) {
                    k = gout_2->nbors[e_k];
                    local_ret.cycle1 += (local_star1[k] * (local_star1[k] - 1)) / 2;
                    local_ret.cycle2 += (local_star1[k] * (local_star2_1[k] + local_star2_2[k]));
                    local_star1[k] = 0;
                }

                for (e_k = gin_2->offsets[j]; e_k < gin_2->offsets[j+1]; ++e_k) {
                    k = gin_2->nbors[e_k];
                    if (k > i) {
                        local_ret.cycle1 += (local_star1[k] * (local_star1[k] - 1)) / 2;
                        local_ret.cycle2 += (local_star1[k] * (local_star2_1[k] + local_star2_2[k]));
                    }
                    local_star1[k] = 0;
                }

                for (e_k = gout->offsets[j]; e_k < gout->offsets[j+1]; ++e_k) {
                    k = gout->nbors[e_k];
                    local_ret.cycle3 += (local_star2_1[k] * local_star2_2[k]);
                    local_star2_1[k] = 0;
                    local_star2_2[k] = 0;
                }

                for (e_k = gin->offsets[j]; e_k < gin->offsets[j+1]; ++e_k) {
                    k = gin->nbors[e_k];
                    if (k > i) {local_ret.cycle3 += (local_star2_1[k] * local_star2_2[k]);}
                    local_star2_1[k] = 0;
                    local_star2_2[k] = 0;
                }
            }
            
            delete[] local_star1;
            delete[] local_star2_1;
            delete[] local_star2_2;
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
            ret.cycle1 += local_ret.cycle1;
            ret.cycle2 += local_ret.cycle2;
            ret.cycle3 += local_ret.cycle3;
        }

    }

    Count chord1 = 0, chord2 = 0, chord3 = 0, chord4 = 0, chord5 = 0, chord6 = 0, chord7 = 0, chord8 = 0;
    Count e_tri1, e_tri2_1, e_tri2_2, e_tri3_1, e_tri3_2, e_tri4;
    for (EdgeIdx e2 = 0; e2 < gout_2->nEdges; ++e2){
        e_tri1 = 0, e_tri2_2 = 0, e_tri3_2 = 0;
        for (Count id = 0; id < num_threads; ++id){
            e_tri1 += all_local_tri1[id][e2];
            e_tri2_2 += all_local_tri2_2[id][e2];
            e_tri3_2 += all_local_tri3_2[id][e2];
        }
        chord1 += e_tri1 * (e_tri1 - 1) / 2;
        chord3 += e_tri2_2 * e_tri1;
        chord4 += e_tri2_2 * (e_tri2_2 - 1) / 2;
        chord6 += e_tri3_2 * e_tri1;
    }

    for (EdgeIdx e1 = 0; e1 < gout->nEdges; ++e1){
        e_tri2_1 = 0, e_tri3_1 = 0, e_tri4 = 0;
        for (Count id = 0; id < num_threads; ++id){
            e_tri2_1 += all_local_tri2_1[id][e1];
            e_tri3_1 += all_local_tri3_1[id][e1];
            e_tri4 += all_local_tri4[id][e1];
        }
        chord2 += e_tri2_1 * (e_tri2_1 - 1) / 2;
        chord5 += e_tri3_1 * e_tri2_1;
        chord7 += e_tri4 * e_tri2_1;
        chord8 += e_tri3_1 * (e_tri3_1 - 1) / 2;
    }

    ret.chord1 += chord1;
    ret.chord2 += chord2;
    ret.chord3 += chord3;
    ret.chord4 += chord4;
    ret.chord5 += chord5;
    ret.chord6 += chord6;
    ret.chord7 += chord7;
    ret.chord8 += chord8;

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
    mcounts[17] = motifcounts.chord1;
    mcounts[18] = motifcounts.chord2;
    mcounts[19] = motifcounts.chord3;
    mcounts[20] = motifcounts.chord4;
    mcounts[21] = motifcounts.chord5;
    mcounts[22] = motifcounts.chord6;
    mcounts[23] = motifcounts.chord7;
    mcounts[24] = motifcounts.chord8;
    
    mcounts[33] = motifcounts.cycle1;
    mcounts[34] = motifcounts.cycle2;
    mcounts[35] = motifcounts.cycle3;
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
    double d2_1 = mcounts[17] - 6 * mcounts[6] - mcounts[7];
    double d2_2 = mcounts[18] - mcounts[7] - 2 * mcounts[8];
    double d2_3 = mcounts[19] - 4 * mcounts[7] - 2 * mcounts[9];
    double d2_4 = mcounts[20] - 4 * mcounts[8] - mcounts[10] - mcounts[9] - 3 * mcounts[12];
    double d2_5 = mcounts[21] - 2 * mcounts[9] - 2 * mcounts[10];
    double d2_8 = mcounts[24] - mcounts[10] - 4 * mcounts[13] - 3 * mcounts[11] - mcounts[14];
    printf("d2-1 : %.1f\n", d2_1);
    printf("d2-2 : %.1f\n", d2_2);
    printf("d2-3 : %.1f\n", d2_3);
    printf("d2-4 : %.1f\n", d2_4);
    printf("d2-5 : %.1f\n", d2_5);
    printf("d2-6 : %.1f\n", mcounts[22] - mcounts[9] - 3 * mcounts[11]);
    printf("d2-7 : %.1f\n", mcounts[23] - 3 * mcounts[12] - mcounts[14]);
    printf("d2-8 : %.1f\n", d2_8);
    printf("d4-1 : %.1f\n", mcounts[33] - d2_1 - d2_2 - 3 * mcounts[6] - mcounts[7] - mcounts[8]);
    printf("d4-2 : %.1f\n", mcounts[34] - d2_3 - d2_5 - 2 * mcounts[7] - mcounts[9] - mcounts[10]);
    printf("d4-3 : %.1f\n", mcounts[35] - d2_4 - d2_8 - 2 * mcounts[8] - mcounts[10] - 2 * mcounts[13]);
}