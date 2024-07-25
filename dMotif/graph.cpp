#include "graph.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>  
#include <map>
#include <vector>

Graph::Graph() : nVertices(0), nEdges(0), out_srcs(nullptr), out_dsts(nullptr), in_srcs(nullptr), in_dsts(nullptr), verticesDegree(nullptr), out_offsets(nullptr), in_offsets(nullptr), adjVector(nullptr) {}

Graph::~Graph() {
    delete[] out_srcs;
    delete[] out_dsts;
    delete[] in_srcs;
    delete[] in_dsts;
    delete[] verticesDegree;
    delete[] out_offsets;
    delete[] in_offsets;
    delete adjVector;
}

void makeCSR(Graph& graph){

    graph.out_offsets = new EdgeIdx[graph.nVertices + 1]();
    std::fill(graph.out_offsets, graph.out_offsets + graph.nVertices+1, 0);

    graph.in_offsets = new EdgeIdx[graph.nVertices + 1]();
    std::fill(graph.in_offsets, graph.in_offsets + graph.nVertices+1, 0);

    for (EdgeIdx i = 0; i < graph.nEdges; ++i){
        auto out_src = graph.out_srcs[i];
        auto in_src = graph.in_srcs[i];
        graph.out_offsets[out_src+1]++;
        graph.in_offsets[in_src+1]++; 
    }
    
    for (EdgeIdx i = 0; i < graph.nVertices+1; ++i){
        graph.out_offsets[i+1] += graph.out_offsets[i]; 
        graph.in_offsets[i+1] += graph.in_offsets[i]; 
    }
    

    // Already Sorted by VertexID (MTX FIle case)
    // If you want to input a some other files (such as .edges..), some manipulation will be needed. (Like a sortById)
}

bool compareDegree(Pair firstPair, Pair nextPair) {
    if (firstPair.second != nextPair.second)
        return firstPair.second < nextPair.second;
    return firstPair.first < nextPair.first;
}

Graph Graph::renameByDegreeOrder() const {
    Graph ret;
    ret.nVertices = nVertices;
    ret.nEdges = nEdges;
    ret.out_srcs = new VertexIdx[nEdges]();
    ret.out_dsts = new VertexIdx[nEdges]();
    ret.in_srcs = new VertexIdx[nEdges]();
    ret.in_dsts = new VertexIdx[nEdges]();
    ret.verticesDegree = new VertexIdx[nVertices]();
    std::fill(ret.out_dsts, ret.out_dsts + ret.nEdges, 0);
    std::fill(ret.in_dsts, ret.in_dsts + ret.nEdges, 0);

    Pair* deg_info = new Pair[nVertices];
    VertexIdx* mapping = new VertexIdx[nVertices];
    VertexIdx* inverse = new VertexIdx[nVertices];

    // Construct array of pairs, storing old vertex label and degree
    for (VertexIdx i = 0; i < nVertices; i++) {
        deg_info[i].first = i;
        deg_info[i].second = verticesDegree[i];
    }

    // Sort the pairs by degree (if degree is same, sort by old vertex label)
    std::sort(deg_info, deg_info + nVertices, compareDegree);

    for (VertexIdx i = 0; i < nVertices; i++) {
        mapping[deg_info[i].first] = i;
        inverse[i] = deg_info[i].first;
    }

    EdgeIdx current_out = 0;
    EdgeIdx current_in = 0;
    for (VertexIdx new_label = 0; new_label < nVertices; new_label++) {
        VertexIdx old_label = inverse[new_label];
        ret.verticesDegree[new_label] = verticesDegree[old_label];
        for (const auto& neighbor : (*adjVector)[old_label]) {
            VertexIdx new_nbr = mapping[neighbor];
            if (new_nbr <= new_label) {
                ret.in_dsts[current_in] = new_nbr;
                ret.in_srcs[current_in] = new_label;
                current_in++;
            }
            else{
                ret.out_dsts[current_out] = new_nbr;
                ret.out_srcs[current_out] = new_label;
                current_out++;
            }  
        }
    }

    printf("Make CSR...\n");
    makeCSR(ret);
    ret.sortById();
    return ret;
}

void Graph::sortById() const {
    for (VertexIdx i = 0; i < nVertices; i++) {
        EdgeIdx out_start = out_offsets[i];
        EdgeIdx out_end = out_offsets[i+1];
        if (out_start < out_end) { 
            std::sort(out_dsts + out_start, out_dsts + out_end);
        }

        EdgeIdx in_start = in_offsets[i];
        EdgeIdx in_end = in_offsets[i+1];
        if (in_start < in_end) { 
            std::sort(in_dsts + in_start, in_dsts + in_end);
        }
    }
}


EdgeIdx Graph::getEdgeBinary(VertexIdx v1, VertexIdx v2) const
{
    if (v1 >= nVertices)
        return -1;
    EdgeIdx low = out_offsets[v1];
    EdgeIdx high = out_offsets[v1+1]-1;
    EdgeIdx mid;

    while(low <= high)
    {
        mid = (low+high)/2;
        if (out_dsts[mid] == v2)
            return mid;
        if (out_dsts[mid] > v2)
            high = mid-1;
        else
            low = mid+1;
    }
    return -1;
}


void read_mtx(const std::string& filename, Graph& graph) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (getline(file, line)) {
        if (line[0] != '%') break;
    }

    std::istringstream iss(line);
    int64_t n, m;
    iss >> n >> n >> m;
    graph.nVertices = n;
    graph.nEdges = m;
    graph.verticesDegree = new VertexIdx[graph.nVertices]();
    std::fill(graph.verticesDegree, graph.verticesDegree + graph.nVertices, 0);

    int64_t u, v;
    graph.adjVector = new std::map<VertexIdx, std::vector<VertexIdx>>();

    while (file >> u >> v) {
        u--, v--;
        (*graph.adjVector)[u].push_back(v);
        (*graph.adjVector)[v].push_back(u);
        graph.verticesDegree[u]++;
        graph.verticesDegree[v]++;
    }
    file.close();
}