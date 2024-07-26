#pragma once

#include <cstdint>
#include <map>
#include <vector>
#include <string>

using VertexIdx = int64_t;
using EdgeIdx = int64_t;

struct Graph {
    VertexIdx nVertices;
    EdgeIdx nEdges;
    EdgeIdx nEdges2;

    // distance 1 case
    VertexIdx* out_srcs;
    VertexIdx* out_dsts;
    VertexIdx* in_srcs;
    VertexIdx* in_dsts;
    VertexIdx* verticesDegree;
    EdgeIdx* out_offsets;
    EdgeIdx* in_offsets;

    // distance 2 case
    VertexIdx* out_srcs2;
    VertexIdx* out_dsts2;
    VertexIdx* in_srcs2;
    VertexIdx* in_dsts2;
    VertexIdx* verticesDegree2;
    EdgeIdx* out_offsets2;
    EdgeIdx* in_offsets2;

    std::map<VertexIdx, std::vector<VertexIdx>>* adjVector;
    std::map<VertexIdx, std::vector<VertexIdx>>* adjVector2;

    Graph();
    ~Graph();

    Graph renameByDegreeOrder() const;
    void createE2attr();
    void sortById(int etype) const;
    EdgeIdx getEdgeBinary(VertexIdx v1, VertexIdx v2, VertexIdx* found) const;
};

struct Pair {
    VertexIdx first;
    VertexIdx second;
};

bool compareDegree(Pair firstPair, Pair nextPair);
void makeCSR(Graph& graph, int etype);
void read_mtx(const std::string& filename, Graph& graph);