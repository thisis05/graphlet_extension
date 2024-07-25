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
    VertexIdx* out_srcs;
    VertexIdx* out_dsts;
    VertexIdx* in_srcs;
    VertexIdx* in_dsts;
    VertexIdx* verticesDegree;
    EdgeIdx* out_offsets;
    EdgeIdx* in_offsets;
    std::map<VertexIdx, std::vector<VertexIdx>>* adjVector;

    Graph();
    ~Graph();
    Graph renameByDegreeOrder() const;
    void sortById() const;
    EdgeIdx getEdgeBinary(VertexIdx v1, VertexIdx v2) const;
};

struct Pair {
    VertexIdx first;
    VertexIdx second;
};

bool compareDegree(Pair firstPair, Pair nextPair);
void makeCSR(Graph& graph);
void read_mtx(const std::string& filename, Graph& graph);