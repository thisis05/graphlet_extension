#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <string>
#include <chrono>
#include <omp.h>
#include <limits>
#include <iomanip>
#include <atomic>
#include <algorithm>

#define PRINT_CSR

using namespace std;
using namespace std::chrono;

using VertexIdx = int64_t;
using EdgeIdx = int64_t;
using Count = int64_t;

struct Graph {
    VertexIdx nVertices; // number of vertices in the graph
    EdgeIdx nEdges; // number of edges in this list
    VertexIdx* srcs; // array of source vertices
    VertexIdx* dsts; // array of destination vertices
    VertexIdx* verticesDegree;
    EdgeIdx* offsets; // vertex v's edges are [offset[v], offset[v + 1])
    map<VertexIdx, vector<VertexIdx>>* adjVector;

    Graph() : nVertices(0), nEdges(0), srcs(nullptr), dsts(nullptr), verticesDegree(nullptr), offsets(nullptr), adjVector(nullptr) {}
    ~Graph() {
        delete[] srcs;
        delete[] dsts;
        delete[] verticesDegree;
        delete[] offsets;
        delete adjVector;
    }
    Graph renameByDegreeOrder() const; // sort by degree of vertex. for all i < j, the degree of i is less than that of j.
    void sortById() const;
};

struct Pair {
    VertexIdx first;
    VertexIdx second;
};

bool compareDegree(Pair firstPair, Pair nextPair) {
    if (firstPair.second != nextPair.second)
        return firstPair.second < nextPair.second;
    return firstPair.first < nextPair.first;
}

void makeCSR(Graph& graph){

    graph.offsets = new EdgeIdx[graph.nVertices + 1]();
    fill(graph.offsets, graph.offsets + graph.nVertices+1, 0);
    VertexIdx cv = 0; //currenct vertex id

    for (EdgeIdx i = 0; i < graph.nEdges; ++i){
        auto src = graph.srcs[i];
        graph.offsets[src+1]++; 
    }
    for (EdgeIdx i = 0; i < graph.nVertices+1; ++i){
        graph.offsets[i+1] += graph.offsets[i]; 
    }

    // Already Sorted by VertexID (MTX FIle case)
    // If you want to input a some other files (such as .edges..), some manipulation will be needed. (Like a sortById)
}

Graph Graph::renameByDegreeOrder() const {
    Graph ret;
    ret.nVertices = nVertices;
    ret.nEdges = nEdges;
    ret.srcs = new VertexIdx[nEdges]();
    ret.dsts = new VertexIdx[nEdges]();
    ret.verticesDegree = new VertexIdx[nVertices]();
    fill(ret.dsts, ret.dsts + nEdges, 0);

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

    VertexIdx current = 0;
    for (VertexIdx new_label = 0; new_label < nVertices; new_label++) {
        VertexIdx old_label = inverse[new_label];
        ret.verticesDegree[new_label] = verticesDegree[old_label];
        for (const auto& neighbor : (*adjVector)[old_label]) {
            VertexIdx new_nbr = mapping[neighbor];
            if (new_nbr <= new_label) {continue;}
            ret.dsts[current] = new_nbr;
            ret.srcs[current] = new_label;
            current++;
        }
    }
    ret.nEdges = current;

    printf("Make CSR...\n");
    makeCSR(ret);
    ret.sortById();
    return ret;
}

void Graph::sortById() const {
    for (VertexIdx i = 0; i < nVertices; i++) {
        EdgeIdx start = offsets[i];
        EdgeIdx end = offsets[i+1];
        if (start < end) { 
            std::sort(dsts + start, dsts + end);
        }
    }
}



void read_mtx(const string& filename, Graph& graph) {
    ifstream file(filename);
    if (!file) {
        cerr << "Failed to open file: " << filename << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        if (line[0] != '%') break;
    }

    istringstream iss(line);
    int64_t n, m;
    iss >> n >> n >> m;
    graph.nVertices = n;
    graph.nEdges = m;
    graph.verticesDegree = new VertexIdx[graph.nVertices]();
    fill(graph.verticesDegree, graph.verticesDegree + graph.nVertices, 0);

    int64_t u, v;
    graph.adjVector = new map<VertexIdx, vector<VertexIdx>>();

    while (file >> u >> v) {
        u--, v--;
        (*graph.adjVector)[u].push_back(v);
        (*graph.adjVector)[v].push_back(u);
        graph.verticesDegree[u]++;
        graph.verticesDegree[v]++;
    }
    file.close();
}


int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }
    string filename = argv[1];
    Graph g;

    printf("Read Graph...\n");
    read_mtx(filename, g);

    printf("Convert to DAG...\n");
    Graph dag = g.renameByDegreeOrder();

    #ifdef PRINT_CSR
    printf("Offsets: ");
    for (VertexIdx i = 0; i <= dag.nVertices; ++i) {
        printf("%lld ", dag.offsets[i]);
    }
    printf("\n");

    printf("Sources: ");
    for (EdgeIdx i = 0; i < dag.nEdges; ++i) {
        printf("%lld ", dag.srcs[i]);
    }
    printf("\n");

    printf("Destinations: ");
    for (EdgeIdx i = 0; i < dag.nEdges; ++i) {
        printf("%lld ", dag.dsts[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < dag.nVertices; ++i) {
        printf("%lld ", dag.verticesDegree[i]);
    }
    printf("\n");
    #endif
    
    
    return 0;
}