#include "graph.h"
#include "counter.h"
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace std::chrono;

//#define PRINT_CSR

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

   
    auto start_time = high_resolution_clock::now();

    // 1. Count 3-size d-Motifs
    double mcounts[6];
    countThree(&dag, mcounts);

    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    double seconds = duration.count() / 1000.0; // Convert milliseconds to seconds
    std::cout << fixed << setprecision(3) << "Execution time: " << seconds << " seconds" << endl;

    #ifdef PRINT_CSR
    printf("\nEdge 1 ------------------------------------------\n\n");
    printf("Offsets (Out): ");
    for (VertexIdx i = 0; i <= dag.nVertices; ++i) {
        printf("%lld ", dag.out_offsets[i]);
    }
    printf("\n");

    printf("Sources (Out): ");
    for (EdgeIdx i = 0; i < dag.nEdges; ++i) {
        printf("%lld ", dag.out_srcs[i]);
    }
    printf("\n");

    printf("Destinations (Out): ");
    for (EdgeIdx i = 0; i < dag.nEdges; ++i) {
        printf("%lld ", dag.out_dsts[i]);
    }
    printf("\n");

        printf("Offsets (In): ");
    for (VertexIdx i = 0; i <= dag.nVertices; ++i) {
        printf("%lld ", dag.in_offsets[i]);
    }
    printf("\n");

    printf("Sources (In): ");
    for (EdgeIdx i = 0; i < dag.nEdges; ++i) {
        printf("%lld ", dag.in_srcs[i]);
    }
    printf("\n");

    printf("Destinations (In): ");
    for (EdgeIdx i = 0; i < dag.nEdges; ++i) {
        printf("%lld ", dag.in_dsts[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < dag.nVertices; ++i) {
        printf("%lld ", dag.verticesDegree[i]);
    }
    printf("\n");

    printf("\nEdge 2 ------------------------------------------\n\n");
    printf("Offsets (Out 2): ");
    for (VertexIdx i = 0; i <= dag.nVertices; ++i) {
        printf("%lld ", dag.out_offsets2[i]);
    }
    printf("\n");

    printf("Sources (Out 2): ");
    for (EdgeIdx i = 0; i < dag.nEdges2; ++i) {
        printf("%lld ", dag.out_srcs2[i]);
    }
    printf("\n");

    printf("Destinations (Out 2): ");
    for (EdgeIdx i = 0; i < dag.nEdges2; ++i) {
        printf("%lld ", dag.out_dsts2[i]);
    }
    printf("\n");
    
    printf("Offsets (In 2): ");
    for (VertexIdx i = 0; i <= dag.nVertices; ++i) {
        printf("%lld ", dag.in_offsets2[i]);
    }
    printf("\n");

    printf("Sources (In 2): ");
    for (EdgeIdx i = 0; i < dag.nEdges2; ++i) {
        printf("%lld ", dag.in_srcs2[i]);
    }
    printf("\n");

    printf("Destinations (In 2): ");
    for (EdgeIdx i = 0; i < dag.nEdges2; ++i) {
        printf("%lld ", dag.in_dsts2[i]);
    }
    printf("\n");

    printf("Degree 2: ");
    for (VertexIdx i = 0; i < dag.nVertices; ++i) {
        printf("%lld ", dag.verticesDegree2[i]);
    }
    #endif
    
    return 0;
}