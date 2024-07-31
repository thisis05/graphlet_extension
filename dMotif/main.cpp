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
    printf("Read Graph\n");
    read_mtx(filename, g);

    printf("Loaded graph\n");
    CGraph pre_cg = makeCSR(g);
    printf("Converted to CSR\n");

    printf("Relabeling graph\n");
    CGraph cg = pre_cg.renameByDegreeOrder();
    cg.sortById();

    printf("Creating DAG\n");
    CDAG dag = degreeOrdered(&cg);
    (dag.outlist).sortById();
    (dag.inlist).sortById();

    auto start_time1 = high_resolution_clock::now();

    printf("Get Edge set with distance 2\n");
    CGraph cg_2 = cg.getE2();
    cg_2.sortById();

    CDAG dag_2 = degreeOrdered2(&cg_2, &cg);
    (dag_2.outlist).sortById();
    (dag_2.inlist).sortById();

    auto end_time1 = high_resolution_clock::now();
    auto duration1 = duration_cast<milliseconds>(end_time1 - start_time1);
    double seconds1 = duration1.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for get E2:  %.3f\n", seconds1);
 
    // 1. Count 3-size d-Motifs
    double mcounts3[6];
    printf("Count d-Motifs (3-size)\n");
    auto start_time2 = high_resolution_clock::now();
    countThree(&cg, &dag, &cg_2, &dag_2, mcounts3);
    auto end_time2 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(end_time2 - start_time2);
    double seconds2 = duration2.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for Counting Motifs (3-size): %.3f\n", seconds2);
    mEquation3(mcounts3);

    // 2. Count 4-size d-Motifs
    double mcounts4[40];
    printf("Count d-Motifs (4-size)\n");
    auto start_time3 = high_resolution_clock::now();
    countFour(&cg, &dag, &cg_2, &dag_2, mcounts4);
    auto end_time3 = high_resolution_clock::now();
    auto duration3 = duration_cast<milliseconds>(end_time3 - start_time3);
    double seconds3 = duration3.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for Counting Motifs: %.3f\n", seconds3);
    mEquation4(mcounts4);

    printf("Total Execution time for 3-size: %.3f\n", seconds1 + seconds2);
    printf("Total Execution time for 4-size: %.3f\n", seconds1 + seconds3);

    printf("# of Edge (1): %lld", cg.nEdges / 2);
    printf("\n# of Edge (2): %lld", cg_2.nEdges / 2);
    printf("\nMax Degree (1): %lld", cg.maxDegree);
    printf("\nMax Degree (2): %lld\n", cg_2.maxDegree);

    #ifdef PRINT_CSR

    printf("CGraph ______________________________________________________\n");
    printf("Offsets ");
    for (VertexIdx i = 0; i <= cg.nVertices; ++i) {
        printf("%lld ", cg.offsets[i]);
    }
    printf("\n");


    printf("Destinations : ");
    for (EdgeIdx i = 0; i < cg.nEdges; ++i) {
        printf("%lld ", cg.nbors[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < cg.nVertices; ++i) {
        printf("%lld ", cg.degree(i));
    }
    printf("\n");
    
    printf("Offsets ");
    for (VertexIdx i = 0; i <= cg_2.nVertices; ++i) {
        printf("%lld ", cg_2.offsets[i]);
    }
    printf("\n");

    printf("Destinations : ");
    for (EdgeIdx i = 0; i < cg_2.nEdges; ++i) {
        printf("%lld ", cg_2.nbors[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < cg_2.nVertices; ++i) {
        printf("%lld ", cg_2.degree(i));
    }
    printf("\n");

    printf("DAG (Out) ______________________________________________________\n");
    printf("Offsets ");
    for (VertexIdx i = 0; i <= dag.outlist.nVertices; ++i) {
        printf("%lld ", dag.outlist.offsets[i]);
    }
    printf("\n");


    printf("Destinations : ");
    for (EdgeIdx i = 0; i < dag.outlist.nEdges; ++i) {
        printf("%lld ", dag.outlist.nbors[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < dag.outlist.nVertices; ++i) {
        printf("%lld ", dag.outlist.degree(i));
    }
    printf("\n");

    printf("Offsets ");
    for (VertexIdx i = 0; i <= dag_2.outlist.nVertices; ++i) {
        printf("%lld ", dag_2.outlist.offsets[i]);
    }
    printf("\n");

    printf("Destinations : ");
    for (EdgeIdx i = 0; i < dag_2.outlist.nEdges; ++i) {
        printf("%lld ", dag_2.outlist.nbors[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < dag_2.outlist.nVertices; ++i) {
        printf("%lld ", dag_2.outlist.degree(i));
    }
    printf("\n");
    #endif
    
    return 0;
}