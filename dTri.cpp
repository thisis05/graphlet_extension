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
#include "vertex.h"
#include "edge.h"

using namespace std;
using namespace std::chrono;

tuple<VertexSet, EdgeSet, EdgeSet, unordered_map<int, vector<int>>, unordered_map<int, vector<int>>> readAndSet(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Failed to open file: " << filename << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        if (line[0] != '%') break;
    }

    int n, m;
    istringstream iss(line);
    iss >> n >> n >> m;

    printf("Read File ... \n");
    unordered_map<int, vector<int>> adjVector;
    unordered_map<int, vector<int>> adjVector2;

    int u, v;
    VertexSet vertexSet(n);
    EdgeSet edgeSet1;
    EdgeSet edgeSet2;

    while (file >> u >> v) {
        u--, v--;
        adjVector[u].push_back(v);
        adjVector[v].push_back(u);
    }
    file.close();

    // 정점의 차수를 adjVector의 크기로 설정
    for (int i = 0; i < n; i++) {
        vertexSet.addVertex(i);
        vertexSet.setDegree1(i, adjVector[i].size());
    }

    for (int u = 0; u < n; u++) {
        for (int v : adjVector[u]) {
            if (u < v)
                edgeSet1.addEdge(u, v, vertexSet.getDegree1(u), vertexSet.getDegree1(v));
        }
    }
    edgeSet1.sortEdgesByDegree();

    printf("Calculate distance ... \n");
    for (int start_node = 0; start_node < n; start_node++) {
        unordered_set<int> alreadyNeighbors(adjVector[start_node].begin(), adjVector[start_node].end());
        for (const auto& second_node : adjVector[start_node]) {
            for (const auto& final_node : adjVector[second_node]) {
                if (start_node != final_node && alreadyNeighbors.find(final_node) == alreadyNeighbors.end()) {
                    adjVector2[start_node].push_back(final_node);
                }
            }
        }
        vertexSet.setDegree2(start_node, adjVector2[start_node].size());
    }

    for (int u = 0; u < n; u++) {
        for (int v : adjVector2[u]) {
            if (u < v)
                edgeSet2.addEdge(u, v, vertexSet.getDegree2(u), vertexSet.getDegree2(v));
        }
    }
    edgeSet2.sortEdgesByDegree();

    return make_tuple(vertexSet, edgeSet1, edgeSet2, adjVector, adjVector2);
}

void getTri1(VertexSet& vset, Edge& e, unordered_map<int, vector<int>>& adj1, unordered_map<int, vector<int>>& adj2, 
             long long& Tri4, long long& Tri2, long long& Path3) {
    int u = e.src;
    int v = e.dest;
    
    int deg2_u = 0; int deg2_v= 0;

    unordered_map<int, bool> X1;
    unordered_map<int, bool> X2;

    for (auto w : adj1[u]) {
        if (w == v) continue;
        X1[w] = true;
    }
    for (auto w : adj1[v]) {
        if (w == u) continue;
        if (X1[w] == true) {
            Tri4++;
        }
    }

    for (auto w : adj2[u]) {
        if (w == v) continue;
        X2[w] = true;
    }
    for (auto w : adj2[v]) {
        if (w == u) continue;
        if (X2[w] == true) {
            Tri2++;
        }
    }

    vset.setTri4(u, Tri4);
    vset.setTri4(v, Tri4);
    e.setTri4(Tri4);

    vset.setTri2(u, Tri2);
    vset.setTri2(v, Tri2);
    e.setTri2(Tri2);

    deg2_u = vset.getDegree2(u);
    deg2_v = vset.getDegree2(v);

    //4-6-3
    Path3 = (deg2_u-1)*(deg2_v -1);


}

void getTri2(VertexSet& vset, Edge& e, unordered_map<int, vector<int>>& adj1, unordered_map<int, vector<int>>& adj2,
             vector<int>& star, vector<int>& tri1, vector<int>& tri2, unordered_map<int, int>& X1, unordered_map<int, int>& X2,
             long long& Tri3, long long& Tri1, long long& Path1, long long& Path2, long long& Path4) {
    int u = e.src;
    int v = e.dest;
    
    int deg1_u = 0; int deg1_v= 0; int deg2_u = 0; int deg2_v= 0;


    for (auto w : adj1[u]) {
        if (w == v) continue;
        X1[w] = 1;

    }
    for (auto w : adj1[v]) {
        if (w == u) continue;
        if (X1[w] == 1) {
            Tri3++;
            X1[w] = 2;
            tri1.push_back(w);
        }
        else{
            X1[w] = 3;
        }
    }


    for (auto w : adj2[u]) {
        if (w == v) continue;
        X2[w] = 1;
    }
    for (auto w : adj2[v]) {
        if (w == u) continue;
        if (X2[w] == 1) {
            Tri1++;
            X2[w] =2;
            tri2.push_back(w);
        }
        else{
            X2[w] = 3;
            star.push_back(w);
        }
    }

    vset.setTri3(u, Tri3);
    vset.setTri3(v, Tri3);
    e.setTri3(Tri3);

    vset.setTri1(u, Tri1);
    vset.setTri1(v, Tri1);
    e.setTri1(Tri1);


    deg1_u = vset.getDegree1(u); deg2_u = vset.getDegree2(u);
    deg1_v = vset.getDegree1(v); deg2_v = vset.getDegree2(v);

    //4-6-1, 4-6-2, 4-6-4
    Path1 = (deg2_u-1)*(deg2_v -1);
    Path2 = (deg2_u-1)*(deg1_v -1) + (deg1_u-1)*(deg2_v -1);
    Path4 = (deg1_u-1)*(deg1_v -1);

}

void getCycle(unordered_map<int, int>& X, vector<int>& star, unordered_map<int, vector<int>>& adj, long long& cycle){
    for (auto w : star){
        for (auto r : adj[w]){
            if (X[r] == 3)
                cycle +=1;
        }
    }

}


map<string, long long> countMotifs(VertexSet& vset, EdgeSet& e1, EdgeSet& e2, unordered_map<int, vector<int>>& adj1, unordered_map<int, vector<int>>& adj2) {
    map<string, long long> motifCounts = {
        {"Tri1", 0}, {"Tri2", 0}, {"Tri3", 0}, {"Tri4", 0},
        //clique
        {"4-1-1", 0}, {"4-1-2", 0}, {"4-1-3", 0}, {"4-1-4", 0}, {"4-1-5", 0}, {"4-1-6", 0}, {"4-1-7", 0}, {"4-1-8", 0}, {"4-1-9", 0}, {"4-1-10", 0},
        //chordalcycle
        {"4-2-1", 0}, {"4-2-2", 0}, {"4-2-3", 0}, {"4-2-4", 0}, {"4-2-5", 0}, {"4-2-6", 0}, {"4-2-7", 0},
        //tailedtriangle
        {"4-3-1", 0}, {"4-3-2", 0}, {"4-3-3", 0}, {"4-3-4", 0}, {"4-3-5", 0}, {"4-3-6", 0},
        //cycle
        {"4-4-1", 0}, {"4-4-2", 0}, {"4-4-3", 0},
        //star
        {"4-5-1", 0}, {"4-5-2", 0},
        //path
        {"4-6-1", 0}, {"4-6-2", 0}, {"4-6-3", 0}, {"4-6-4", 0}
    };

    int max_num_workers = omp_get_max_threads() - 1;
    omp_set_num_threads(max_num_workers);
    int total_e1 = e1.size();
    int total_e2 = e2.size();
    printf("Get Max Thread : %d\n", max_num_workers);
    vector<map<string, long long>> thread_motifCounts(max_num_workers, motifCounts);

    int progress1 = 0;
    int next_progress1_update = total_e1 / 10;

    int progress2 = 0;
    int next_progress2_update = total_e2 / 10;
    // start parallel
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
    
        #pragma omp for schedule(dynamic)
        for (Edge edge : e1.getEdges()) {
            long long Tri4 =0; long long Tri2 = 0; 
            long long Path3=0;
            getTri1(vset, edge, adj1, adj2, Tri4, Tri2, Path3);
            thread_motifCounts[thread_id]["Tri4"] += Tri4;
            thread_motifCounts[thread_id]["Tri2"] += Tri2;
            thread_motifCounts[thread_id]["4-6-3"] += Path3;

            #pragma omp atomic
            progress1++;
            if (progress1 == next_progress1_update) {
                #pragma omp critical
                {
                    if (progress1 >= next_progress1_update) {
                        printf("Edge 1 Progress: %d / %d (%d%%)\n", progress1, total_e1, (progress1 * 100) / total_e1);
                        next_progress1_update += total_e1 / 10;  
                    }
                }
            }
        }

        #pragma omp for schedule(dynamic)
        for (Edge edge : e2.getEdges()) {
            long long Tri3 =0; long long Tri1 = 0; 
            long long Path1=0; long long Path2=0; long long Path4=0;
            vector<int> star; vector<int> tri1; vector<int> tri2;
            unordered_map<int, int> X1; unordered_map<int, int> X2;
            getTri2(vset, edge, adj1, adj2, star, tri1, tri2, X1, X2, Tri3, Tri1, Path1, Path2, Path4);
            long long Cycle1 = 0;
            getCycle(X2, star, adj2, Cycle1);
            thread_motifCounts[thread_id]["Tri3"] += Tri3;
            thread_motifCounts[thread_id]["Tri1"] += Tri1;
            thread_motifCounts[thread_id]["4-6-1"] += Path1;
            thread_motifCounts[thread_id]["4-6-2"] += Path2;
            thread_motifCounts[thread_id]["4-6-4"] += Path4;
            thread_motifCounts[thread_id]["4-4-1"] += Cycle1;

            #pragma omp atomic
            progress2++;
            if (progress2 == next_progress2_update) {
                #pragma omp critical
                {
                    if (progress2 >= next_progress2_update) {
                        printf("Edge 2 Progress: %d / %d (%d%%)\n", progress2, total_e2, (progress2 * 100) / total_e2);
                        next_progress2_update += total_e2 / 10; 
                    }
                }
            }
        }
    } // end parallel

    for (int i = 0; i < max_num_workers; i++) {
        for (const auto& motif : thread_motifCounts[i]) {
            motifCounts[motif.first] += motif.second;
        }
    }

    return motifCounts;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }
    auto start_time = high_resolution_clock::now();

    string filename = argv[1];
    tuple<VertexSet, EdgeSet, EdgeSet, unordered_map<int, vector<int>>, unordered_map<int, vector<int>>> result = readAndSet(filename);

    VertexSet vset = get<0>(result);
    EdgeSet eset1 = get<1>(result);
    EdgeSet eset2 = get<2>(result);
    unordered_map<int, vector<int>> adj1 = get<3>(result);
    unordered_map<int, vector<int>> adj2 = get<4>(result);

    map <string, long long> results = countMotifs(vset, eset1, eset2, adj1, adj2);

    results["Tri1"] /= 3; results["Tri2"] /= 3; results["Tri3"] /= 3; results["Tri4"] /= 3;
    results["4-4-1"] /= 4; results["4-4-3"] /= 2;
    results["4-6-1"] -= 3*results["Tri1"]; results["4-6-2"] -= 6*results["Tri2"];
    results["4-6-3"] -= 3*results["Tri2"]; results["4-6-4"] -= 3*results["Tri3"];

    
    cout << "Motif Counts:" << endl;
    for (const auto& motif : results) {
        cout << "\"" << motif.first << "\"" << " : " << motif.second << "," << endl;
    }
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    double seconds = duration.count() / 1000.0; // Convert milliseconds to seconds

    cout << fixed << setprecision(3) << "Execution time: " << seconds << " seconds" << endl;

    // // 정점들 출력
    // vset.print();

    // // edgeSet1 출력
    // std::cout << "EdgeSet1:" << std::endl;
    // eset1.print();

    // // edgeSet2 출력
    // std::cout << "EdgeSet2:" << std::endl;
    // eset2.print();

    // 간선들의 개수 출력
    std::cout << "Number of edges in EdgeSet1: " << eset1.size() << std::endl;
    std::cout << "Number of edges in EdgeSet2: " << eset2.size() << std::endl;


    return 0;
}
