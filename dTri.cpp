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

tuple<VertexSet, EdgeSet, EdgeSet, unordered_map<int, unordered_set<int>>, unordered_map<int, unordered_set<int>>> readAndSet(const string& filename) {
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
    unordered_map<int, unordered_set<int>> adjVector;
    unordered_map<int, unordered_set<int>> adjVector2;

    int u, v;
    VertexSet vertexSet(n);
    EdgeSet edgeSet1;
    EdgeSet edgeSet2;

    while (file >> u >> v) {
        u--, v--;
        adjVector[u].insert(v);
        adjVector[v].insert(u);
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
                    adjVector2[start_node].insert(final_node);
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

void getTri(VertexSet& vset, Edge& e, unordered_map<int, unordered_set<int>>& adj1, unordered_map<int, unordered_set<int>>& adj2,
            unordered_set<int>& star_u1, unordered_set<int>& star_u2, unordered_set<int>& star_v1, unordered_set<int>& star_v2, 
            vector<int>& tri1, vector<int>& tri2, 
            unordered_map<int, int>& X1, unordered_map<int, int>& X2,
            long long& Tri1, long long& Tri2, long long& Path1, long long& Path2, long long& Path3, 
            const int& type) {

    // type 1: edge 1 -> distance 1 : Tri4 / edge 1 -> distance 2 : Tri2
    // type 2: edge 2 -> distance 1 : Tri3 / edge 2 -> distance 2 : Tri1

    int u = e.src;
    int v = e.dest;

    for (auto w : adj1[u]) {
        if (w == v) continue;
        X1[w] = 1;
        star_u1.insert(w);
    }

    for (auto w : adj2[u]) {
        if (w == v) continue;
        X2[w] = 1;
        star_u2.insert(w);
    }

    for (auto w : adj1[v]) {
        if (w == u) continue;
        if (X2[w] == 1){
            star_u2.erase(w);
            continue;
        }
        if (X1[w] == 1) {
            Tri1++;
            X1[w] = 2;
            tri1.push_back(w);
            star_u1.erase(w);
        }
        else {
            star_v1.insert(w);
            X1[w] = 3;     
        }
    }
    for (auto w : adj2[v]) {
        if (w == u) continue;
        if (X1[w] == 1){
            star_u1.erase(w);
            continue;
        }
        if (X2[w] == 1) {
            X2[w] = 2;
            Tri2++;
            tri2.push_back(w);
            star_u2.erase(w);
        }
        else {
            star_v2.insert(w);
            X2[w] = 3;
        }
    }

    if (type == 1){
        vset.setTri4(u, Tri1);
        vset.setTri4(v, Tri1);
        e.setTri4(Tri1);

        vset.setTri2(u, Tri2);
        vset.setTri2(v, Tri2);
        e.setTri2(Tri2);

        Path1 = star_u2.size() * star_v2.size();
    }
    else {
        vset.setTri3(u, Tri1);
        vset.setTri3(v, Tri1);
        e.setTri3(Tri1);

        vset.setTri1(u, Tri2);
        vset.setTri1(v, Tri2);
        e.setTri1(Tri2);
        
        Path1 = star_u2.size() * star_v2.size();
        Path2 = (star_u2.size() * star_v1.size()) + (star_u1.size() * star_v2.size());
        Path3 = star_u1.size() * star_v1.size();
    }
}


void getCycle(unordered_map<int, int> X, unordered_map<int, int> X2, unordered_set<int>& star, unordered_map<int, unordered_set<int>>& adj, long long& cycle, Edge& edge, int type){
    for (int w : star){
        for (int r : adj[w]){
            if (X[r] == 3) {
                cycle +=1;
            }
        }
        X[w] = 0;
    }

}

void getClique(unordered_map<int, int> X, vector<int>& tri, unordered_map<int, unordered_set<int>>& adj, long long& clique){
    for (auto w : tri){
        for (auto r : adj[w]){
            if (X[r] == 2)
                clique +=1;
        }
        X[w] = 0;
    }

}


map<string, long long> countMotifs(VertexSet& vset, EdgeSet& e1, EdgeSet& e2, unordered_map<int, unordered_set<int>>& adj1, unordered_map<int, unordered_set<int>>& adj2) {
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
            long long Path3 = 0; 
            long long Star = 0;
            unordered_set<int> star_u1; unordered_set<int> star_u2; unordered_set<int> star_v1; unordered_set<int> star_v2;
            vector<int> tri4_set; vector<int> tri2_set;
            unordered_map<int, int> X1; unordered_map<int, int> X2;
            unordered_map<int, bool> star_check;
            getTri(vset, edge, adj1, adj2, star_u1, star_u2, star_v1, star_v2, tri4_set, tri2_set, X1, X2, Tri4, Tri2, Path3, Path3, Path3, 1);
            long long Cycle2 = 0; long long Cycle3 = 0;
            long long Clique1 = 0;
            getCycle(X2, X1, star_u2, adj1, Cycle3, edge, 3);
            getCycle(X2, X1, star_u2, adj2, Cycle2, edge, 2);
            thread_motifCounts[thread_id]["Tri4"] += Tri4;
            thread_motifCounts[thread_id]["Tri2"] += Tri2;
            thread_motifCounts[thread_id]["4-6-3"] += Path3;
            thread_motifCounts[thread_id]["4-4-2"] += Cycle2;
            thread_motifCounts[thread_id]["4-4-3"] += Cycle3;
            thread_motifCounts[thread_id]["4-5-2"] += Star;

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
            long long Star = 0;
            unordered_set<int> star_u1; unordered_set<int> star_u2; unordered_set<int> star_v1; unordered_set<int> star_v2;
            vector<int> tri3_set; vector<int> tri1_set;
            unordered_map<int, int> X1; unordered_map<int, int> X2;
            unordered_map<int, bool> star_check;
            getTri(vset, edge, adj1, adj2, star_u1, star_u2, star_v1, star_v2, tri3_set, tri1_set, X1, X2, Tri3, Tri1, Path1, Path2, Path4, 2);
            long long Cycle1 = 0;
            long long Clique1 = 0; long long Clique2 = 0;
            getCycle(X2, X1, star_u2, adj2, Cycle1, edge, 1);
            getClique(X2, tri1_set, adj2, Clique1);
            getClique(X2, tri1_set, adj1, Clique2);
            thread_motifCounts[thread_id]["Tri3"] += Tri3;
            thread_motifCounts[thread_id]["Tri1"] += Tri1;
            thread_motifCounts[thread_id]["4-6-1"] += Path1;
            thread_motifCounts[thread_id]["4-6-2"] += Path2;
            thread_motifCounts[thread_id]["4-6-4"] += Path4;
            thread_motifCounts[thread_id]["4-4-1"] += Cycle1;
            thread_motifCounts[thread_id]["4-1-1"] += Clique1;
            thread_motifCounts[thread_id]["4-1-2"] += Clique2;

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
    tuple<VertexSet, EdgeSet, EdgeSet, unordered_map<int, unordered_set<int>>, unordered_map<int, unordered_set<int>>> result = readAndSet(filename);

    VertexSet vset = get<0>(result);
    EdgeSet eset1 = get<1>(result);
    EdgeSet eset2 = get<2>(result);
    unordered_map<int, unordered_set<int>> adj1 = get<3>(result);
    unordered_map<int, unordered_set<int>> adj2 = get<4>(result);

    map <string, long long> results = countMotifs(vset, eset1, eset2, adj1, adj2);

    //graphlet equation

    results["Tri1"] /= 3; results["Tri4"] /= 3;
    results["4-4-1"] /= 4; results["4-4-3"] /= 2;
    results["4-1-1"] /= 6;
    results["4-6-1"] = results["4-6-1"] - 4*results["4-4-1"] - results["4-4-2"];
    results["4-6-2"] = results["4-6-2"] - 2*results["4-4-2"];
    results["4-6-3"] = results["4-6-3"] - 2*results["4-4-3"] - results["4-4-2"];
    results["4-6-4"] = results["4-6-4"] - 2*results["4-4-3"];

    
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
