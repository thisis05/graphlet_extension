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
#include <filesystem>
#include <iomanip>
#include <atomic> 
#include "vertex.h"
#include "edge.h"

namespace fs = std::filesystem;
using namespace std;
using namespace std::chrono;
int max_num_workers = omp_get_max_threads()-1;

tuple<VertexSet, EdgeSet, EdgeSet, vector<vector<int>>, vector<vector<int>>> readAndSet(const string& filename) {
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
    vector<vector<int>> adjVector(n);
    vector<vector<int>> adjVector2(n);

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
    //edgeSet1.sortEdgesByDegree();

    printf("Calculate distance ... \n");
    
    atomic<int> completedNodes(0);

    #pragma omp parallel for
    for (int start_node = 0; start_node < n; start_node++) {
        unordered_set<int> alreadyNeighbors;
        for (const auto& second_node : adjVector[start_node]) {
            for (const auto& final_node : adjVector[second_node]) {
                if (start_node != final_node && !binary_search(adjVector[start_node].begin(), adjVector[start_node].end(), final_node)) {
                    #pragma omp critical
                    {
                        if (alreadyNeighbors.insert(final_node).second)
                            adjVector2[start_node].push_back(final_node);
                    }
                }
            }
        }
        vertexSet.setDegree2(start_node, adjVector2[start_node].size());

        // 진행 상황 출력
        int completed = completedNodes.fetch_add(1) + 1;
        if (completed % (n / 10) == 0 || completed == n) {
            int thread_id = omp_get_thread_num();
            int num_threads = omp_get_num_threads();
            #pragma omp critical
            {
                cout << "Completed Nodes : " << completed << " / " << n << endl;
            }
        }
    }


    for (int u = 0; u < n; u++) {
        for (int v : adjVector2[u]) {
            if (u < v){
                edgeSet2.addEdge(u, v, vertexSet.getDegree2(u), vertexSet.getDegree2(v));
                //printf("%d %d\n", u, v);
            }
        }
    }
    //edgeSet2.sortEdgesByDegree();

    return make_tuple(vertexSet, edgeSet1, edgeSet2, adjVector, adjVector2);
}

void removeVec(vector<int>& vec, int value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end()) {
        std::iter_swap(it, vec.end() - 1); // 삭제할 원소와 마지막 원소를 교체
        vec.pop_back(); // 마지막 원소를 삭제
    }
}

void getTri1(VertexSet& vset, Edge& e, vector<vector<int>>& adj1, vector<vector<int>>& adj2,
            unordered_set<int>& star2_u, unordered_set<int>& star2_v, 
            unordered_set<int>& tri2, unordered_set<int>& tri3_1, unordered_set<int>& tri3_2, unordered_set<int>& tri4, 
            vector<int>& X1, vector<int>& X2,
            long long& Tri2,  long long& Tri4,
            long long& Path3,
            long long& Star2,
            long long& Tailed3, long long& Tailed6, long long& Tailed8,
            long long& Chord2, long long& Chord5, long long& Chord7, long long& Chord8,
            long long& Clique10) {

    // edge : original edge
    // star3_u : temp for original star version of motif

    int u = e.src;
    int v = e.dest;

    for (auto w : adj1[u]) {
        if (w == v) continue;
        X1[w] = 1;
        //star3_u.insert(w);
    }

    for (auto w : adj2[u]) {
        if (w == v) continue;
        X2[w] = 1;
        star2_u.insert(w);
    }

    for (auto w : adj1[v]) {
        if (w == u) continue;
        if (X2[w] == 1){
            star2_u.erase(w);
            tri3_1.insert(w);
            X1[w] = 3;
            continue;
        }
        if (X1[w] == 1) {
            Tri4++;
            X1[w] = 4;
            tri4.insert(w);
            //star3_u.erase(w);
        }
        else {
            //star3_v.insert(w);
            X1[w] = 2;     
        }
    }
    for (auto w : adj2[v]) {
        if (w == u) continue;
        if (X1[w] == 1){
            //star3_u.erase(w);
            tri3_2.insert(w);
            X2[w] = 3;
            continue;
        }
        if (X2[w] == 1) {
            X2[w] = 4;
            Tri2++;
            tri2.insert(w);
            star2_u.erase(w);
        }
        else {
            star2_v.insert(w);
            X2[w] = 2;
        }
    }

    Path3 = star2_u.size() * star2_v.size();
    Star2 = (star2_u.size() * (star2_u.size()-1) / 2) + (star2_v.size() * (star2_v.size() - 1) / 2);
    Tailed3 = (star2_u.size() + star2_v.size()) * tri2.size();
    Tailed6 = (star2_u.size() * tri3_2.size()) + (star2_v.size() * tri3_1.size());
    Tailed8 = (star2_u.size() + star2_v.size()) * tri4.size();
    Chord2 = tri2.size() * (tri2.size()-1) / 2;
    Chord5 = (tri3_1.size() + tri3_2.size()) * tri2.size();
    Chord7 = tri4.size() * tri2.size();
    Chord8 = tri3_1.size() * tri3_2.size();
    Clique10 = tri4.size() * (tri4.size() - 1) / 2;

        
    //     Path1 = star_u2.size() * star_v2.size();
    //     Path2 = (star_u2.size() * star_v1.size()) + (star_u1.size() * star_v2.size());
    //     
    //     Star1 = (star_u2.size() - 1) * (star_u2.size() -1) + (star_v2.size() - 1) * (star_v2.size() -1);
    // }
}

void getTri2(VertexSet& vset, Edge& e, vector<vector<int>>& adj1, vector<vector<int>>& adj2,
            unordered_set<int>& star1_u, unordered_set<int>& star1_v, unordered_set<int>& star2_u, unordered_set<int>& star2_v, 
            unordered_set<int>& tri1, unordered_set<int>& tri2_1, unordered_set<int>& tri2_2, unordered_set<int>& tri3, 
            vector<int>& X1, vector<int>& X2,
            long long& Tri1, long long& Tri3,
            long long& Star1,
            long long& Path1, long long& Path2, long long& Path4,
            long long& Tailed1, long long& Tailed2, long long& Tailed4, long long& Tailed5, long long& Tailed7,
            long long& Chord1, long long& Chord3, long long& Chord4, long long& Chord6,
            long long& Clique8) {

    // edge : extended edge

    int u = e.src;
    int v = e.dest;

    for (auto w : adj1[u]) {
        if (w == v) continue;
        X1[w] = 1;
        star2_u.insert(w);
    }

    for (auto w : adj2[u]) {
        if (w == v) continue;
        X2[w] = 1;
        star1_u.insert(w);
    }

    for (auto w : adj1[v]) {
        if (w == u) continue;
        if (X2[w] == 1){
            star1_u.erase(w);
            tri2_1.insert(w);
            X1[w] = 3;
            continue;
        }
        if (X1[w] == 1) {
            Tri3++;
            X1[w] = 4;
            tri3.insert(w);
            star2_u.erase(w);
        }
        else {
            star2_v.insert(w);
            X1[w] = 2;     
        }
    }
    for (auto w : adj2[v]) {
        if (w == u) continue;
        if (X1[w] == 1){
            tri2_2.insert(w);
            star2_u.erase(w);
            X2[w] = 3;
            continue;
        }
        if (X2[w] == 1) {
            X2[w] = 4;
            Tri1++;
            tri1.insert(w);
            star1_u.erase(w);
        }
        else {
            star1_v.insert(w);
            X2[w] = 2;
        }
    }

    Star1 = (star1_u.size()  * (star1_u.size() -1))/2 + (star1_v.size()  * (star1_v.size() -1))/2;
    Path1 = star1_u.size() * star1_v.size();
    Path2 = (star2_u.size() * star1_v.size()) + (star1_u.size() * star2_v.size());
    Path4 = star2_u.size() * star2_v.size();
    Tailed1 = (star1_u.size() + star1_v.size()) * tri1.size();
    Tailed2 = (star2_u.size() + star2_v.size()) * tri1.size();
    Tailed4 = (star1_u.size() * tri2_1.size()) + (star1_v.size() * tri2_2.size());
    Tailed5 = (star2_u.size() * tri2_1.size()) + (star2_v.size() * tri2_2.size());
    Tailed7 = (star1_u.size() + star1_v.size()) * tri3.size();
    Chord1 = tri1.size() * (tri1.size() -1) / 2;
    Chord3 = tri1.size() * (tri2_1.size() + tri2_2.size());
    Chord4 = tri2_1.size() * tri2_2.size();
    Chord6 = tri3.size() * tri1.size();
    Clique8 = tri3.size() * (tri3.size() -1) / 2;

           
    //     Path2 = (star_u2.size() * star_v1.size()) + (star_u1.size() * star_v2.size());
    //     Path3 = star_u1.size() * star_v1.size();
    //     Star1 = (star_u2.size() - 1) * (star_u2.size() -1) + (star_v2.size() - 1) * (star_v2.size() -1);
    // 
}



void getCycle(vector<int> X, unordered_set<int>& star, vector<vector<int>>& adj, long long& cycle){
    for (int w : star){
        for (int r : adj[w]){
            if (X[r] == 2) {
                cycle +=1;
            }
        }
        X[w] = 0;
    }
}


void getClique(vector<int> X, int value, unordered_set<int>& tri, vector<vector<int>>& adj, long long& clique){
    for (auto w : tri){
        for (auto r : adj[w]){
            if (X[r] == value)
                clique +=1;
        }
        X[w] = 0;
    }
}



map<string, long long> countMotifs(VertexSet& vset, EdgeSet& e1, EdgeSet& e2, vector<vector<int>>& adj1, vector<vector<int>>& adj2) {
    map<string, long long> motifCounts = {
        {"Tri1", 0}, {"Tri2", 0}, {"Tri3", 0}, {"Tri4", 0},
        //clique
        {"4-1-1", 0}, {"4-1-2", 0}, {"4-1-3", 0}, {"4-1-4", 0}, {"4-1-5", 0}, {"4-1-6", 0}, {"4-1-7", 0}, {"4-1-8", 0}, {"4-1-9", 0}, {"4-1-10", 0}, {"4-1-11", 0},
        //chordalcycle
        {"4-2-1", 0}, {"4-2-2", 0}, {"4-2-3", 0}, {"4-2-4", 0}, {"4-2-5", 0}, {"4-2-6", 0}, {"4-2-7", 0}, {"4-2-8", 0},
        //tailedtriangle
        {"4-3-1", 0}, {"4-3-2", 0}, {"4-3-3", 0}, {"4-3-4", 0}, {"4-3-5", 0}, {"4-3-6", 0},
        //cycle
        {"4-4-1", 0}, {"4-4-2", 0}, {"4-4-3", 0},
        //star
        {"4-5-1", 0}, {"4-5-2", 0},
        //path
        {"4-6-1", 0}, {"4-6-2", 0}, {"4-6-3", 0}, {"4-6-4", 0}
    };
    int n = adj1.size();

    int max_num_workers = omp_get_max_threads() - 1;
    omp_set_num_threads(max_num_workers);
    int total_e1 = e1.size();
    int total_e2 = e2.size();
    printf("Get Max Thread : %d\n", max_num_workers);
    vector<map<string, long long>> thread_motifCounts(max_num_workers, motifCounts);

    int progress1 = 0;
    int next_progress1_update = total_e1 / 10;

    int progress2 = 0;
    //int next_progress2_update = total_e2 / 20;
    // start parallel
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
    
        #pragma omp for schedule(dynamic) 
        for (Edge edge : e1.getEdges()) {
            unordered_set<int> star2_u; unordered_set<int> star2_v;
            unordered_set<int> tri2_set; unordered_set<int> tri3_1_set; unordered_set<int> tri3_2_set; unordered_set<int> tri4_set; 
            vector<int> X1(n); vector<int> X2(n);
            long long Tri2 = 0; long long Tri4 = 0;
            long long Path3 = 0;
            long long Star2 = 0;
            long long Tailed3 = 0; long long Tailed6 = 0; long long Tailed8 = 0;
            long long Chord2 = 0; long long Chord5 = 0; long long Chord7 = 0; long long Chord8 = 0;
            long long Clique10 = 0;
            getTri1(vset, edge, adj1, adj2, star2_u, star2_v, tri2_set, tri3_1_set, tri3_2_set, tri4_set, X1, X2, Tri2, Tri4, Path3, Star2, Tailed3, Tailed6, Tailed8, Chord2, Chord5, Chord7, Chord8, Clique10);
            long long Cycle3 = 0;
            
            long long Clique3 = 0; long long Clique7 = 0; long long Clique9 = 0; long long Clique11 = 0;
            getCycle(X2, star2_u, adj1, Cycle3);
            getClique(X2, 4, tri2_set, adj1, Clique3);
            getClique(X2, 4, tri4_set, adj2, Clique7);
            getClique(X2, 4, tri4_set, adj1, Clique9);
            getClique(X1, 4, tri4_set, adj1, Clique11);
            thread_motifCounts[thread_id]["Tri2"] += Tri2;
            thread_motifCounts[thread_id]["Tri4"] += Tri4;
            thread_motifCounts[thread_id]["4-1-3"] += Clique3;
            thread_motifCounts[thread_id]["4-1-7"] += Clique7;
            thread_motifCounts[thread_id]["4-1-9"] += Clique9;
            thread_motifCounts[thread_id]["4-1-10"] += Clique10;
            thread_motifCounts[thread_id]["4-1-11"] += Clique11;
            thread_motifCounts[thread_id]["4-2-2"] += Chord2;
            thread_motifCounts[thread_id]["4-2-5"] += Chord5;
            thread_motifCounts[thread_id]["4-2-7"] += Chord7;
            thread_motifCounts[thread_id]["4-2-8"] += Chord8;
            thread_motifCounts[thread_id]["4-3-3"] += Tailed3;
            thread_motifCounts[thread_id]["4-3-6"] += Tailed6;
            thread_motifCounts[thread_id]["4-3-8"] += Tailed8;
            thread_motifCounts[thread_id]["4-4-3"] += Cycle3;
            thread_motifCounts[thread_id]["4-5-2"] += Star2;
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
            unordered_set<int> star1_u; unordered_set<int> star1_v; unordered_set<int> star2_u; unordered_set<int> star2_v;
            unordered_set<int> tri1_set; unordered_set<int> tri2_1_set; unordered_set<int> tri2_2_set; unordered_set<int> tri3_set;
            vector<int> X1(n); vector<int> X2(n);
            long long Tri1 = 0; long long Tri3 = 0; 
            long long Star1 = 0;
            long long Path1 = 0; long long Path2 = 0; long long Path4 = 0;
            long long Tailed1 = 0; long long Tailed2 = 0; long long Tailed4 = 0; long long Tailed5 = 0; long long Tailed7 = 0;
            long long Chord1 = 0; long long Chord3 = 0; long long Chord4 = 0; long long Chord6 = 0;
            long long Clique8 = 0;
            getTri2(vset, edge, adj1, adj2, star1_u, star1_v, star2_u, star2_v, tri1_set, tri2_1_set, tri2_2_set, tri3_set, X1, X2, Tri1, Tri3, Star1, Path1, Path2, Path4, Tailed1, Tailed2, Tailed4, Tailed5, Tailed7, Chord1, Chord3, Chord4, Chord6, Clique8);
            long long Cycle1 = 0; long long Cycle2 = 0;
            long long Clique1 = 0; long long Clique2 = 0; long long Clique4 = 0; long long Clique5 = 0; long long Clique6 = 0; long long Clique10 = 0;
            getCycle(X2, star1_u, adj2, Cycle1);
            getCycle(X2, star1_u, adj1, Cycle2);
            getClique(X2, 4, tri1_set, adj2, Clique1);
            getClique(X2, 4, tri1_set, adj1, Clique2);
            getClique(X2, 4, tri3_set, adj2, Clique4);
            getClique(X1, 4, tri2_1_set, adj2, Clique5);
            getClique(X1, 4, tri2_2_set, adj2, Clique5);
            getClique(X2, 4, tri3_set, adj1, Clique6);
            thread_motifCounts[thread_id]["Tri1"] += Tri1;
            thread_motifCounts[thread_id]["Tri3"] += Tri3;
            thread_motifCounts[thread_id]["4-1-1"] += Clique1;
            thread_motifCounts[thread_id]["4-1-2"] += Clique2;
            thread_motifCounts[thread_id]["4-1-4"] += Clique4;
            thread_motifCounts[thread_id]["4-1-5"] += Clique5;
            thread_motifCounts[thread_id]["4-1-6"] += Clique6;
            thread_motifCounts[thread_id]["4-1-8"] += Clique8;
            thread_motifCounts[thread_id]["4-2-1"] += Chord1;
            thread_motifCounts[thread_id]["4-2-3"] += Chord3;
            thread_motifCounts[thread_id]["4-2-4"] += Chord4;
            thread_motifCounts[thread_id]["4-2-6"] += Chord6;
            thread_motifCounts[thread_id]["4-3-1"] += Tailed1;
            thread_motifCounts[thread_id]["4-3-2"] += Tailed2;
            thread_motifCounts[thread_id]["4-3-4"] += Tailed4;
            thread_motifCounts[thread_id]["4-3-5"] += Tailed5;
            thread_motifCounts[thread_id]["4-3-7"] += Tailed7;
            thread_motifCounts[thread_id]["4-4-1"] += Cycle1;
            thread_motifCounts[thread_id]["4-4-2"] += Cycle2;
            thread_motifCounts[thread_id]["4-5-1"] += Star1;
            thread_motifCounts[thread_id]["4-6-1"] += Path1;
            thread_motifCounts[thread_id]["4-6-2"] += Path2;
            thread_motifCounts[thread_id]["4-6-4"] += Path4;
            


            #pragma omp atomic
            progress2++;
            if (progress2 % 10000 == 0) {
                #pragma omp critical
                {
                    printf("Progress: %d / %d\n", progress2, total_e2);
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
    string directory_path = argv[1];
    map<string, vector<long long>> all_motif_counts;
    int file_count = 0;

    for (const auto& entry : fs::directory_iterator(directory_path)) {
        printf("Random Graph %d \n", file_count+1);
        if (entry.is_regular_file() && entry.path().extension() == ".mtx") {

            tuple<VertexSet, EdgeSet, EdgeSet, vector<vector<int>>,  vector<vector<int>>> result = readAndSet(entry.path().string());    
            VertexSet vset = get<0>(result);
            EdgeSet eset1 = get<1>(result);
            EdgeSet eset2 = get<2>(result);
            vector<vector<int>> adj1 = get<3>(result);
            vector<vector<int>> adj2 = get<4>(result);

            map <string, long long> results = countMotifs(vset, eset1, eset2, adj1, adj2);
            
            //graphlet equation
            results["Tri1"] /= 3; results["Tri4"] /= 3;

            // for cycle
            results["4-4-1"] /= 4; results["4-4-3"] /= 2;

            // for path
            results["4-6-1"] = results["4-6-1"] - 4*results["4-4-1"] - results["4-4-2"];
            results["4-6-3"] = results["4-6-3"] - 2*results["4-4-3"] - results["4-4-2"];
            results["4-6-2"] = results["4-6-2"] - 2*results["4-4-2"];
            results["4-6-4"] = results["4-6-4"] - 2*results["4-4-3"];

            //for Clique
            results["4-1-1"] /= 6; results["4-1-3"] /= 2;  results["4-1-5"] /= 2; results["4-1-6"] /= 3; results["4-1-7"] /= 3; results["4-1-11"] /= 6;
            results["4-1-10"] = results["4-1-10"] - 6*results["4-1-11"];
            results["4-1-8"] = (results["4-1-8"] - results["4-1-10"])/2;
            

            //for Chordalcycle
            results["4-2-1"] = results["4-2-1"] - 6*results["4-1-1"] - results["4-1-2"];
            results["4-2-2"] = results["4-2-2"] - results["4-1-2"] - 2 * results["4-1-3"];
            results["4-2-3"] = results["4-2-3"] - 4 * results["4-1-2"] - 2 * results["4-1-4"];
            results["4-2-4"] = results["4-2-4"] - 4 * results["4-1-3"] - results["4-1-5"];
            results["4-2-5"] = results["4-2-5"] - 2 * results["4-1-4"] - 2 * results["4-1-5"];
            results["4-2-6"] = results["4-2-6"] - results["4-1-4"] - 3 * results["4-1-6"];
            results["4-2-7"] = results["4-2-7"] - 3 * results["4-1-7"] - results["4-1-9"];
            results["4-2-8"] = results["4-2-8"] - results["4-1-5"] - 4 * results["4-1-8"];


            //for TailedTriangle
            results["4-3-1"] = (results["4-3-1"] - 4 * results["4-2-1"] - results["4-2-3"])/2;
            results["4-3-2"] = (results["4-3-2"] - results["4-2-3"] - 2 * results["4-2-6"])/2;
            results["4-3-3"] = (results["4-3-3"] - results["4-2-3"] - 2 * results["4-2-4"]);
            results["4-3-4"] = (results["4-3-4"] - results["4-2-3"]) / 2;
            results["4-3-5"] = (results["4-3-5"] - 2 * results["4-2-4"]) / 2;
            results["4-3-6"] = (results["4-3-6"] - results["4-2-5"] - 2 * results["4-2-8"]) / 2;
            results["4-3-7"] = (results["4-3-7"] - results["4-2-5"]);
            results["4-3-8"] = (results["4-3-8"] - 2 * results["4-2-7"]) / 2;

            //for Star
            results["4-5-1"] = (results["4-5-1"] - results["4-3-1"] - results["4-3-4"])/3;
            results["4-5-2"] = (results["4-5-2"] - results["4-3-2"] - results["4-3-5"]);

            for (const auto& motif : results) {
                all_motif_counts[motif.first].push_back(motif.second);
            }
            file_count++;
        }
    }

    if (file_count == 0) {
        cerr << "No .mtx files found in the specified directory." << endl;
        return 1;
    }

    cout << fixed << setprecision(1);
    cout << "Motif Averages:" << endl;
    for (const auto& motif : all_motif_counts) {
        double sum = 0;
        for (long long count : motif.second) {
            sum += count;
        }
        double average = sum / file_count;
        cout << "\"" << motif.first << "\"" << " : " << average << "," << endl;
    }

    return 0;
}
