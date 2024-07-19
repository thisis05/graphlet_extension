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

using namespace std;
using namespace std::chrono;

tuple<vector<pair<int, int>>, vector<pair<int, int>>, vector<vector<int>>, vector<vector<int>>> readAndSet(const string& filename) {
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
    
    while (file >> u >> v) {
        u--, v--;
        adjVector[v].push_back(u);
        adjVector[u].push_back(v);
    }
    file.close();
    
    vector<pair<int, int>> edgeSet1; 
    for (int u = 0; u < n; u++) {
        for (int v : adjVector[u]) {
            if (u < v)
                edgeSet1.emplace_back(u, v);
        }
    }

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


    vector<pair<int, int>> edgeSet2; 
    for (int u = 0; u < n; u++) {
        for (int v : adjVector2[u]) {
            if (u < v)
                edgeSet2.emplace_back(u, v);
        }
    }

    return make_tuple(edgeSet1, edgeSet2, adjVector, adjVector2);
}


void get3size(int& u, int& v, vector<vector<int>>& adj1, vector<vector<int>>& adj2,
            unordered_set<int>& star_u, unordered_set<int>& star_v, unordered_set<int>& star2_u, unordered_set<int>& star2_v, 
            vector<int>& X1, vector<int>& X2,
            long long& Star1, long long& Star2, 
            long long& Tri1, long long& Tri2, long long& Tri3, long long& Tri4) {

    // edge : original edge
    // star3_u : temp for original star version of motif

    for (auto w : adj1[u]) {
        if (w == v) continue;
        X1[w] = 1;
        Tri3++;
        star_u.insert(w);
    }

    for (auto w : adj1[v]) {
        if (w == u) continue;
        if (X1[w] == 1) {
            Tri3--;
            Tri4++;
            X1[w] = 3;
            star_u.erase(w);
        }
        else {
            Tri3++;
            star_v.insert(w);
            X1[w] = 2;     
        }
    }

    for (auto w : adj2[u]) {
        if(X1[w] == 2) continue;
        X2[w] = 1;
        Star2++;
        star2_u.insert(w);
    }

    for (auto w : adj2[v]) {
        if (X1[w] == 1) continue;
        if (X2[w] == 1) {
            X2[w] = 3;
            Tri2++;
            Star2--;
            star2_u.erase(w);
        }
        else {
            Star2++;
            star2_v.insert(w);
            X2[w] = 2;
        }
    }

    Star1 = star_u.size() * star2_v.size();
    Tri1 = (star_u.size() * (star_u.size() - 1) / 2) + (star_v.size() * (star_v.size() - 1) / 2); 
}

void getTri1(vector<int> X, unordered_set<int>& star, int value, vector<vector<int>>& adj, long long& Tri){
    for (int w : star){
        for (int r : adj[w]){
            if (X[r] == value) {
                Tri +=1;
            }
        }
        X[w] = 0;
    }
}



map<string, long long> countMotifs(vector<pair<int, int>>& e1, vector<pair<int, int>>& e2, vector<vector<int>>& adj1, vector<vector<int>>& adj2) {
    map<string, long long> motifCounts = {
        {"Star1", 0}, {"Star2", 0},
        {"Tri1", 0}, {"Tri2", 0}, {"Tri3", 0}, {"Tri4", 0}
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
        for (int i = 0; i < e1.size(); i++){
            int u = e1[i].first;
            int v = e1[i].second;
            unordered_set<int> star_u; unordered_set<int> star_v; unordered_set<int> star2_u; unordered_set<int> star2_v;
            vector<int> X1(n); vector<int> X2(n);
            long long Star1 = 0; long long Star2 = 0;
            long long Tri1 = 0; long long Tri2 = 0; long long Tri3 = 0; long long Tri4 = 0;
            get3size(u, v, adj1, adj2, star_u, star_v, star2_u, star2_v, X1, X2, Star1, Star2, Tri1, Tri2, Tri3, Tri4);
            long long Tri1_plus = 0; long long Tri1_minus = 0; 
            getTri1(X2, star2_u, 2, adj1, Tri1_plus);
            printf("%lld\n", Tri1_plus);
            getTri1(X1, star_u, 1, adj1, Tri1_minus);
            getTri1(X1, star_v, 1, adj1, Tri1_minus);
            //getTri1(X2, star2_v, 1, adj2, Tri1_2);
            Tri1 = Tri1 + (Tri1_plus / 6) - Tri1_minus;
            thread_motifCounts[thread_id]["Star1"] += Star1;
            thread_motifCounts[thread_id]["Star2"] += Star2;
            thread_motifCounts[thread_id]["Tri1"] += Tri1;
            thread_motifCounts[thread_id]["Tri2"] += Tri2;
            thread_motifCounts[thread_id]["Tri3"] += Tri3;
            thread_motifCounts[thread_id]["Tri4"] += Tri4;
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
    tuple<vector<pair<int, int>>, vector<pair<int, int>>, vector<vector<int>>,  vector<vector<int>>> result = readAndSet(filename);

    vector<pair<int, int>> eset1 = get<0>(result);
    vector<pair<int, int>> eset2 = get<1>(result);
    vector<vector<int>> adj1 = get<2>(result);
    vector<vector<int>> adj2 = get<3>(result);

    map <string, long long> results = countMotifs(eset1, eset2, adj1, adj2);

    //graphlet equation
    results["Tri1"] /= 3; results["Tri3"] /= 2; results["Tri4"] /= 3;

    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    double seconds = duration.count() / 1000.0; // Convert milliseconds to seconds

    cout << fixed << setprecision(3) << "Execution time: " << seconds << " seconds" << endl;

    std::cout << "Number of edges in EdgeSet1: " << eset1.size() << std::endl;
    std::cout << "Number of edges in EdgeSet2: " << eset2.size() << std::endl;
    cout << "Motif Counts:" << endl;
    for (const auto& motif : results) {
        cout << "\"" << motif.first << "\"" << " : " << motif.second << "," << endl;
    }

    return 0;
}
