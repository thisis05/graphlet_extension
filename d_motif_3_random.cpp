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

namespace fs = std::filesystem;
using namespace std;
using namespace std::chrono;
int max_num_workers = omp_get_max_threads()-1;

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
        adjVector[u].push_back(v);
        adjVector[v].push_back(u);
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


void getTri1(int& u, int& v, vector<vector<int>>& adj1, vector<vector<int>>& adj2,
            unordered_set<int>& star2_u, unordered_set<int>& star2_v, 
            unordered_set<int>& tri2, unordered_set<int>& tri3_1, unordered_set<int>& tri3_2, unordered_set<int>& tri4, 
            vector<int>& X1, vector<int>& X2,
            long long& Star2_small,
            long long& Tri2,  long long& Tri4) {

    // edge : original edge
    // star3_u : temp for original star version of motif

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

    Star2_small = star2_u.size() + star2_v.size();

}

void getTri2(int& u, int& v, vector<vector<int>>& adj1, vector<vector<int>>& adj2,
            unordered_set<int>& star1_u, unordered_set<int>& star1_v, unordered_set<int>& star2_u, unordered_set<int>& star2_v, 
            unordered_set<int>& tri1, unordered_set<int>& tri2_1, unordered_set<int>& tri2_2, unordered_set<int>& tri3, 
            vector<int>& X1, vector<int>& X2,
            long long& Star1_small,
            long long& Tri1, long long& Tri3) {


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

    Star1_small = star1_u.size() + star1_v.size();
}


map<string, long long> countMotifs(vector<pair<int, int>>& e1, vector<pair<int, int>>& e2, vector<vector<int>>& adj1, vector<vector<int>>& adj2) {
    map<string, long long> motifCounts = {
        {"Star1", 0}, {"Star2", 0},
        {"Tri1", 0}, {"Tri2", 0}, {"Tri3", 0}, {"Tri4", 0},
    };
    int n = adj1.size();

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
            unordered_set<int> star2_u; unordered_set<int> star2_v;
            unordered_set<int> tri2_set; unordered_set<int> tri3_1_set; unordered_set<int> tri3_2_set; unordered_set<int> tri4_set; 
            vector<int> X1(n); vector<int> X2(n);
            long long Star2_small = 0;
            long long Tri2 = 0; long long Tri4 = 0;
            getTri1(u, v, adj1, adj2, star2_u, star2_v, tri2_set, tri3_1_set, tri3_2_set, tri4_set, X1, X2, Star2_small, Tri2, Tri4);
            thread_motifCounts[thread_id]["Star2"] += Star2_small;
            thread_motifCounts[thread_id]["Tri2"] += Tri2;
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

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < e2.size(); i++){
            int u = e2[i].first;
            int v = e2[i].second;
            unordered_set<int> star1_u; unordered_set<int> star1_v; unordered_set<int> star2_u; unordered_set<int> star2_v;
            unordered_set<int> tri1_set; unordered_set<int> tri2_1_set; unordered_set<int> tri2_2_set; unordered_set<int> tri3_set;
            vector<int> X1(n); vector<int> X2(n);
            long long Star1_small = 0;
            long long Tri1 = 0; long long Tri3 = 0; 
            getTri2(u, v, adj1, adj2, star1_u, star1_v, star2_u, star2_v, tri1_set, tri2_1_set, tri2_2_set, tri3_set, X1, X2, Star1_small, Tri1, Tri3);
            thread_motifCounts[thread_id]["Star1"] += Star1_small;
            thread_motifCounts[thread_id]["Tri1"] += Tri1;
            thread_motifCounts[thread_id]["Tri3"] += Tri3;


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

            tuple<vector<pair<int, int>>, vector<pair<int, int>>, vector<vector<int>>,  vector<vector<int>>> result = readAndSet(entry.path().string());

            vector<pair<int, int>> eset1 = get<0>(result);
            vector<pair<int, int>> eset2 = get<1>(result);
            vector<vector<int>> adj1 = get<2>(result);
            vector<vector<int>> adj2 = get<3>(result);

            map <string, long long> results = countMotifs(eset1, eset2, adj1, adj2);

            //graphlet equation
            results["Star1"] /= 2;
            results["Tri1"] /= 3; results["Tri4"] /= 3;

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
