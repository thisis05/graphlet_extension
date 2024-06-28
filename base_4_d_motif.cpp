#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <set>
#include <string>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

const int threshold = 2;

tuple<vector<vector<int>>, vector<vector<int>>> readGraph(const string& filename) {
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
    vector<vector<int>> adjMatrix(n, vector<int>(n,0));
    vector<vector<int>> adjVector(n);
    int u, v;
    while (file >> u >> v) {
        u--, v--;
        adjMatrix[u][v] = 1;
        adjMatrix[v][u] = 1;  
        adjVector[u].push_back(v);
        adjVector[v].push_back(u);
    }
    file.close();


    printf("Calculate Distance ... \n");
    vector<vector<int>> dist = adjMatrix;

    for (int start_node = 0; start_node < n; start_node++) {
        vector<int> second = adjVector[start_node];
        for (const auto& second_node : second) {
            vector<int> final = adjVector[second_node];
            for (const auto& final_node : final) {
                if (start_node != final_node && dist[start_node][final_node] != 1){
                    dist[start_node][final_node] = 2;
                    adjMatrix[start_node][final_node] = 1;
                }
                // else{
                //     dist[start_node][final_node] = 0;
                // }
            }
        }
    }

    return make_tuple(dist, adjMatrix);
}


void ExtendSubgraph(set<int>& subgraph, set<int> extension, const vector<unordered_set<int>>& adjList, const vector<vector<int>>& adjMatrix, vector<set<int>>& result, const int& v) {

    if (subgraph.size() == 4) {
        result.push_back(subgraph);
        return;
    }
    while (!extension.empty()) {
        int w = *extension.begin();
        extension.erase(extension.begin());
        set<int> newExtension = extension;
        for (int u : adjList[w]) {
            if (u > v) {
                bool isNeighbor = false;
                for (int subgraphNode : subgraph) {
                    if (adjMatrix[subgraphNode][u] > 0) {
                        isNeighbor = true;
                        break;
                    }
                }
                if (!isNeighbor) {
                    newExtension.insert(u);
                }
            }
        }
        
        subgraph.insert(w);
        //extension에 삽입할 때 더 빨라질 수 있는 방법이 있을까? 일단 여기까지
        ExtendSubgraph(subgraph, newExtension, adjList, adjMatrix, result, v);
        subgraph.erase(w);
    }
    return;
}

map<string, int> countMotifs(vector<vector<int>> dist, const vector<vector<int>>& adjMatrix) {
    
    int n = adjMatrix.size();

    vector<unordered_set<int>> adjList(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (adjMatrix[i][j] > 0 && i != j) { 
                // 거리가 2 이하이고 자기 자신이 아닌 경우
                adjList[i].insert(j);
            }
        }
    }

    map <string, int> motifCounts = {
        {"4-1", 0}, //clique
        {"4-2", 0}, //chordalcycle
        {"4-3", 0}, //tailedtriangle
        {"4-4", 0}, //cycle
        {"4-5", 0}, //star
        {"4-6", 0},  //path
        //clique
        {"4-1-1", 0}, {"4-1-2", 0}, {"4-1-3", 0}, {"4-1-4", 0}, {"4-1-5", 0}, {"4-1-6", 0}, {"4-1-7", 0}, {"4-1-8", 0}, {"4-1-9", 0}, {"4-1-10", 0},
        //chordalcycle
        {"4-2-1", 0}, {"4-2-2", 0}, {"4-2-3", 0}, {"4-2-4", 0}, {"4-2-5", 0}, {"4-2-6", 0}, {"4-2-7", 0}, {"4-2-0", 0},
        //tailedtriangle
        {"4-3-1", 0}, {"4-3-2", 0}, {"4-3-3", 0}, {"4-3-4", 0}, {"4-3-5", 0}, {"4-3-6", 0}, {"4-3-0", 0}, 
        //cycle
        {"4-4-1", 0}, {"4-4-2", 0}, {"4-4-3", 0}, {"4-4-0", 0}, 
        //star
        {"4-5-1", 0}, {"4-5-2", 0}, {"4-5-0", 0},
        //path
        {"4-6-1", 0}, {"4-6-2", 0}, {"4-6-3", 0}, {"4-6-4", 0}, {"4-6-0", 0}
    };

    int max_num_workers = omp_get_max_threads();
    vector<map<string, int>> thread_motifCounts(max_num_workers, motifCounts);
    int progress = 1;
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        #pragma omp for schedule(dynamic)
        // 1st Loop: Iterate over each key in the map
        for (int a = 0; a < n; a++) {
            #pragma omp atomic
            progress++;
            if (progress % 50 == 0){
                printf("Progress : %d / %d \n", progress, n);
            }

            set<int> extension;
            for (int b : adjList[a]) {
                if (b > a) {
                    extension.insert(b);
                }
            }
            if (extension.size() == 0){
                continue;
            }

            set<int> subgraph = {a};
            vector<set<int>> result;
            //auto start = high_resolution_clock::now();
            ExtendSubgraph(subgraph, extension, adjList, adjMatrix, result, a);
            // auto end = high_resolution_clock::now();
            // auto duration = duration_cast<milliseconds>(end - start);
            // cout << "ExtendSubgraph time: " << duration.count() << " milliseconds" << endl;

            // auto start2 = high_resolution_clock::now();
            for (const auto& finalpath : result) {
                auto it = finalpath.begin();
                int a = *it++;
                int b = *it++;
                int c = *it++;
                int d = *it;

                //printf("%d %d %d %d \n", a, b, c, d);
                int degree = (adjMatrix[a][b] + adjMatrix[a][c] + adjMatrix[a][d] +
                                adjMatrix[b][c] + adjMatrix[b][d] + adjMatrix[c][d]);

                int total_distance = dist[a][b]+ dist[a][c] + dist[a][d] + dist[b][c] + dist[b][d] + dist[c][d];
                int dist_ab = dist[a][b];
                int dist_ac = dist[a][c];
                int dist_ad = dist[a][d];
                int dist_bc = dist[b][c];
                int dist_bd = dist[b][d];
                int dist_cd = dist[c][d]; 

                //printf("degree : %d and total distance : %d\n", degree, total_distance);

                if (degree == 6) {
                    thread_motifCounts[thread_id]["4-1"]++;
                    if (total_distance == 12) 
                        thread_motifCounts[thread_id]["4-1-1"]++;
                    else if (total_distance == 11) 
                        thread_motifCounts[thread_id]["4-1-2"]++;
                    else if (total_distance == 10) {
                        if ((dist_ab == 2 && dist_ac == 2 && dist_ad == 2) ||
                            (dist_ab == 2 && dist_bc == 2 && dist_bd == 2)||
                            (dist_ac == 2 && dist_bc == 2 && dist_cd == 2) ||
                            (dist_ad == 2 && dist_bd == 2 && dist_cd == 2)) {
                            thread_motifCounts[thread_id]["4-1-4"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-1-3"]++;
                        }   
                    }
                    else if (total_distance == 9) {
                        if ((dist_ab == 1 && dist_ac == 1 && dist_ad == 1) ||
                            (dist_ab == 1 && dist_bc == 1 && dist_bd == 1)||
                            (dist_ac == 1 && dist_bc == 1 && dist_cd == 1) ||
                            (dist_ad == 1 && dist_bd == 1 && dist_cd == 1)) {
                            //printf("%d %d %d %d\n", a, b, c, d);
                            thread_motifCounts[thread_id]["4-1-6"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-1-5"]++;
                        }   
                    }
                    else if (total_distance == 8) {
                        if ((dist_ab == 1 && dist_ac == 1 && dist_ad == 1) ||
                            (dist_ab == 1 && dist_bc == 1 && dist_bd == 1)||
                            (dist_ac == 1 && dist_bc == 1 && dist_cd == 1) ||
                            (dist_ad == 1 && dist_bd == 1 && dist_cd == 1)) {
                            thread_motifCounts[thread_id]["4-1-8"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-1-7"]++;
                        }   
                    }
                    else if (total_distance == 7) {
                        thread_motifCounts[thread_id]["4-1-9"]++;
                    }
                    else{
                        thread_motifCounts[thread_id]["4-1-10"]++;
                    }
                } else if (degree == 5) {
                    thread_motifCounts[thread_id]["4-2"]++;
                    if (total_distance == 10) {
                        thread_motifCounts[thread_id]["4-2-1"]++;
                    }
                    else if (total_distance == 9) {
                        if ((dist_ab == 2 && dist_ac == 2 && dist_ad == 2) ||
                            (dist_ab == 2 && dist_bc == 2 && dist_bd == 2)||
                            (dist_ac == 2 && dist_bc == 2 && dist_cd == 2) ||
                            (dist_ad == 2 && dist_bd == 2 && dist_cd == 2)) {
                            thread_motifCounts[thread_id]["4-2-3"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-2-2"]++;
                        } 
                    }
                    else if (total_distance == 8) {
                        
                        if ((dist_ab + dist_ac + dist_ad == 4) ||
                            (dist_ab + dist_bc + dist_bd == 4) ||
                            (dist_ac + dist_bc + dist_cd == 4) ||
                            (dist_ad + dist_bd + dist_cd == 4)) {
                            thread_motifCounts[thread_id]["4-2-5"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-2-4"]++;
                        } 
                    }
                    else if (total_distance == 7) {
                         if ((dist_ab + dist_ac + dist_ad == 3) ||
                            (dist_ab + dist_bc + dist_bd == 3) ||
                            (dist_ac + dist_bc + dist_cd == 3) ||
                            (dist_ad + dist_bd + dist_cd == 3)) {
                            thread_motifCounts[thread_id]["4-2-7"]++;
                        }
                        else{
                            
                            thread_motifCounts[thread_id]["4-2-6"]++;
                        } 
                    }
                    else {
                        thread_motifCounts[thread_id]["4-2-0"]++;
                    }
                } else if (degree == 4) {
                    if ((adjMatrix[a][b] && adjMatrix[a][c] && adjMatrix[a][d]) ||
                        (adjMatrix[b][a] && adjMatrix[b][c] && adjMatrix[b][d]) ||
                        (adjMatrix[c][a] && adjMatrix[c][b] && adjMatrix[c][d]) ||
                        (adjMatrix[d][a] && adjMatrix[d][b] && adjMatrix[d][c])) {
                        thread_motifCounts[thread_id]["4-3"]++;
                        if (total_distance == 8) {
                            thread_motifCounts[thread_id]["4-3-1"]++;
                        }
                        else if (total_distance == 7) {
                            if ((dist_ab + dist_ac + dist_ad == 3) ||
                                (dist_ab + dist_bc + dist_bd == 3) ||
                                (dist_ac + dist_bc + dist_cd == 3) ||
                                (dist_ad + dist_bd + dist_cd == 3)) {
                                thread_motifCounts[thread_id]["4-3-3"]++;
                            }
                            else{
                                thread_motifCounts[thread_id]["4-3-2"]++;
                            } 
                        }
                        else if (total_distance == 6) {
                            if  ((dist_ab + dist_ac + dist_ad == 5) ||
                                (dist_ab + dist_bc + dist_bd == 5) ||
                                (dist_ac + dist_bc + dist_cd == 5) ||
                                (dist_ad + dist_bd + dist_cd == 5)) {
                                thread_motifCounts[thread_id]["4-3-5"]++;
                            }
                            else{
                                thread_motifCounts[thread_id]["4-3-4"]++;
                            } 
                        }
                        else if (total_distance == 5) {
                            thread_motifCounts[thread_id]["4-3-6"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-3-0"]++;
                        }
                    } else {
                        thread_motifCounts[thread_id]["4-4"]++;
                        if (total_distance == 8) {
                            thread_motifCounts[thread_id]["4-4-1"]++;
                        }
                        else if (total_distance == 7) {
                            thread_motifCounts[thread_id]["4-4-2"]++;
                        }
                        else if (total_distance == 6) {
                            thread_motifCounts[thread_id]["4-4-3"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-4-0"]++;
                        }
                    }
                } else if (degree == 3) {
                    // printf("%d %d %d |", adjMatrix[a][b], adjMatrix[a][c], adjMatrix[a][d]);
                    // printf("%d %d %d |", adjMatrix[b][a], adjMatrix[b][c], adjMatrix[b][d]);
                    // printf("%d %d %d |", adjMatrix[c][a], adjMatrix[c][b], adjMatrix[c][d]);
                    // printf("%d %d %d \n", adjMatrix[d][a], adjMatrix[d][b], adjMatrix[d][c]);
                    if ((adjMatrix[a][b] && adjMatrix[a][c] && adjMatrix[a][d]) ||
                        (adjMatrix[b][a] && adjMatrix[b][c] && adjMatrix[b][d]) ||
                        (adjMatrix[c][a] && adjMatrix[c][b] && adjMatrix[c][d]) ||
                        (adjMatrix[d][a] && adjMatrix[d][b] && adjMatrix[d][c])) {
                        thread_motifCounts[thread_id]["4-5"]++;
                        if (total_distance == 6) {
                            thread_motifCounts[thread_id]["4-5-1"]++;
                        }
                        else if (total_distance == 5) {
                            thread_motifCounts[thread_id]["4-5-2"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-5-0"]++;
                        }
                    } else{
                        thread_motifCounts[thread_id]["4-6"]++;
                        if (total_distance == 6) {
                            thread_motifCounts[thread_id]["4-6-1"]++;
                        }
                        else if (total_distance == 5) {
                            if ((dist[a][b] + dist[a][c] + dist[a][d] == 1) ||
                                (dist[b][a] + dist[b][c] + dist[b][d] == 1) ||
                                (dist[c][a] + dist[c][b] + dist[c][d] == 1) ||
                                (dist[d][a] + dist[d][b] + dist[d][c] == 1)) {
                                thread_motifCounts[thread_id]["4-6-2"]++;
                            }
                            else{
                                thread_motifCounts[thread_id]["4-6-3"]++;
                            } 
                        }
                        else if (total_distance == 4) {
                            thread_motifCounts[thread_id]["4-6-4"]++;
                        }
                        else{
                            thread_motifCounts[thread_id]["4-6-0"]++;
                        }
                    }
                }
            }
            // auto end2 = high_resolution_clock::now();
            // auto duration2 = duration_cast<milliseconds>(end2 - start2);
            //cout << "MotifEnumeration time: " << duration2.count() << " milliseconds" << endl;
        }
    }
    for (int i = 0; i < max_num_workers; i++) {
        for (const auto& motif : thread_motifCounts[i]) {
            motifCounts[motif.first] += motif.second;
        }
    }
    return motifCounts;
}


void printAdjMatrix(vector<vector<int>> adjMatrix) {
    int n = adjMatrix.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
}

void printDistMatrix(vector<vector<int>> dist) {
    int n = dist.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << dist[i][j] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    string filename = argv[1];
    tuple<vector<vector<int>>, vector<vector<int>>> result = readGraph(filename);
    vector<vector<int>> dist = get<0>(result);
    vector<vector<int>> adjMatrix = get<1>(result);

    // cout << "Adjacency Matrix:" << endl;
    // printAdjMatrix(adjMatrix);
    // cout << "Dist Matrix:" << endl;
    // printDistMatrix(dist);
    map <string, int> results = countMotifs(dist, adjMatrix);

    cout << "Motif Counts:" << endl;
    for (const auto& motif : results) {
        cout << "\"" << motif.first << "\"" << " : " << motif.second << "," << endl;
    }

    return 0;
}