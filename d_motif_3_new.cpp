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

vector<long long> vertices;
vector<int> edges;
map<int, vector<long long>> edgeSet;
unordered_map<int, unordered_map<int, short>> adjMat;
unsigned long long Star1, Star2;
unsigned long long Tri1, Tri2, Tri3, Tri4;
int n, m, max_degree;

void get_max_degree(){
    int degree = 0;
    max_degree = 0;
    for (long long v=0; v<n; v++) {
        degree = vertices[v+1] - vertices[v];
        if (max_degree < degree)  max_degree = degree;
    }
}

void readAndSet(const string& filename) {
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
    iss >> n >> n >> m;

    printf("Read File ... \n");
    map<int, vector<int>> adjVector;  // Use vector<vector<int>> for adjVector

    int u, v;
    
    while (file >> u >> v) {
        u--, v--;
        adjVector[v].push_back(u);
        adjVector[u].push_back(v);
        adjMat[u][v] = 1;
        adjMat[v][u] = 1;
        edgeSet[v].push_back(u);
    }
    file.close();

    printf("Preprocessing ... \n");
    vertices.push_back(0);
    for (int i=0; i < adjVector.size(); i++) {
        edges.insert(edges.end(), adjVector[i].begin(),adjVector[i].end());
        vertices.push_back(edges.size());
    }

    adjVector.clear();
    get_max_degree();
}


void get3size(long long & u, long long & v,  unsigned long long & tri_count, 
		vector<long long> & wedge, unsigned long long & wedge_count, unordered_map<int, unordered_map<int, short>>& adjMat) {
	
    wedge_count += vertices[u+1] - vertices[u] - 1;

    for (long long j = vertices[v]; j < vertices[v+1]; ++j) {
		long long w = edges[j];
		if (w==u) { continue; }
		if (adjMat[u][w] == 1) {
			tri_count++;
		}
		else {
            if (adjMat[u][w] == 0){
                wedge.push_back(w);
                adjMat[u][w] = 2;
            }
            wedge_count++;
		}
	}
}

void wedgebasedCount(long long& u, vector<long long> & wedge, unsigned long long & tri1_count, unsigned long long & tri2_count,
    unsigned long long& star1_count,  unordered_map<int, unordered_map<int, short>>& adjMat){
    
    for (long long w1 : wedge){
        for (long long w2 : wedge){
            if (w1 == w2) continue;
            short value = adjMat[w1][w2];
            if (value == 1){
                tri2_count++;
            }
            else if (value == 2){
                tri1_count++;
            }
            else {
                star1_count++;
            }
        }
    }
}


void countMotifs() {
    
    long long v, u, w;
    vector<long long> wedge(max_degree+1, 0);
    int max_num_workers = omp_get_max_threads() - 1;
    omp_set_num_threads(max_num_workers);
    printf("Get Max Thread : %d\n", max_num_workers);
    vector<vector<unsigned long long>> thread_motifCounts(max_num_workers, vector<unsigned long long>(6,0));
    int wedge_idx = 0;
    int progress1 = 0;
    int next_progress1_update = n / 10;
    
        
    // start parallel
    #pragma omp parallel for schedule(dynamic, 64) \
        firstprivate(adjMat) private(v,u,w,wedge)
    
    for (long long u = 0; u < n; u++){
        unsigned long long wedge_count = 0, tri_count = 0;
        int thread_id = omp_get_thread_num();
        for (long long j = vertices[u]; j < vertices[u+1]; ++j) {
            long long v = edges[j];
            get3size(v, u, tri_count, wedge, wedge_count, adjMat);
        }
        unsigned long long tri1_count = 0, tri2_count = 0, star1_count = 0;
        //wedgebasedCount(u, wedge, tri1_count, tri2_count, star1_count, adjMat);
        thread_motifCounts[thread_id][0] += star1_count;
        thread_motifCounts[thread_id][2] += tri1_count;
        thread_motifCounts[thread_id][3] += tri2_count;
        thread_motifCounts[thread_id][4] += wedge_count - tri_count;
        thread_motifCounts[thread_id][5] += tri_count;
        #pragma omp atomic
        progress1++;
        if (progress1 == next_progress1_update) {
            #pragma omp critical
            {
                if (progress1 >= next_progress1_update) {
                    printf("Edge Progress: %d / %d (%d%%)\n", progress1, n, (progress1 * 100) / n);
                    next_progress1_update += n / 10;  
                }
            }
        }
    } // end parallel

    vector<unsigned long long> motifCount(6, 0);
    for (int i = 0; i < max_num_workers; ++i) {
        motifCount[0] += thread_motifCounts[i][0];
        motifCount[1] += thread_motifCounts[i][1];
        motifCount[2] += thread_motifCounts[i][2];
        motifCount[3] += thread_motifCounts[i][3];
        motifCount[4] += thread_motifCounts[i][4];
        motifCount[5] += thread_motifCounts[i][5];
    }

    Star1 = motifCount[0];
    Star2 = motifCount[1];
    Tri1 = motifCount[2];
    Tri2 = motifCount[3];
    Tri3 = motifCount[4];
    Tri4 = motifCount[5];
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }
    auto start_time = high_resolution_clock::now();

    string filename = argv[1];
    readAndSet(filename);
    countMotifs();

    //graphlet equation
    Tri3 /= 4; Tri4 /= 6;

    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    double seconds = duration.count() / 1000.0; // Convert milliseconds to seconds

    std::cout << fixed << setprecision(3) << "Execution time: " << seconds << " seconds" << endl;

    std::cout << "Number of vertices : " << n << std::endl;
    std::cout << "Number of edges : " << m << std::endl;
    std::cout << "Max degree : " << max_degree << std::endl;

    std::cout << "Motif Counts:" << endl;
    std::cout << "\"Star1\"" << " : " << Star1 << "," << endl;
    std::cout << "\"Star2\"" << " : " << Star2 << "," << endl;
    std::cout << "\"Tri1\"" << " : " << Tri1 << "," << endl;
    std::cout << "\"Tri2\"" << " : " << Tri2 << "," << endl;
    std::cout << "\"Tri3\"" << " : " << Tri3 << "," << endl;
    std::cout << "\"Tri4\"" << " : " << Tri4 << endl;

    return 0;
}
