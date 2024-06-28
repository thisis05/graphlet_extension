#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>

using namespace std;

vector<vector<int>> readGraph(const string& filename) {
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

    vector<vector<int>> adjMatrix(n, vector<int>(n, 0));

    int u, v;
    while (file >> u >> v) {
        u--, v--;
        adjMatrix[u][v] = adjMatrix[v][u] = 1;
    }

    file.close();
    return adjMatrix;
}

unordered_map<string, int> countMotifs(const vector<vector<int>>& adjMatrix) {
    int n = adjMatrix.size();
    unordered_map<string, int> motifCounts = {
        {"4-1", 0}, //clique
        {"4-2", 0}, //chordalcycle
        {"4-3", 0}, //tailedtriangle
        {"4-4", 0}, //cycle
        {"4-5", 0}, //star
        {"4-6", 0}  //path
    };

    for (int a = 0; a < n; a++) {
        printf("Node %d / %d is done. \n", a, n);
        for (int b = a + 1; b < n; b++) {
            for (int c = b + 1; c < n; c++) {
                for (int d = c + 1; d < n; d++) {
                    int degree = adjMatrix[a][b] + adjMatrix[a][c] + adjMatrix[a][d] +
                                 adjMatrix[b][c] + adjMatrix[b][d] + adjMatrix[c][d];
                    if (degree == 6) {
                        motifCounts["4-1"]++;

                    } else if (degree == 5) {
                        motifCounts["4-2"]++;

                    } else if (degree == 4) {
                        if ((adjMatrix[a][b] && adjMatrix[a][c] && adjMatrix[a][d]) ||
                            (adjMatrix[b][a] && adjMatrix[b][c] && adjMatrix[b][d]) ||
                            (adjMatrix[c][a] && adjMatrix[c][b] && adjMatrix[c][d]) ||
                            (adjMatrix[d][a] && adjMatrix[d][b] && adjMatrix[d][c])) {
                            motifCounts["4-3"]++;

                        } else {
                            motifCounts["4-4"]++;

                        }
                    } else if (degree == 3) {
                        if ((adjMatrix[a][b] && adjMatrix[a][c] && adjMatrix[a][d]) ||
                            (adjMatrix[b][a] && adjMatrix[b][c] && adjMatrix[b][d]) ||
                            (adjMatrix[c][a] && adjMatrix[c][b] && adjMatrix[c][d]) ||
                            (adjMatrix[d][a] && adjMatrix[d][b] && adjMatrix[d][c])) {
                            motifCounts["4-5"]++;
                        } else{
                            int degree_a = adjMatrix[a][b] + adjMatrix[a][c] + adjMatrix[a][d];
                            int degree_b = adjMatrix[b][a] + adjMatrix[b][c] + adjMatrix[b][d];
                            int degree_c = adjMatrix[c][a] + adjMatrix[c][b] + adjMatrix[c][d];
                            int degree_d = adjMatrix[d][a] + adjMatrix[d][b] + adjMatrix[d][c];

                            if (degree_a != 0 && degree_b != 0 && degree_c != 0 && degree_d != 0) {
                                motifCounts["4-6"]++;
                            }
                        }
                    }
                }
            }
        }
    }

    return motifCounts;
}


void printAdjMatrix(const vector<vector<int>>& adjMatrix) {
    int n = adjMatrix.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << adjMatrix[i][j] << " ";
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
    vector<vector<int>> adjMatrix = readGraph(filename);

    // cout << "Adjacency Matrix:" << endl;
    // printAdjMatrix(adjMatrix);
    unordered_map<string, int> results = countMotifs(adjMatrix);

    cout << "Motif Counts:" << endl;
    for (const auto& motif : results) {
        cout << motif.first << ": " << motif.second << endl;
    }

    return 0;
}