#include "graphlet_core.h"  // 라이브러리 헤더 파일 포함
#include <iostream>
#include <string>

using namespace graphlet;
using namespace std;

int main() {
    // 'graphlet_core' 클래스 인스턴스 생성 및 그래프 파일 읽기
    const string filename = "./sample_graph.csv";
    graphlet_core G(filename);

    // 그래프 모티프 분석
    G.graphlet_decomposition();

    std::cout << "Graphlet decomposition completed successfully." << std::endl;

    return 0;
}
