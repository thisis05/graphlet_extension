#include <iostream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <sstream>
#include <vector>
#include <map>

// 명령어를 실행하고 결과를 문자열로 반환하는 함수
std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

// 출력 파싱 및 합산을 위한 구조체
struct MotifCounts {
    int total_4_clique = 0;
    int total_4_chordcycle = 0;
    int total_4_tailed_tris = 0;
    int total_4_cycle = 0;
    int total_3_star = 0;
    int total_4_path = 0;
};

// 결과 문자열을 파싱하여 MotifCounts 구조체에 합산하는 함수
void parseAndAddOutput(const std::string& output, MotifCounts& counts) {
    std::istringstream stream(output);
    std::string line;
    bool in_section = false;

    while (std::getline(stream, line)) {
        if (line.find("total_4_clique") != std::string::npos) {
            in_section = true;
        }

        if (in_section) {
            std::istringstream linestream(line);
            std::string key;
            char eq;
            int value;
            if (linestream >> key >> eq >> value) {
                if (key == "total_4_clique") counts.total_4_clique += value;
                else if (key == "total_4_chordcycle") counts.total_4_chordcycle += value;
                else if (key == "total_4_tailed_tris") counts.total_4_tailed_tris += value;
                else if (key == "total_4_cycle") counts.total_4_cycle += value;
                else if (key == "total_3_star") counts.total_3_star += value;
                else if (key == "total_4_path") counts.total_4_path += value;
            } else {
                std::cout << "Failed to parse line: " << line << std::endl; // 파싱 실패 시 출력
            }

            if (line.find("total_4_path") != std::string::npos) {
                in_section = false;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <base filename>" << std::endl;
        return 1;
    }

    MotifCounts totalCounts;
    std::string baseFilename = argv[1];
    
    // 명령어와 인자를 정의합니다.
    std::string commandBase = "./pgd/pgd -f ./datasets/" + baseFilename + "_random/" + baseFilename + "_r";

    // 명령어를 여러 번 실행합니다.
    for (int i = 1; i <= 10; ++i) { // r1부터 r10까지 실행
        std::string command = commandBase + std::to_string(i) + ".mtx";
        std::cout << "Running command: " << command << std::endl;
        std::string output = exec(command.c_str());
        std::cout << "Output: " << output << std::endl;

        // 출력 파싱 및 합산
        parseAndAddOutput(output, totalCounts);
    }

    // 결과 출력
    std::cout << "Total Motif Counts:" << std::endl;
    std::cout << "\"" << "4-1" << "\" : " << totalCounts.total_4_clique << ","<< std::endl;
    std::cout << "\"" << "4-2" << "\" : " << totalCounts.total_4_chordcycle << "," << std::endl;
    std::cout << "\"" << "4-3" << "\" : " << totalCounts.total_4_tailed_tris << "," << std::endl;
    std::cout << "\"" << "4-4" << "\" : " << totalCounts.total_4_cycle << ","<< std::endl;
    std::cout << "\"" << "4-5" << "\" : " << totalCounts.total_3_star << ","<< std::endl;
    std::cout << "\"" << "4-6" << "\" : " << totalCounts.total_4_path << ","<< std::endl;

    return 0;
}
