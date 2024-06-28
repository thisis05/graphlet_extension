#ifndef EDGE_H
#define EDGE_H

#include <vector>
#include <algorithm>
#include <iostream>

// Edge 클래스
class Edge {
public:
    int src; // 시작 정점 (source vertex)
    int dest; // 도착 정점 (destination vertex)
    int degree; // 간선의 차수 (양 끝 정점의 차수 합)
    int Tri1; // Tri1 멤버 변수 추가
    int Tri2; // Tri2 멤버 변수 추가
    int Tri3;
    int Tri4;

    // 기본 생성자
    Edge();

    // 정점과 차수를 사용하는 생성자
    Edge(int s, int d, int sDegree, int dDegree);

    void setTri1(int tri);
    int getTri1() const;

    void setTri2(int tri);
    int getTri2() const;

    void setTri3(int tri);
    int getTri3() const;

    void setTri4(int tri);
    int getTri4() const;

    void print() const;
};

// EdgeSet 클래스
class EdgeSet {
private:
    std::vector<Edge> edges; // 간선들의 벡터

public:
    // 간선을 추가하는 함수
    void addEdge(int src, int dest, int srcDegree, int destDegree);

    // 모든 간선을 출력하는 함수
    void print() const;

    // 간선들의 개수를 반환하는 함수
    size_t size() const;

    // 간선들에 접근할 수 있는 함수
    const std::vector<Edge>& getEdges() const;

    // 간선을 차수에 따라 정렬하는 함수
    void sortEdgesByDegree();
};

#endif // EDGE_H
