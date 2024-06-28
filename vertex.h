#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>
#include <unordered_map>

// Vertex 클래스 정의
class Vertex {
public:
    int id;
    int degree1;
    int degree2;
    int Tri1; 
    int Tri2;
    int Tri3;
    int Tri4;

    Vertex();
    Vertex(int id);

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

// VertexSet 클래스 정의
class VertexSet {
private:
    std::unordered_map<int, Vertex> vertices;

public:
    VertexSet(int numVertices); // 전체 노드 개수를 받는 생성자
    void addVertex(int id);
    void setDegree1(int id, int degree);
    void setDegree2(int id, int degree);
    int getDegree1(int id) const;
    int getDegree2(int id) const;
    void print() const;

    void setTri1(int id, int tri);
    void setTri2(int id, int tri);
    void setTri3(int id, int tri);
    void setTri4(int id, int tri);
    int getTri1(int id) const;
    int getTri2(int id) const;
    int getTri3(int id) const;
    int getTri4(int id) const;
};

#endif // VERTEX_H
