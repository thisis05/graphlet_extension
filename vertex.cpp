#include "vertex.h"

// Vertex 클래스 정의

Vertex::Vertex() : id(0), degree1(0), degree2(0), Tri1(0), Tri2(0), Tri3(0), Tri4(0){}

Vertex::Vertex(int id) : id(id), degree1(0), degree2(0), Tri1(0), Tri2(0), Tri3(0), Tri4(0) {}

void Vertex::setTri1(int tri) {
    Tri1 = tri;
}

int Vertex::getTri1() const {
    return Tri1;
}

void Vertex::setTri2(int tri) {
    Tri2 = tri;
}

int Vertex::getTri2() const {
    return Tri2;
}

void Vertex::setTri3(int tri) {
    Tri3 = tri;
}

int Vertex::getTri3() const {
    return Tri3;
}

void Vertex::setTri4(int tri) {
    Tri4 = tri;
}

int Vertex::getTri4() const {
    return Tri4;
}

void Vertex::print() const {
    std::cout << "Vertex(" << id << ", degree1: " << degree1 << ", degree2: " << degree2 << ", Tri1: " << Tri1 << ", Tri2: " << Tri2 << ")" << std::endl;
}

// VertexSet 클래스 정의

VertexSet::VertexSet(int numVertices) {
    vertices.reserve(numVertices); // 해시 테이블의 크기를 예약
}

void VertexSet::addVertex(int id) {
    vertices[id] = Vertex(id);
}

void VertexSet::setDegree1(int id, int degree) {
    vertices[id].degree1 = degree;
}

void VertexSet::setDegree2(int id, int degree) {
    vertices[id].degree2 = degree;
}

void VertexSet::setTri1(int id, int tri) {
    vertices[id].setTri1(tri);
}

void VertexSet::setTri2(int id, int tri) {
    vertices[id].setTri2(tri);
}

void VertexSet::setTri3(int id, int tri) {
    vertices[id].setTri3(tri);
}

void VertexSet::setTri4(int id, int tri) {
    vertices[id].setTri4(tri);
}

int VertexSet::getDegree1(int id) const {
    return vertices.at(id).degree1;
}

int VertexSet::getDegree2(int id) const {
    return vertices.at(id).degree2;
}

int VertexSet::getTri1(int id) const {
    return vertices.at(id).getTri1();
}

int VertexSet::getTri2(int id) const {
    return vertices.at(id).getTri2();
}

int VertexSet::getTri3(int id) const {
    return vertices.at(id).getTri3();
}

int VertexSet::getTri4(int id) const {
    return vertices.at(id).getTri4();
}

void VertexSet::print() const {
    for (const auto& pair : vertices) {
        pair.second.print();
    }
}
