#include "edge.h"

Edge::Edge() : src(0), dest(0), degree(0), Tri1(0), Tri2(0), Tri3(0), Tri4(0) {}

Edge::Edge(int s, int d, int sDegree, int dDegree) : src(s), dest(d), degree(sDegree + dDegree), Tri1(0), Tri2(0), Tri3(0), Tri4(0) {}

void Edge::setTri1(int tri) {
    Tri1 = tri;
}

int Edge::getTri1() const {
    return Tri1;
}

void Edge::setTri2(int tri) {
    Tri2 = tri;
}

int Edge::getTri2() const {
    return Tri2;
}

void Edge::setTri3(int tri) {
    Tri3 = tri;
}

int Edge::getTri3() const {
    return Tri3;
}

void Edge::setTri4(int tri) {
    Tri4 = tri;
}

int Edge::getTri4() const {
    return Tri4;
}

void Edge::print() const {
    std::cout << "Edge(" << src << ", " << dest << ", degree: " << degree << ")" << std::endl;
}

// EdgeSet 클래스 정의

void EdgeSet::addEdge(int src, int dest, int srcDegree, int destDegree) {
    edges.emplace_back(src, dest, srcDegree, destDegree);
}

void EdgeSet::print() const {
    for (const auto& edge : edges) {
        edge.print();
    }
}

size_t EdgeSet::size() const {
    return edges.size();
}

const std::vector<Edge>& EdgeSet::getEdges() const {
    return edges;
}

void EdgeSet::sortEdgesByDegree() {
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        return a.degree > b.degree;
    });
}
