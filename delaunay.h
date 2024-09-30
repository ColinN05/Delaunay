#pragma once
#include <vector>
#include <unordered_set>

struct Vertex
{
    double x=0.0, y=0.0;
};

bool operator==(const Vertex& v1, const Vertex& v2);

struct Edge
{
    Vertex a,b;

    struct HashFunction
    {
        // need this to store Edge type in unordered_set
        size_t operator()(const Edge& e) const
        {
            size_t h1 = std::hash<int>()(e.a.x);
            size_t h2 = std::hash<int>()(e.a.y);
            size_t h3 = std::hash<int>()(e.b.x);
            size_t h4 = std::hash<int>()(e.b.y);
            return h1 ^ h2 ^ h3 ^ h4;
        }
    };
};

bool operator==(const Edge& e1, const Edge& e2);

struct Triangle
{
    Triangle(Vertex v1, Vertex v2, Vertex v3);
    Vertex a, b, c;
    Edge edges[3];

    struct HashFunction
    {
        // need this to store Triangle type in unordered_set
        size_t operator()(const Triangle& t) const
        {
            size_t h1 = std::hash<int>()(t.a.x);
            size_t h2 = std::hash<int>()(t.a.y);
            size_t h3 = std::hash<int>()(t.b.x);
            size_t h4 = std::hash<int>()(t.b.y);
            size_t h5 = std::hash<int>()(t.c.x);
            size_t h6 = std::hash<int>()(t.c.y);
            return h1 ^ h2 ^ h3 ^ h4 ^ h5 ^ h6;
        }
    };
};

bool operator==(const Triangle& t1, const Triangle t2);

std::vector<Triangle> bowyerWatsonTriangulation(const std::vector<Vertex>& vertices);
bool vertexInCircumcircle(const Vertex& vertex, const Triangle& triangle);
double determinant3x3(double A[3][3]);