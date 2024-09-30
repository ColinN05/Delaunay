#include "delaunay.h"
#include <cmath>

Triangle::Triangle(Vertex v1, Vertex v2, Vertex v3)
    : a(v1), b(v2), c(v3)
{
    edges[0] = {a,b};
    edges[1] = {a,c};
    edges[2] = {b,c};
}

bool operator==(const Triangle& t1, const Triangle t2)
{
    // Check if the vertices are equal in any order
    return 
       (t1.a == t2.a && t1.b == t2.b && t1.c == t2.c)
    || (t1.a == t2.a && t1.b == t2.c && t1.c == t2.b)
    || (t1.a == t2.b && t1.b == t2.a && t1.c == t2.c)
    || (t1.a == t2.b && t1.b == t2.c && t1.c == t2.a)
    || (t1.a == t2.c && t1.b == t2.a && t1.c == t2.b)
    || (t1.a == t2.c && t1.b == t2.b && t1.c == t2.a);
}

bool operator==(const Vertex& v1, const Vertex& v2)
{
    return v1.x == v2.x && v1.y == v2.y;
}

bool operator==(const Edge& e1, const Edge& e2)
{
    return (e1.a == e2.a && e1.b == e2.b)
    || (e1.a == e2.b && e1.b == e2.a);
}

std::vector<Triangle> bowyerWatsonTriangulation(const std::vector<Vertex>& vertices)
{
    std::unordered_set<Triangle, Triangle::HashFunction> triangulation; // store the triangles in a set rather than vector for fast insertion/deletion

    // Generate super triangle and add to triangulation
    double maxSquareDistance = 1.0, squareDistance;
    for (const Vertex& vertex : vertices)
    {
        squareDistance = vertex.x * vertex.x + vertex.y * vertex.y;
        if (squareDistance > maxSquareDistance)
        {
            maxSquareDistance = squareDistance;
        } 
    }
    double superRadius = 1.25 * sqrt(3) * sqrt(maxSquareDistance);
    Vertex superA, superB, superC;
    superA = {superRadius, 0.0};
    superB = {-0.5 * superRadius, 0.5 * sqrt(3) * superRadius};
    superC = {-0.5 * superRadius, -0.5 * sqrt(3) * superRadius};
    triangulation.insert(Triangle(superA, superB, superC));

    // Insertion of new vertex
    for (const Vertex& vertex : vertices)
    {
        // Triangles containing the new point in their circumcircles
        std::unordered_set<Triangle, Triangle::HashFunction> badTriangles;

        // Identify triangles that are now invalid due to insertion
        for (const Triangle& triangle : triangulation)
        {
            if (vertexInCircumcircle(vertex, triangle))
            {
                badTriangles.insert(triangle);
            }
        }

        // Edges of bad triangles that are not shared by two triangles
        std::unordered_set<Edge, Edge::HashFunction> polygon;
        for (const Triangle& badTriangle1 : badTriangles)
        {
            for (const Edge& edge1 : badTriangle1.edges)
            {
                bool outerEdge = true;
                for (const Triangle& badTriangle2 : badTriangles)
                {
                    if (badTriangle1 == badTriangle2) continue;
                    for (const Edge& edge2 : badTriangle2.edges)
                    {
                        if (edge1 == edge2) 
                        {
                            outerEdge = false;
                            goto END_EDGE_CHECK;
                        }
                    }
                }
                END_EDGE_CHECK:
                if (outerEdge)
                {
                    polygon.insert(edge1);
                }
            }
        }

        for (const Triangle& badTriangle : badTriangles)
        {
            triangulation.erase(badTriangle);
        }

        for (const Edge& edge : polygon)
        {
            Triangle newTri(vertex, edge.a, edge.b);
            triangulation.insert(newTri);
        }
    }

    // Erase triangles containing a vertex of the super triangle
    std::vector<decltype(triangulation)::key_type> outerTriangleKeys;
    for (const Triangle& t : triangulation)
    {
        if 
        (
            t.a == superA || t.b == superA || t.c == superA ||
            t.a == superB || t.b == superB || t.c == superB ||
            t.a == superC || t.b == superC || t.c == superC
        )
        {
            outerTriangleKeys.emplace_back(t);
        }
    }
    for (const decltype(triangulation)::key_type& key : outerTriangleKeys)
    {
        triangulation.erase(key);
    }

    std::vector<Triangle> triangulationVector = {};

    triangulationVector.insert(triangulationVector.end(), triangulation.begin(), triangulation.end());
    return triangulationVector;
}


bool vertexInCircumcircle(const Vertex& vertex, const Triangle& triangle)
{
    // orientation > 0 iff vertices of triangle are ordered counterclockwise, < 0 iff ordered clockwise
    double orientation = (triangle.b.x-triangle.a.x) * (triangle.c.y-triangle.a.y) - (triangle.c.x-triangle.a.x) * (triangle.b.y-triangle.a.y);

    Vertex a, b, c;
    if (orientation > 0)
    {
        a = triangle.a;
        b = triangle.b;
        c = triangle.c;
    }
    else
    {
        a = triangle.b;
        b = triangle.a;
        c = triangle.c;
    }
    const Vertex& d = vertex;
    double matrix[3][3] = 
    {
        {a.x-d.x, a.y-d.y, (a.x-d.x) * (a.x-d.x) + (a.y - d.y) * (a.y - d.y)},
        {b.x-d.x, b.y-d.y, (b.x-d.x) * (b.x-d.x) + (b.y - d.y) * (b.y - d.y)},
        {c.x-d.x, c.y-d.y, (c.x-d.x) * (c.x-d.x) + (c.y - d.y) * (c.y - d.y)}
    };
    return determinant3x3(matrix) > 0.0;
}

double determinant3x3(double A[3][3])
{
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
    - A[1][0] * (A[0][1] * A[2][2] - A[0][2] * A[2][1])
    + A[2][0] * (A[0][1] * A[1][2] - A[0][2] * A[1][1]);
}