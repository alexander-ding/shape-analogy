#include "arap.h"
#include "mesh.h"

#include <limits>
#include <QtConcurrent>

using namespace std;
using namespace Eigen;

#define EPSILON 1e-6

// TODO: use qtconcurrent blockingmap
void ARAP::push(Vector3f position, Vector3f direction) {
    if (!this->m_isValid) this->build();

    MatrixXf vertices = this->m_mesh.getVertices();
    Vector3f planeNormal = direction.normalized();

    for (size_t i = 0; i< vertices.cols(); i++) {
        // checking if point is behind the plane
        Vector3f vertex = vertices.col(i);
        Vector3f planeToVertex = vertex - position;
        float dotProduct = planeToVertex.dot(planeNormal);
        if (dotProduct >= 0) {
            continue;
        }

        Vector3f pointOnPlane = vertex - dotProduct * planeNormal;

        vertices.col(i) = pointOnPlane;
    }

    this->m_mesh.setVertices(vertices);
}

Intersection ARAP::intersectTriangle(const Vector3f& start, const Vector3f& ray, const VectorXi& face) {
    Intersection result;
    result.hit = false;

    MatrixXf& vertices = this->m_mesh.getVertices();

    Vector3f v0 = vertices.col(face[0]);
    Vector3f v1 = vertices.col(face[1]);
    Vector3f v2 = vertices.col(face[2]);

    // Compute edges
    Vector3f edge1 = v1 - v0;
    Vector3f edge2 = v2 - v0;

    // Compute normal
    Vector3f normal = edge1.cross(edge2).normalized();

    if (normal.dot(ray) > 0) return result;

    // Compute intersection
    Vector3f h = ray.cross(edge2);
    float a = edge1.dot(h);

    if (a > -EPSILON && a < EPSILON) // Ray is parallel to triangle
        return result;

    float f = 1.0f / a;
    Vector3f s = start - v0;
    float u = f * s.dot(h);

    if (u < 0.0 || u > 1.0)
        return result;

    Vector3f q = s.cross(edge1);
    float v = f * ray.dot(q);

    if (v < 0.0 || u + v > 1.0)
        return result;

    float t = f * edge2.dot(q);
    if (t > EPSILON) {
        result.hit = true;
        result.position = start + t * ray;
        result.normal = normal;
        result.t = t;
    }

    return result;
}

Intersection ARAP::intersectMesh(Vector3f start, Vector3f ray) {
    bool hit = false;
    Intersection intersection;

    float minT = numeric_limits<float>::max();
    ray = ray.normalized();

    MatrixXi& faces = this->m_mesh.getFaces();
    for (size_t i=0; i < faces.cols(); i++) {
        Intersection testIntersect = intersectTriangle(start, ray, faces.col(i));
        if (testIntersect.hit && testIntersect.t < minT) {
            intersection = testIntersect;
        }
    }

    return intersection;
}
