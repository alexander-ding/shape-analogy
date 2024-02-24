#pragma once

#include <vector>
#include <list>
#include "Eigen/StdVector"
#include "Eigen/Dense"

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix2f);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3f);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3i);

struct HalfEdge;
struct Vertex;
struct Edge;
struct Face;

typedef std::list<HalfEdge*>::iterator HalfEdgeIt;
typedef std::list<Vertex*>::iterator VertexIt;
typedef std::list<Edge*>::iterator EdgeIt;
typedef std::list<Face*>::iterator FaceIt;

struct HalfEdge
{
    HalfEdgeIt id;
    HalfEdge* twin;
    HalfEdge* next;
    Vertex* vertex;
    Edge* edge;
    Face* face;
};

struct Vertex
{
    VertexIt id;
    HalfEdge* halfedge;

    Eigen::Vector3f coords;
};

struct Edge
{
    EdgeIt id;
    HalfEdge* halfedge;
};

struct Face
{
    FaceIt id;
    HalfEdge* halfedge;
};

class Mesh
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    void initFromVectors(const std::vector<Eigen::Vector3f> &vertices,
                         const std::vector<Eigen::Vector3i> &faces);
    void loadFromFile(const std::string &filePath);
    void saveToFile(const std::string &filePath);
    void validate();
private:
    std::list<HalfEdge*> m_halfedges;
    std::list<Vertex*> m_vertices;
    std::list<Edge*> m_edges;
    std::list<Face*> m_faces;


    HalfEdge* addHalfEdge();
    Vertex* addVertex();
    Edge* addEdge();
    Face* addFace();
};

