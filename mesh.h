#pragma once

#include <vector>
#include <unordered_set>
#include "Eigen/StdVector"
#include "Eigen/Dense"

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix2f);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3f);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3i);

struct HalfEdge;
struct Vertex;
struct Edge;
struct Face;

class Indexable {
public:
    virtual size_t getIndex() const = 0;
};


struct IndexablePointerHash {
    std::size_t operator()(Indexable* const &i) const {
        // Use the address pointed to by the iterator as the basis for hashing
        return std::hash<size_t>{}(i->getIndex());
    }
};

struct IndexablePointerEq {
    std::size_t operator()(Indexable* const &i1, Indexable* const &i2) const {
        // Use the address pointed to by the iterator as the basis for hashing
        return i1->getIndex() == i2->getIndex();
    }
};

template<typename T>
using IndexablePointerSet = std::unordered_set<T*, IndexablePointerHash, IndexablePointerEq>;

struct HalfEdge : public Indexable
{
    size_t index;
    HalfEdge* twin;
    HalfEdge* next;
    Vertex* vertex;
    Edge* edge;
    Face* face;

    size_t getIndex() const override {
        return this->index;
    }
};

struct Vertex : public Indexable
{
    size_t index;
    HalfEdge* halfedge;

    Eigen::Vector3f coords;

    size_t getIndex() const override {
        return this->index;
    }
};

struct Edge : public Indexable
{
    size_t index;
    HalfEdge* halfedge;

    size_t getIndex() const override {
        return this->index;
    }
};

struct Face : public Indexable
{
    size_t index;
    HalfEdge* halfedge;

    size_t getIndex() const override {
        return this->index;
    }
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

    size_t degree(const Vertex* vertex);
    IndexablePointerSet<HalfEdge> vertexHalfEdges(const Vertex* vertex);
    IndexablePointerSet<Vertex> vertexNeighbors(const Vertex* vertex);
    std::array<Vertex*, 4> edgeNeighbors(const Edge* edge);
    float edgeLength(const Edge* edge);
    std::array<Vertex*, 3> faceVertices(const Face* face);

    Edge* flipEdge(Edge* edge);
    Vertex* splitEdge(Edge* edge);
    Vertex* collapseEdge(Edge* edge);

    void loopSubdivide();
    void quadricErrorSimplify(size_t nEdges);
    void isotropicRemesh(float w = 0.5);

    IndexablePointerSet<Vertex>& getVertices();
    IndexablePointerSet<Edge>& getEdges();
    IndexablePointerSet<Face>& getFaces();

    Mesh();
private:
    IndexablePointerSet<HalfEdge> m_halfedges;
    IndexablePointerSet<Vertex> m_vertices;
    IndexablePointerSet<Edge> m_edges;
    IndexablePointerSet<Face> m_faces;
    size_t m_nextIndex;

    std::tuple<Vertex*, Edge*, Edge*> splitEdgeExtraInfo(Edge* edge);
    std::tuple<Vertex*, Vertex*, Edge*, Edge*, Edge*, Edge*, Edge*> collapseEdgeExtraInfo(Edge* edge);

    HalfEdge* addHalfEdge();
    Vertex* addVertex();
    Edge* addEdge();
    Face* addFace();

    void removeHalfEdge(HalfEdge* halfedge);
    void removeVertex(Vertex* vertex);
    void removeEdge(Edge* edge);
    void removeFace(Face* face);
};

