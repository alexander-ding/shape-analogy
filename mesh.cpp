#include "mesh.h"

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <fstream>
#include <cassert>

#include <QFileInfo>
#include <QString>

#define TINYOBJLOADER_IMPLEMENTATION
#include "util/tiny_obj_loader.h"


using namespace Eigen;
using namespace std;

struct VertexPointerHash {
    template<typename T>
    std::size_t operator()(T it) const {
        // Use the address pointed to by the iterator as the basis for hashing
        return std::hash<const Vertex*>{}(it);
    }
};

struct PairVertexPointerHash {
    template<typename T>
    std::size_t operator()(const T& it_pair) const {
        auto hash1 = std::hash<const Vertex*>{}(it_pair.first);
        auto hash2 = std::hash<const Vertex*>{}(it_pair.second);
        return hash1 ^ (hash2 << 1);
    }
};

void Mesh::initFromVectors(const vector<Vector3f> &vertices,
                           const vector<Vector3i> &faces)
{
    vector<Vertex*> allocatedVertices;
    for (const auto &coords : vertices) {
        Vertex* vertex = this->addVertex();
        vertex->coords = coords;
        allocatedVertices.push_back(vertex);
    }
    unordered_map<pair<Vertex*, Vertex*>, HalfEdge*, PairVertexPointerHash> allocatedHalfEdges;
    for (const auto &indices : faces) {
        Face* face = this->addFace();
        Vertex* vertexA = allocatedVertices[indices[0]];
        Vertex* vertexB = allocatedVertices[indices[1]];
        Vertex* vertexC = allocatedVertices[indices[2]];
        HalfEdge* halfEdgeAB = this->addHalfEdge();
        HalfEdge* halfEdgeBC = this->addHalfEdge();
        HalfEdge* halfEdgeCA = this->addHalfEdge();
        face->halfedge = halfEdgeAB;
        vertexA->halfedge = halfEdgeAB;
        vertexB->halfedge = halfEdgeBC;
        vertexC->halfedge = halfEdgeCA;
        halfEdgeAB->vertex = vertexA;
        halfEdgeBC->vertex = vertexB;
        halfEdgeCA->vertex = vertexC;
        halfEdgeAB->face = face;
        halfEdgeBC->face = face;
        halfEdgeCA->face = face;
        halfEdgeAB->next = halfEdgeBC;
        halfEdgeBC->next = halfEdgeCA;
        halfEdgeCA->next = halfEdgeAB;
        auto halfEdgeBAIt = allocatedHalfEdges.find(make_pair(vertexB, vertexA));
        auto halfEdgeCBIt = allocatedHalfEdges.find(make_pair(vertexC, vertexB));
        auto halfEdgeACIt = allocatedHalfEdges.find(make_pair(vertexA, vertexC));

        if (halfEdgeBAIt != allocatedHalfEdges.end()) {
            auto halfEdgeBA = halfEdgeBAIt->second;
            halfEdgeAB->twin = halfEdgeBA;
            halfEdgeBA->twin = halfEdgeAB;
            Edge* edge = this->addEdge();
            edge->halfedge = halfEdgeAB;
            halfEdgeAB->edge = edge;
            halfEdgeBA->edge = edge;
        }
        if (halfEdgeCBIt != allocatedHalfEdges.end()) {
            auto halfEdgeCB = halfEdgeCBIt->second;
            halfEdgeBC->twin = halfEdgeCB;
            halfEdgeCB->twin = halfEdgeBC;
            Edge* edge = this->addEdge();
            edge->halfedge = halfEdgeBC;
            halfEdgeBC->edge = edge;
            halfEdgeCB->edge = edge;
        }
        if (halfEdgeACIt != allocatedHalfEdges.end()) {
            auto halfEdgeAC = halfEdgeACIt->second;
            halfEdgeCA->twin = halfEdgeAC;
            halfEdgeAC->twin = halfEdgeCA;
            Edge* edge = this->addEdge();
            edge->halfedge = halfEdgeCA;
            halfEdgeCA->edge = edge;
            halfEdgeAC->edge = edge;
        }
        allocatedHalfEdges.emplace(make_pair(vertexA, vertexB), halfEdgeAB);
        allocatedHalfEdges.emplace(make_pair(vertexB, vertexC), halfEdgeBC);
        allocatedHalfEdges.emplace(make_pair(vertexC, vertexA), halfEdgeCA);
    }
    this->validate();
}

void Mesh::loadFromFile(const string &filePath)
{
    tinyobj::attrib_t attrib;
    vector<tinyobj::shape_t> shapes;
    vector<tinyobj::material_t> materials;

    QFileInfo info(QString(filePath.c_str()));
    string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err,
                                info.absoluteFilePath().toStdString().c_str(), (info.absolutePath().toStdString() + "/").c_str(), true);
    if (!err.empty()) {
        cerr << err << endl;
    }

    if (!ret) {
        cerr << "Failed to load/parse .obj file" << endl;
        return;
    }

    vector<Vector3f> vertices;
    vector<Vector3i> faces;


    for (size_t s = 0; s < shapes.size(); s++) {
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            unsigned int fv = shapes[s].mesh.num_face_vertices[f];

            Vector3i face;
            for (size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

                face[v] = idx.vertex_index;
            }
            faces.push_back(face);

            index_offset += fv;
        }
    }
    for (size_t i = 0; i < attrib.vertices.size(); i += 3) {
        vertices.emplace_back(attrib.vertices[i], attrib.vertices[i + 1], attrib.vertices[i + 2]);
    }
    initFromVectors(vertices, faces);
    cout << "Loaded " << faces.size() << " faces and " << vertices.size() << " vertices" << endl;
}

void Mesh::saveToFile(const string &filePath)
{
    this->validate();
    ofstream outfile;
    outfile.open(filePath);

    unordered_map<Vertex*, size_t, VertexPointerHash> vertexToId;
    // Write vertices
    size_t i = 0;
    for (Vertex* v : this->m_vertices)
    {
        outfile << "v " << v->coords[0] << " " << v->coords[1] << " " << v->coords[2] << endl;
        vertexToId.emplace(v, i);
        i++;
    }

    // Write faces
    for (Face* f : this->m_faces)
    {
        Vertex* vA = f->halfedge->vertex;
        Vertex* vB = f->halfedge->next->vertex;
        Vertex* vC = f->halfedge->next->next->vertex;
        size_t indexA = vertexToId.find(vA)->second;
        size_t indexB = vertexToId.find(vB)->second;
        size_t indexC = vertexToId.find(vC)->second;

        outfile << "f " << (indexA+1) << " " << (indexB+1) << " " << (indexC+1) << endl;
    }

    outfile.close();
}

void Mesh::validate()
{
    unordered_set<HalfEdge*> validHalfEdges;
    for (HalfEdge* halfedge : this->m_halfedges) {

        validHalfEdges.emplace(halfedge);
    }

    unordered_set<Vertex*> validVertices;
    for (Vertex* vertex : this->m_vertices) {
        validVertices.emplace(vertex);
    }

    unordered_set<Face*> validFaces;
    for (Face* face : this->m_faces) {
        validFaces.emplace(face);
    }

    unordered_set<Edge*> validEdges;
    for (Edge* edge : this->m_edges) {
        validEdges.emplace(edge);
    }

    for (const HalfEdge* halfedge : this->m_halfedges) {
        // every halfedge should be a valid pointer
        assert(halfedge != nullptr);
        // every halfedge's ID backpointer should be correct
        assert(*(halfedge->id) == halfedge);
        // every halfedge should link back to itself
        assert(validHalfEdges.contains(halfedge->next));
        assert(validHalfEdges.contains(halfedge->next->next));
        assert(validHalfEdges.contains(halfedge->next->next->next));
        assert(halfedge->next->next->next == halfedge);
        // every halfedge should have a twin that points back to itself
        assert(validHalfEdges.contains(halfedge->twin));
        assert(halfedge->twin->twin == halfedge);
        // every halfedge's twin shouldn't the halfedge itself
        assert(halfedge->twin != halfedge);
        // every halfedge should have valid pointers for all fields
        assert(validEdges.contains(halfedge->edge));
        assert(validVertices.contains(halfedge->vertex));
        assert(validFaces.contains(halfedge->face));
    }

    for (const Vertex* vertex : this->m_vertices) {
        // every vertex should be a valid pointer
        assert(vertex != nullptr);
        // every vertex's ID backpointer should be correct
        assert(*(vertex->id) == vertex);
        // every vertex's halfedge should point back to it
        assert(validHalfEdges.contains(vertex->halfedge));
        assert(vertex->halfedge->vertex == vertex);
    }

    for (const Edge* edge : this->m_edges) {
        // every edge should be a valid pointer
        assert(edge != nullptr);
        // every edge's ID backpointer should be correct
        assert(*(edge->id) == edge);
        // every edge's halfedge and its twin halfedge should point back to the edge
        assert(validHalfEdges.contains(edge->halfedge));
        assert(edge->halfedge->edge == edge);
        assert(validHalfEdges.contains(edge->halfedge->twin));
        assert(edge->halfedge->twin->edge == edge);
    }

    for (const Face* face : this->m_faces) {
        // every face should be a valid pointer
        assert(face != nullptr);
        // every face's ID backpointer should be correct
        assert(*(face->id) == face);
        // every face's halfedge and its next two halfedges should point back to the face
        assert(validHalfEdges.contains(face->halfedge));
        assert(face->halfedge->face == face);
        assert(validHalfEdges.contains(face->halfedge->next));
        assert(face->halfedge->next->face == face);
        assert(validHalfEdges.contains(face->halfedge->next->next));
        assert(face->halfedge->next->next->face == face);
    }
}

HalfEdge* Mesh::addHalfEdge()
{
    HalfEdge* halfedge = new HalfEdge {};
    this->m_halfedges.push_back(halfedge);
    HalfEdgeIt it = this->m_halfedges.end();
    it--;
    halfedge->id = it;
    return halfedge;
}

Vertex* Mesh::addVertex()
{
    Vertex* vertex = new Vertex {};
    this->m_vertices.push_back(vertex);
    VertexIt it = this->m_vertices.end();
    it--;
    vertex->id = it;
    return vertex;
}

Edge* Mesh::addEdge()
{
    Edge* edge = new Edge {};
    this->m_edges.push_back(edge);
    EdgeIt it = this->m_edges.end();
    it--;
    edge->id = it;
    return edge;
}

Face* Mesh::addFace()
{
    Face* face = new Face {};
    this->m_faces.push_back(face);
    FaceIt it = this->m_faces.end();
    it--;
    face->id = it;
    return face;
}


