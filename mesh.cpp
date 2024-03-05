#include "mesh.h"

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <utility>
#include <fstream>
#include <cassert>

#include <QFileInfo>
#include <QString>
#include "Eigen/Dense"

#define TINYOBJLOADER_IMPLEMENTATION
#include "util/tiny_obj_loader.h"


using namespace Eigen;
using namespace std;

bool isClose(float value, float target, float delta=0.0001) {
    return abs(value - target) < delta;
}

struct PairVertexPointerHash {
    template<typename T>
    std::size_t operator()(const T& it_pair) const {
        auto hash1 = IndexablePointerHash{}(it_pair.first);
        auto hash2 = IndexablePointerHash{}(it_pair.second);
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
        allocatedHalfEdges[make_pair(vertexA, vertexB)] = halfEdgeAB;
        allocatedHalfEdges[make_pair(vertexB, vertexC)] = halfEdgeBC;
        allocatedHalfEdges[make_pair(vertexC, vertexA)] = halfEdgeCA;
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

    unordered_map<Vertex*, size_t, IndexablePointerHash> vertexToId;
    // Write vertices
    size_t i = 0;
    for (Vertex* v : this->m_vertices)
    {
        outfile << "v " << v->coords[0] << " " << v->coords[1] << " " << v->coords[2] << endl;
        vertexToId[v] = i;
        i++;
    }

    // Write faces
    for (Face* f : this->m_faces)
    {
        Vertex* vA = f->halfedge->vertex;
        Vertex* vB = f->halfedge->next->vertex;
        Vertex* vC = f->halfedge->next->next->vertex;
        size_t indexA = vertexToId[vA];
        size_t indexB = vertexToId[vB];
        size_t indexC = vertexToId[vC];

        outfile << "f " << (indexA+1) << " " << (indexB+1) << " " << (indexC+1) << endl;
    }

    outfile.close();
}

void Mesh::validate()
{
    for (const HalfEdge* halfedge : this->m_halfedges) {
        // every halfedge should be a valid pointer
        assert(halfedge != nullptr);
        // every halfedge should link back to itself
        assert(this->m_halfedges.contains(halfedge->next));
        assert(this->m_halfedges.contains(halfedge->next->next));
        assert(this->m_halfedges.contains(halfedge->next->next->next));
        assert(halfedge->next->next->next == halfedge);
        // every halfedge should have a twin that points back to itself
        assert(this->m_halfedges.contains(halfedge->twin));
        assert(halfedge->twin->twin == halfedge);
        // every halfedge's twin shouldn't the halfedge itself
        assert(halfedge->twin != halfedge);
        // every halfedge should have valid pointers for all fields
        assert(this->m_edges.contains(halfedge->edge));
        assert(this->m_vertices.contains(halfedge->vertex));
        assert(this->m_faces.contains(halfedge->face));
    }

    for (const Vertex* vertex : this->m_vertices) {
        // every vertex should be a valid pointer
        assert(vertex != nullptr);
        // every vertex's halfedge should point back to it
        assert(this->m_halfedges.contains(vertex->halfedge));
        assert(vertex->halfedge->vertex == vertex);
    }

    for (const Edge* edge : this->m_edges) {
        // every edge should be a valid pointer
        assert(edge != nullptr);
        // every edge's halfedge and its twin halfedge should point back to the edge
        assert(this->m_halfedges.contains(edge->halfedge));
        assert(edge->halfedge->edge == edge);
        assert(this->m_halfedges.contains(edge->halfedge->twin));
        assert(edge->halfedge->twin->edge == edge);
    }

    for (const Face* face : this->m_faces) {
        // every face should be a valid pointer
        assert(face != nullptr);
        // every face's halfedge and its next two halfedges should point back to the face
        assert(this->m_halfedges.contains(face->halfedge));
        assert(face->halfedge->face == face);
        assert(this->m_halfedges.contains(face->halfedge->next));
        assert(face->halfedge->next->face == face);
        assert(this->m_halfedges.contains(face->halfedge->next->next));
        assert(face->halfedge->next->next->face == face);
    }
}

size_t Mesh::degree(const Vertex *vertex)
{
    HalfEdge* h = vertex->halfedge;
    size_t count = 0;
    do {
        count++;
        h = h->twin->next;
    } while (h != vertex->halfedge);
    return count;
}

IndexablePointerSet<HalfEdge> Mesh::vertexHalfEdges(const Vertex *vertex)
{
    IndexablePointerSet<HalfEdge> halfedges;
    HalfEdge* h = vertex->halfedge;
    do {
        halfedges.insert(h);
        h = h->twin->next;
    } while (h != vertex->halfedge);
    return halfedges;
}

IndexablePointerSet<Vertex> Mesh::vertexNeighbors(const Vertex *vertex)
{
    IndexablePointerSet<Vertex> neighbors;
    HalfEdge* h = vertex->halfedge;
    do {
        neighbors.insert(h->twin->vertex);
        h = h->twin->next;
    } while (h != vertex->halfedge);
    return neighbors;
}

std::array<Vertex *, 4> Mesh::edgeNeighbors(const Edge *edge)
{
    HalfEdge* h = edge->halfedge;
    return {h->vertex, h->twin->vertex, h->next->next->vertex, h->twin->next->twin->vertex};
}

inline float Mesh::edgeLength(const Edge *edge)
{
    HalfEdge* h = edge->halfedge;
    return (h->vertex->coords - h->twin->vertex->coords).norm();
}

std::array<Vertex *, 3> Mesh::faceVertices(const Face *face)
{
    HalfEdge* h = face->halfedge;
    return {h->vertex, h->next->vertex, h->next->next->vertex};
}

Edge* Mesh::flipEdge(Edge *edge)
{
    // validate edge flip: don't flip when at least one of the vertices has degree 3
    if (this->degree(edge->halfedge->vertex) == 3 || this->degree(edge->halfedge->twin->vertex) == 3) {
        return nullptr;
    }

    HalfEdge* bc = edge->halfedge;
    HalfEdge* cb = edge->halfedge->twin;
    HalfEdge* ca = bc->next;
    HalfEdge* ac = ca->twin;
    HalfEdge* ab = ca->next;
    HalfEdge* ba = ab->twin;
    HalfEdge* bd = cb->next;
    HalfEdge* db = bd->twin;
    HalfEdge* dc = bd->next;
    HalfEdge* cd = dc->twin;

    Edge* bcEdge = edge;
    Edge* caEdge = ca->edge;
    Edge* abEdge = ab->edge;
    Edge* dcEdge = dc->edge;
    Edge* bdEdge = bd->edge;

    Vertex* a = ac->vertex;
    Vertex* b = bc->vertex;
    Vertex* c = cd->vertex;
    Vertex* d = db->vertex;

    Face* abc = ab->face;
    Face* bdc = bd->face;

    // bc -> ad; cb -> da
    HalfEdge* ad = bc;
    HalfEdge* da = cb;
    // bcEdge -> adEdge
    Edge* adEdge = bcEdge;
    // abc -> abd
    Face* abd = abc;
    // bdc -> adc
    Face* adc = bdc;

    ad->next = dc;
    ad->vertex = a;
    ad->edge = adEdge;
    ad->face = adc;
    da->next = ab;
    da->vertex = d;
    da->edge = adEdge;
    da->face = abd;

    // ac->next doesn't change
    ac->vertex = a;
    ac->edge = caEdge;
    // ac->face doesn't change
    ca->next = ad;
    ca->vertex = c;
    ca->edge = caEdge;
    ca->face = adc;

    // cd->next doesn't change
    cd->vertex = c;
    cd->edge = dcEdge;
    // cd->face doesn't change
    dc->next = ca;
    dc->vertex = d;
    dc->edge = dcEdge;
    dc->face = adc;

    // ba->next doesn't change
    ba->vertex = b;
    ba->edge = abEdge;
    // ba->face doesn't change
    ab->next = bd;
    ab->vertex = a;
    ab->edge = abEdge;
    ab->face = abd;

    // db->next doesn't change
    db->vertex = d;
    db->edge = bdEdge;
    // db->face doesn't change
    bd->next = da;
    bd->vertex = b;
    bd->edge = bdEdge;
    bd->face = abd;

    b->halfedge = bd;
    c->halfedge = ca;

    adEdge->halfedge = ad;

    abc->halfedge = ab;
    bdc->halfedge = dc;

    return adEdge;
}

Vertex* Mesh::splitEdge(Edge *edge)
{
    auto [v, _, __] = this->splitEdgeExtraInfo(edge);
    return v;
}

Vertex *Mesh::collapseEdge(Edge *edge)
{
    auto [m, d, acEdge, cbEdge, daEdge, bdEdge, cdEdge] = this->collapseEdgeExtraInfo(edge);
    if (m == nullptr) {
        return nullptr;
    }
    this->removeVertex(d);
    this->removeEdge(daEdge);
    this->removeEdge(bdEdge);
    this->removeEdge(cdEdge);
    return m;
}

void Mesh::loopSubdivide()
{
    unordered_map<Vertex*, IndexablePointerSet<Vertex>, IndexablePointerHash> vertexNeighbors;
    for (Vertex* vertex : this->m_vertices) {
        vertexNeighbors[vertex] = this->vertexNeighbors(vertex);
    }

    std::vector<Edge*> edges(this->m_edges.begin(), this->m_edges.end());
    std::vector<std::array<Vertex*, 4>> edgeNeighbors;
    edgeNeighbors.reserve(edges.size());
    for (Edge* edge : edges) {
        edgeNeighbors.push_back(this->edgeNeighbors(edge));
    }

    std::vector<Edge*> newEdges;

    for (size_t i = 0; i < edges.size(); i++) {
        std::array<Vertex*, 4>& neighbors = edgeNeighbors[i];
        auto [vertex, newEdgeA, newEdgeB] = this->splitEdgeExtraInfo(edges[i]);
        vertex->coords = neighbors[0]->coords * (3.f / 8.f) + neighbors[1]->coords * (3.f / 8.f)
            + neighbors[2]->coords * (1.f / 8.f) + neighbors[3]->coords * (1.f / 8.f);
        newEdges.push_back(newEdgeA);
        newEdges.push_back(newEdgeB);
    }

    for (Edge* edge : newEdges) {
        auto aIt = vertexNeighbors.find(edge->halfedge->vertex);
        auto bIt = vertexNeighbors.find(edge->halfedge->twin->vertex);
        if ((aIt == vertexNeighbors.end() && bIt != vertexNeighbors.end()) ||
            (aIt != vertexNeighbors.end() && bIt == vertexNeighbors.end())) {
            this->flipEdge(edge);
        }
    }

    unordered_map<Vertex*, Eigen::Vector3f, IndexablePointerHash> vertexCoords;
    for (Vertex* vertex : this->m_vertices) {
        auto neighborsIt = vertexNeighbors.find(vertex);
        if (neighborsIt == vertexNeighbors.end()) {
            continue;
        }
        auto [_, neighbors] = *neighborsIt;
        size_t n = neighbors.size();
        float u = n == 3 ?
                      (3.f / 16.f) :
                      (1.f / n * (5.f/8.f - std::pow(3.f/8.f + 1.f/4.f * std::cos(2.f * M_PI / n), 2.f)));
        Eigen::Vector3f coords = (1.f - n * u) * vertex->coords;
        for (Vertex* neighbor : neighbors) {
            coords += neighbor->coords * u;
        }
        vertexCoords[vertex] = coords;
    }

    for (auto& [vertex, coords] : vertexCoords) {
        vertex->coords = coords;
    }
}

inline float evaluateError(const Eigen::Matrix4f& Q, const Eigen::Vector3f v) {
    Eigen::Vector4f vPrime = v.homogeneous();
    float error = vPrime.transpose() * Q * vPrime;
    return error;
}

std::pair<Eigen::Vector3f, float> findMinimalErrorPoint(Eigen::Matrix4f& Q, Vertex* v1, Vertex* v2) {
    Eigen::Matrix4f QOver = Q;
    QOver(3, 0) = 0;
    QOver(3, 1) = 0;
    QOver(3, 2) = 0;
    QOver(3, 3) = 1;
    float determinant = QOver.determinant();
    if (!isClose(determinant, 0)) {
        Eigen::Vector4f v = QOver.inverse() * Eigen::Vector4f {0, 0, 0, 1};
        Eigen::Vector3f coords = Eigen::Vector3f {v(0), v(1), v(2)};
        return std::make_pair(coords, evaluateError(Q, coords));
    }
    // TODO: find along segment
    std::cout << "cannot invert" << std::endl;
    // choose between v1, v2, or midpoint
    Eigen::Vector3f midpoint = (v1->coords + v2->coords) / 2.f;
    float v1Error = evaluateError(Q, v1->coords);
    float v2Error = evaluateError(Q, v2->coords);
    float midpointError = evaluateError(Q, midpoint);
    if (v1Error < v2Error) {
        return v1Error < midpointError ? std::make_pair(v1->coords, v1Error) : std::make_pair(midpoint, midpointError);
    } else {
        return v2Error < midpointError ? std::make_pair(v2->coords, v2Error) : std::make_pair(midpoint, midpointError);
    }
}

void Mesh::quadricErrorSimplify(size_t nEdges)
{
    // compute Q for each triangle
    unordered_map<Face*, Eigen::Matrix4f, IndexablePointerHash> faceQs;
    for (Face* face : this->m_faces) {
        std::array<Vertex*, 3> vertices = this->faceVertices(face);
        // use cross product to find normal
        Eigen::Vector3f normal = (vertices[1]->coords - vertices[0]->coords)
                                     .cross(vertices[2]->coords - vertices[0]->coords)
                                     .normalized();
        float d = (-vertices[0]->coords).dot(normal);
        Eigen::Vector4f v {normal[0], normal[1], normal[2], d};
        Eigen::Matrix4f Q = v * v.transpose();
        faceQs[face] = Q;
    }

    unordered_map<Vertex*, Eigen::Matrix4f, IndexablePointerHash> Qs;
    for (Vertex* vertex : this->m_vertices) {
        IndexablePointerSet<HalfEdge> neighbors = this->vertexHalfEdges(vertex);
        Eigen::Matrix4f totalQ = Eigen::Matrix4f::Zero();
        for (HalfEdge* edge : neighbors) {
            totalQ += faceQs[edge->face];
        }
        Qs[vertex] = totalQ;
    }

    using ErrorsMap = unordered_map<Edge*, tuple<Eigen::Vector3f, Eigen::Matrix4f, float>, IndexablePointerHash>;
    ErrorsMap edgeErrors;
    for (Edge* edge : this->m_edges) {
        Vertex* a = edge->halfedge->vertex;
        Vertex* b = edge->halfedge->twin->vertex;
        Eigen::Matrix4f Q = Qs[a] + Qs[b];
        auto [coords, error] = findMinimalErrorPoint(Q, a, b);
        edgeErrors[edge] = std::make_tuple(coords, Q, error);
    }

    auto cmp = [](pair<Edge*, ErrorsMap*> left, pair<Edge*, ErrorsMap*> right)
    {
        ErrorsMap* edgeErrors = left.second;
        auto leftError = get<2>((*edgeErrors)[left.first]);
        auto rightError = get<2>((*edgeErrors)[right.first]);

        return (leftError == rightError)
                   ? (left.first < right.first)
                   : (leftError < rightError);
    };
    // unordered_Map: Edge* -> error
    // ordered set: Edge*
    set<pair<Edge*, ErrorsMap*>, decltype(cmp)> candidateEdges;
    for (Edge* edge : this->m_edges) {
        candidateEdges.insert(std::make_pair(edge, &edgeErrors));
    }

    for (size_t i = 0; i < nEdges; i++) {
        bool foundEdge = false;

        for (auto [edge, _] : candidateEdges) {
            auto vertexA = edge->halfedge->vertex->coords;
            auto vertexB = edge->halfedge->twin->vertex->coords;
            auto [newVertex, d, acEdge, cbEdge, daEdge, bdEdge, cdEdge] = this->collapseEdgeExtraInfo(edge);
            if (newVertex == nullptr) {
                continue;
            }

            // remove went through
            foundEdge = true;

            // update affected vertices
            auto [coords, newVertexQ, edgeError] = edgeErrors[edge];
            newVertex->coords = coords;
            Qs[newVertex] = newVertexQ;
            Qs.erase(d);
            this->removeVertex(d);

            // remove all edges that need to be deleted
            vector<Edge*> edgesToRemove {daEdge, bdEdge, cdEdge};
            for (auto edge : edgesToRemove) {
                candidateEdges.erase(make_pair(edge, &edgeErrors));
            }
            for (auto edge : edgesToRemove) {
                edgeErrors.erase(edge);
                this->removeEdge(edge);
            }

            // recompute all edges that are affected
            IndexablePointerSet<HalfEdge> neighbors = this->vertexHalfEdges(newVertex);
            for (auto h : neighbors) {
                candidateEdges.erase(make_pair(h->edge, &edgeErrors));
            }

            for (auto h : neighbors) {
                Vertex* otherVertex = h->twin->vertex;
                Eigen::Matrix4f otherVertexQ = Qs[otherVertex];
                Eigen::Matrix4f Q = newVertexQ + otherVertexQ;
                auto [coords, error] = findMinimalErrorPoint(Q, newVertex, otherVertex);

                edgeErrors[h->edge] = make_tuple(coords, Q, error);
            }

            for (auto h : neighbors) {
                candidateEdges.insert(make_pair(h->edge, &edgeErrors));
            }

            assert(candidateEdges.size() == this->m_edges.size());
            assert(edgeErrors.size() == this->m_edges.size());
            assert(Qs.size() == this->m_vertices.size());

            // return;
            break;
        }

        if (!foundEdge) {
            // out of edges, unexpectedly; return early
            return;
        }
    }
}

void Mesh::isotropicRemesh(float w)
{
    // Compute the mean edge length of the input.
    float edgeLength = 0;
    for (Edge* edge : this->m_edges) {
        edgeLength += this->edgeLength(edge);
    }
    edgeLength /= this->m_edges.size();

    // Split all edges that are longer than 4/3 L
    float edgeSplitThreshold = edgeLength * 4.f / 3.f;
    std::vector<Edge> edgesToSplit;
    for (Edge* edge : this->m_edges) {
        if (this->edgeLength(edge) > edgeSplitThreshold) {
            edgesToSplit.push_back(*edge);
        }
    }
    for (Edge& edge : edgesToSplit) {
        auto edgeIt = this->m_edges.find(&edge);
        if (edgeIt != this->m_edges.end() && this->edgeLength(*edgeIt) > edgeSplitThreshold) {
            this->splitEdge(*edgeIt);
        }
    }

    // Collapse all edges that are shorter than 4/5 L
    float edgeCollapseThreshold = edgeLength * 4.f / 5.f;
    std::vector<Edge> edgesToCollapse;
    for (Edge* edge : this->m_edges) {
        if (this->edgeLength(edge) < edgeCollapseThreshold) {
            edgesToCollapse.push_back(*edge);
        }
    }
    for (Edge& edge : edgesToCollapse) {
        auto edgeIt = this->m_edges.find(&edge);
        if (edgeIt != this->m_edges.end() && this->edgeLength(*edgeIt) < edgeCollapseThreshold) {
            this->collapseEdge(*edgeIt);
        }
    }

    // Flip all edges that decrease the total deviation from degree 6.
    std::vector<Edge> edgesToFlip;
    for (Edge* edge : this->m_edges) {
        auto neighbors = this->edgeNeighbors(edge);
        size_t aDegree = this->degree(neighbors[0]);
        size_t bDegree = this->degree(neighbors[1]);
        size_t cDegree = this->degree(neighbors[2]);
        size_t dDegree = this->degree(neighbors[3]);
        size_t costNoChange = abs(static_cast<int>(aDegree - 6)) + abs(static_cast<int>(bDegree - 6))
                              + abs(static_cast<int>(cDegree - 6)) + abs(static_cast<int>(dDegree - 6));
        size_t costChange = abs(static_cast<int>(aDegree - 1 - 6)) + abs(static_cast<int>(bDegree - 1 - 6))
                            + abs(static_cast<int>(cDegree + 1 - 6)) + abs(static_cast<int>(dDegree + 1 - 6));
        if (costChange < costNoChange) {
            // flipping edge does not delete or create edge
            // so we can safely do this while looping through this->m_edges
            this->flipEdge(edge);
        }
    }

    // Compute the centroids for all the vertices.
    std::unordered_map<Face*, Eigen::Vector3f, IndexablePointerHash> faceNormals;
    for (Face* face : this->m_faces) {
        std::array<Vertex*, 3> vertices = this->faceVertices(face);
        // use cross product to find normal
        Eigen::Vector3f normal = (vertices[1]->coords - vertices[0]->coords)
                                     .cross(vertices[2]->coords - vertices[0]->coords)
                                     .normalized();

        faceNormals[face] = normal;
    }
    std::unordered_map<Vertex*, Eigen::Vector3f, IndexablePointerHash> vertexCentroids;
    std::unordered_map<Vertex*, Eigen::Vector3f, IndexablePointerHash> vertexNormals;
    for (Vertex* vertex : this->m_vertices) {
        Eigen::Vector3f centroid = Eigen::Vector3f::Zero();
        Eigen::Vector3f normal = Eigen::Vector3f::Zero();
        auto neighbors = this->vertexHalfEdges(vertex);
        for (HalfEdge* h : neighbors) {
            centroid += h->twin->vertex->coords;
            normal += faceNormals[h->face];
        }
        vertexCentroids[vertex] = centroid / neighbors.size();
        vertexNormals[vertex] = normal / neighbors.size();
    }

    // Move each vertex in the tangent direction toward its centroid.
    for (Vertex* vertex : this->m_vertices) {
        Eigen::Vector3f v = vertexCentroids[vertex] - vertex->coords;
        Eigen::Vector3f n = vertexNormals[vertex];
        v = v - n.dot(v) * n;
        vertex->coords += w * v;
    }
}

IndexablePointerSet<Vertex> &Mesh::getVertices()
{
    return this->m_vertices;
}

IndexablePointerSet<Edge> &Mesh::getEdges()
{
    return this->m_edges;
}

IndexablePointerSet<Face> &Mesh::getFaces()
{
    return this->m_faces;
}

Mesh::Mesh()
{
    this->m_nextIndex = 0;
}

std::tuple<Vertex *, Edge *, Edge *> Mesh::splitEdgeExtraInfo(Edge *edge)
{
    HalfEdge* bc = edge->halfedge;
    HalfEdge* cb = bc->twin;
    HalfEdge* ca = bc->next;
    HalfEdge* ac = ca->twin;
    HalfEdge* ab = ca->next;
    HalfEdge* ba = ab->twin;
    HalfEdge* bd = cb->next;
    HalfEdge* db = bd->twin;
    HalfEdge* dc = bd->next;
    HalfEdge* cd = dc->twin;

    Edge* bcEdge = edge;
    Edge* caEdge = ca->edge;
    Edge* abEdge = ab->edge;
    Edge* dcEdge = dc->edge;
    Edge* bdEdge = bd->edge;

    Vertex* a = ac->vertex;
    Vertex* b = ba->vertex;
    Vertex* c = cd->vertex;
    Vertex* d = db->vertex;

    Face* abc = ab->face;
    Face* bdc = bd->face;

    HalfEdge* ma = this->addHalfEdge();
    HalfEdge* am = this->addHalfEdge();
    HalfEdge* md = this->addHalfEdge();
    HalfEdge* dm = this->addHalfEdge();
    HalfEdge* mb = this->addHalfEdge();
    HalfEdge* bm = this->addHalfEdge();
    HalfEdge* mc = bc;
    HalfEdge* cm = cb;

    Edge* amEdge = this->addEdge();
    Edge* mdEdge = this->addEdge();
    Edge* bmEdge = this->addEdge();
    Edge* mcEdge = bcEdge;
    Vertex* m = this->addVertex();
    Face* abm = this->addFace();
    Face* mbd = this->addFace();
    Face* amc = abc;
    Face* mdc = bdc;

    // ac->next doesn't change
    // ac->edge doesn't change
    // ac->face doesn't change
    // ac->vertex doesn't change
    ca->next = am;
    // ca->edge doesn't change
    ca->face = amc;
    // ca->vertex doesn't change

    am->next = mc;
    am->edge = amEdge;
    am->face = amc;
    am->twin = ma;
    am->vertex = a;
    ma->next = ab;
    ma->edge = amEdge;
    ma->face = abm;
    ma->twin = am;
    ma->vertex = m;

    mc->next = ca;
    mc->edge = mcEdge;
    mc->face = amc;
    mc->twin = cm;
    mc->vertex = m;
    cm->next = md;
    cm->edge = mcEdge;
    cm->face = mdc;
    cm->twin = mc;
    cm->vertex = c;

    md->next = dc;
    md->edge = mdEdge;
    md->face = mdc;
    md->twin = dm;
    md->vertex = m;
    dm->next = mb;
    dm->edge = mdEdge;
    dm->face = mbd;
    dm->twin = md;
    dm->vertex = d;

    // cd->next doesn't change
    // cd->edge doesn't change
    // cd->face doesn't change
    // cd->vertex doesn't change
    dc->next = cm;
    // dc->edge doesn't change
    dc->face = mdc;
    // dc->vertex doesn't change

    // ba->next doesn't change
    // ba->edge doesn't change
    // ba->face doesn't change
    // ba->vertex doesn't change
    ab->next = bm;
    // ab->edge doesn't change
    ab->face = abm;
    // ab->vertex doesn't change

    bm->next = ma;
    bm->edge = bmEdge;
    bm->face = abm;
    bm->twin = mb;
    bm->vertex = b;
    mb->next = bd;
    mb->edge = bmEdge;
    mb->face = mbd;
    mb->twin = bm;
    mb->vertex = m;

    // db->next doesn't change
    // db->edge doesn't change
    // db->face doesn't change
    // db->vertex doesn't change
    bd->next = dm;
    // bd->edge doesn't change
    bd->face = mbd;
    // bd->vertex doesn't change

    caEdge->halfedge = ca;
    amEdge->halfedge = am;
    mcEdge->halfedge = mc;
    mdEdge->halfedge = md;
    dcEdge->halfedge = dc;
    abEdge->halfedge = ab;
    bmEdge->halfedge = bm;
    bdEdge->halfedge = bd;

    a->halfedge = ab;
    b->halfedge = bd;
    c->halfedge = ca;
    d->halfedge = dc;
    m->halfedge = mc;

    abm->halfedge = ab;
    mbd->halfedge = mb;
    amc->halfedge = am;
    mdc->halfedge = md;
    m->coords = (b->coords + c->coords) / 2.f;
    return std::make_tuple(m, amEdge, mdEdge);
}

std::tuple<Vertex *, Vertex *, Edge *, Edge *, Edge *, Edge *, Edge *> Mesh::collapseEdgeExtraInfo(Edge *edge)
{
    HalfEdge* cd = edge->halfedge;
    Vertex* c = cd->vertex;
    HalfEdge* dc = cd->twin;
    Vertex* d = dc->vertex;

    IndexablePointerSet<HalfEdge> cHalfEdges = this->vertexHalfEdges(cd->vertex);
    unordered_map<Vertex*, HalfEdge*, IndexablePointerHash> cNeighbors;
    for (HalfEdge* cHalfEdge : cHalfEdges) {
        cNeighbors[cHalfEdge->twin->vertex] = cHalfEdge;
    }
    IndexablePointerSet<HalfEdge> dHalfEdges = this->vertexHalfEdges(dc->vertex);

    Vertex* a = nullptr;
    Vertex* b = nullptr;
    HalfEdge* ca;
    HalfEdge* ac;
    HalfEdge* da;
    HalfEdge* ad;
    HalfEdge* cb;
    HalfEdge* bc;
    HalfEdge* bd;
    HalfEdge* db;

    std::tuple<Vertex*, Vertex*, Edge*, Edge*, Edge*, Edge*, Edge*> nullValue = std::make_tuple(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);

    for (HalfEdge* dHalfEdge : dHalfEdges) {
        Vertex* dVertex = dHalfEdge->twin->vertex;
        auto cxHalfEdgeIt = cNeighbors.find(dVertex);
        if (cxHalfEdgeIt == cNeighbors.end()) {
            continue;
        }
        if (this->degree(dVertex) == 3) {
            return nullValue;
        }
        auto [_, cxHalfEdge] = *cxHalfEdgeIt;
        if (a == nullptr) {
            a = dVertex;
            ca = cxHalfEdge;
            da = dHalfEdge;
        } else if (b == nullptr) {
            b = dVertex;
            cb = cxHalfEdge;
            db = dHalfEdge;
        } else {
            return nullValue;
        }
    }

    // if we made it here, we must have two shared neighbors
    assert(a != nullptr);
    assert(b != nullptr);
    assert(ca != nullptr);
    assert(da != nullptr);
    assert(cb != nullptr);
    assert(db != nullptr);

    // swap if necessary to canonicalize naming
    if (cd->next->next->vertex == b) {
        swap(a, b);
        swap(ca, cb);
        swap(da, db);
    }
    ac = ca->twin;
    ad = da->twin;
    bc = cb->twin;
    bd = db->twin;

    Edge* cdEdge = edge;
    Edge* acEdge = ac->edge;
    Edge* daEdge = da->edge;
    Edge* cbEdge = cb->edge;
    Edge* bdEdge = bd->edge;
    Face* cbd = cb->face;
    Face* acd = ac->face;

    Vertex* m = c;
    m->coords = (c->coords + d->coords) / 2;
    Edge* amEdge = acEdge;
    Edge* mbEdge = cbEdge;
    HalfEdge* ma = ca;
    HalfEdge* am = ac;
    HalfEdge* mb = cb;
    HalfEdge* bm = bc;

    ma->next = ca->next;
    ma->edge = amEdge;
    ma->face = ca->face;
    ma->vertex = m;

    am->next = ad->next;
    am->edge = amEdge;
    am->face = ad->face;
    am->vertex = a;

    mb->next = db->next;
    mb->edge = mbEdge;
    mb->face = db->face;
    mb->vertex = m;

    bm->next = bc->next;
    bm->edge = mbEdge;
    bm->face = bc->face;
    bm->vertex = b;

    amEdge->halfedge = am;
    mbEdge->halfedge = mb;

    a->halfedge = am;
    m->halfedge = ma;
    b->halfedge = bm;

    // do all the triangles on the adb
    HalfEdge* dbSide = db;
    while (dbSide != da) {
        dbSide->vertex = m;

        dbSide = dbSide->next->next->twin;
    }
    ad->face->halfedge = am;
    ad->next->next->next = am;
    db->face->halfedge = mb;
    db->next->next->next = mb;

    this->removeHalfEdge(da);
    this->removeHalfEdge(ad);
    this->removeHalfEdge(bd);
    this->removeHalfEdge(db);
    this->removeHalfEdge(cd);
    this->removeHalfEdge(dc);
    this->removeFace(acd);
    this->removeFace(cbd);
    return std::make_tuple(m, d, acEdge, cbEdge, cdEdge, daEdge, bdEdge);
}

HalfEdge* Mesh::addHalfEdge()
{
    HalfEdge* halfedge = new HalfEdge {};
    halfedge->index = this->m_nextIndex++;
    this->m_halfedges.insert(halfedge);
    return halfedge;
}

Vertex* Mesh::addVertex()
{
    Vertex* vertex = new Vertex {};
    vertex->index = this->m_nextIndex++;
    this->m_vertices.insert(vertex);
    return vertex;
}

Edge* Mesh::addEdge()
{
    Edge* edge = new Edge {};
    edge->index = this->m_nextIndex++;
    this->m_edges.insert(edge);
    return edge;
}

Face* Mesh::addFace()
{
    Face* face = new Face {};
    face->index = this->m_nextIndex++;
    this->m_faces.insert(face);
    return face;
}

void Mesh::removeHalfEdge(HalfEdge *halfedge)
{
    this->m_halfedges.erase(halfedge);
    delete halfedge;
}

void Mesh::removeVertex(Vertex *vertex)
{
    this->m_vertices.erase(vertex);
    delete vertex;
}

void Mesh::removeEdge(Edge *edge)
{
    this->m_edges.erase(edge);
    delete edge;
}

void Mesh::removeFace(Face *face)
{
    this->m_faces.erase(face);
    delete face;
}


