#include "arap.h"
#include "mesh.h"

#include <iostream>
#include <unordered_set>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <limits>
#include <QtConcurrent>

using namespace std;
using namespace Eigen;



void ARAP::build()
{
    // std::cout << "Building L's" << std::endl;
    const Eigen::MatrixXf& vertices = this->m_mesh.getVertices();

    // zero out adjacency sets
    const Eigen::MatrixXi& faces = this->m_mesh.getFaces();
    for (size_t i = 0; i < m_vertexEdges.size(); i++) {
        for (auto it = this->m_vertexEdges[i].begin(); it != this->m_vertexEdges[i].end(); it++) {
            it->second = 0.f;
        }
    }

    // compute adjacency sets
    for (size_t i = 0; i < faces.cols(); i++) {
        for (int j = 0; j < 3; j++) {
            int v = faces.col(i)[j];
            int edgeA = faces.col(i)[(j + 1) % 3];
            int edgeB = faces.col(i)[(j + 2) % 3];
            Vector3f a = vertices.col(edgeA) - vertices.col(v);
            Vector3f b = vertices.col(edgeB) - vertices.col(v);

            float edgeWeight = std::abs(a.dot(b) / a.cross(b).norm() / 2.f);
            assert(edgeWeight >= 0);
            this->m_vertexEdges[edgeA].emplace(edgeB, 0.f);
            this->m_vertexEdges[edgeB].emplace(edgeA, 0.f);
            this->m_vertexEdges[edgeA][edgeB] += edgeWeight;
            this->m_vertexEdges[edgeB][edgeA] += edgeWeight;
        }
    }

    // compute L
    SparseMatrix<float> L(vertices.cols(), vertices.cols());
    for (size_t i = 0; i < vertices.cols(); i++) {
        for (auto [j, weight] : this->m_vertexEdges[i]) {
            L.coeffRef(i, i) += weight;
            L.coeffRef(i, j) -= weight;
        }
    }

    // apply user constraints to L
    for (int anchor : this->m_anchors) {
        L.row(anchor) *= 0;
        L.col(anchor) *= 0;
        L.coeffRef(anchor, anchor) = 1;
    }

    // precompute the decomposition
    this->m_solver.compute(L);
    assert(this->m_solver.info()==Success);
    this->m_isValid = true;
}

ARAP::ARAP(Mesh mesh) : m_vertexEdges(), m_mesh(mesh), m_prevFrameVertices(mesh.getVertices()) {}


void ARAP::init(Eigen::Vector3f &coeffMin, Eigen::Vector3f &coeffMax)
{
    this->m_isValid = false;
    this->m_vertexEdges.resize(this->m_mesh.getVertices().cols(), std::unordered_map<int, float>());

    MatrixXf& vertices = this->m_mesh.getVertices();
    coeffMin = vertices.rowwise().minCoeff();
    coeffMax = vertices.rowwise().maxCoeff();
}

bool ARAP::invalidate()
{
    bool isValid = this->m_isValid;
    this->m_isValid = false;
    return isValid;
}

inline QVector<int> range(int start, int end)
{
    QVector<int> l(end-start);
    std::iota(l.begin(), l.end(), start);
    return l;
}


// Move an anchored vertex, defined by its index, to targetPosition
void ARAP::move(int vertex, Vector3f targetPosition)
{
    if (!this->m_isValid) this->build();
    // std::cout << "Processing move" << std::endl;

    MatrixXf newVertices = this->m_mesh.getVertices();
    std::vector<Matrix3f> Rs(newVertices.cols(), Matrix3f::Identity());
    MatrixXf b(newVertices.cols(), 3);

    MatrixXf prevVertices;
    float prevError = std::numeric_limits<float>::max();

    QVector<int> vertexRange = range(0, newVertices.cols());
    while (true) {
        prevVertices = newVertices;
        // solve new R's
        // #pragma omp parallel for
        b.setZero();

        QtConcurrent::blockingMap(vertexRange, [&](int i) {
            // for (size_t i = 0; i < this->m_vertices.cols(); i++) {
            Matrix3f S = Matrix3f::Zero();
            for (auto [j, weight] : this->m_vertexEdges[i]) {
                S += weight * (this->m_mesh.getVertices().col(i) - this->m_mesh.getVertices().col(j)) * (newVertices.col(i) - newVertices.col(j)).transpose();
            }
            auto SVD = S.jacobiSvd(ComputeFullU | ComputeFullV);
            Matrix3f R = SVD.matrixV() * SVD.matrixU().transpose();

            if (R.determinant() < 0) {
                // std::cout << "Negative:" << R.determinant() << std::endl;
                size_t smallestCol = 0;
                for (size_t j = 1; j < 3; j++) {
                    if (SVD.singularValues()[j] < SVD.singularValues()[smallestCol]) {
                        smallestCol = j;
                    }
                }
                Matrix3f U(SVD.matrixU());
                U.col(smallestCol) = -U.col(smallestCol);
                R = SVD.matrixV() * U.transpose();
                assert(R.determinant() >= 0);
            }
            Rs[i] = R;
            // }
        });

        QtConcurrent::blockingMap(vertexRange, [&](int i) {
            for (auto [j, weight] : this->m_vertexEdges[i]) {
                b.row(i) += weight / 2.f * ((Rs[i] + Rs[j]) * (this->m_mesh.getVertices().col(i) - this->m_mesh.getVertices().col(j)));
            }
        });

        // sync anchor positions
        newVertices.col(vertex) = targetPosition;

        // account for zero-ed out columns
        for (int anchor : this->m_anchors) {
            for (auto [j, weight] : this->m_vertexEdges[anchor]) {
                b.row(j) += weight * newVertices.col(anchor);
            }
        }

        // set explicit positions for anchors
        for (int anchor : this->m_anchors) {
            b.row(anchor) = newVertices.col(anchor);
        }

        newVertices = this->m_solver.solve(b).transpose();
        assert(this->m_solver.info()==Success);

        float error = (newVertices - prevVertices).norm();
        // std::cout << error << std::endl;
        // threshold
        if (error < 1e-2 || prevError - error < 5e-4) {
            break;
        }
        prevError = error;
    }
    this->m_mesh.setVertices(newVertices);
}

SelectMode ARAP::select(int closest_vertex)
{
    if (closest_vertex == -1) return SelectMode::None;

    bool vertexIsNowSelected = m_anchors.find(closest_vertex) == m_anchors.end();

    if (vertexIsNowSelected) {
        m_anchors.insert(closest_vertex);
    } else {
        m_anchors.erase(closest_vertex);
    }

    return vertexIsNowSelected ? SelectMode::Anchor : SelectMode::Unanchor;
}

bool ARAP::selectWithSpecifiedMode(int closest_vertex, SelectMode mode)
{
    if (closest_vertex == -1) return false;
    switch (mode) {
    case SelectMode::None: {
        return false;
    }
    case SelectMode::Anchor: {
        if (m_anchors.find(closest_vertex) != m_anchors.end()) return false;
        m_anchors.insert(closest_vertex);
        break;
    }
    case SelectMode::Unanchor: {
        if (m_anchors.find(closest_vertex) == m_anchors.end()) return false;
        m_anchors.erase(closest_vertex);
        break;
    }
    }

    return true;
}

int ARAP::getClosestVertex(Vector3f start, Vector3f ray, float threshold)
{
    int closest_vertex = -1;
    int i = 0;
    float dist = numeric_limits<float>::max();
    ParametrizedLine line = ParametrizedLine<float, 3>::Through(start, start + ray);

    MatrixXf& vertices = this->m_mesh.getVertices();
    for (size_t i = 0; i < vertices.cols(); i++) {
        float d = line.distance(vertices.col(i));
        if (d<dist) {
            dist = d;
            closest_vertex = i;
        }
    }

    if (dist >= threshold) closest_vertex = -1;
    // std::cout << closest_vertex << std::endl;

    return closest_vertex;
}

bool ARAP::getAnchorPos(int lastSelected,
                                 Eigen::Vector3f& pos,
                                 Eigen::Vector3f  ray,
                                 Eigen::Vector3f  start)
{
    bool isAnchor = m_anchors.find(lastSelected) != m_anchors.end();
    if (isAnchor) {
        MatrixXf& vertices = this->m_mesh.getVertices();
        Eigen::Vector3f oldPos = vertices.col(lastSelected);
        Eigen::ParametrizedLine line = ParametrizedLine<float, 3>::Through(start, start+ray);
        pos = line.projection(oldPos);
    }
    return isAnchor;
}

void ARAP::undo() {
    this->m_mesh.setVertices(this->m_prevFrameVertices);
}
