#pragma once

#include "mesh.h"
#include <unordered_map>
#include <unordered_set>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>

enum SelectMode
{
    None     = 0,
    Anchor   = 1,
    Unanchor = 2
};

class ARAP
{
private:
    bool m_isValid;
    Mesh m_mesh;
    std::unordered_set<int> m_anchors;
    std::vector<std::unordered_map<int, float>> m_vertexEdges;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> m_solver;

    void build();
public:
    ARAP(Mesh mesh);

    Mesh& getMesh() { return this->m_mesh; }
    std::unordered_set<int>& getAnchors() { return this->m_anchors; }
    void init(Eigen::Vector3f &min, Eigen::Vector3f &max);
    bool invalidate();
    void move(int vertex, Eigen::Vector3f pos);

    // anchor
    SelectMode select(int closest_vertex);
    int getClosestVertex(Vector3f start, Vector3f ray, float threshold);
    bool selectWithSpecifiedMode(int closest_vertex, SelectMode mode);
    bool getAnchorPos(int lastSelected,
                            Eigen::Vector3f& pos,
                            Eigen::Vector3f  ray,
                            Eigen::Vector3f  start);
};
