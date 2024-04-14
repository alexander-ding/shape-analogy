#pragma once

#include <Eigen/Dense>
#include <QString>

using namespace std;
using namespace Eigen;

class Mesh
{
public:
    Mesh(const MatrixXf &vertices, const MatrixXi &faces, bool center = true);
    Mesh(const QString &path, bool center = true);
    ~Mesh();

    void saveToFile(const QString& path);

    MatrixXf &getVertices()
    {
        return this->m_vertices;
    }

    void setVertices(const Eigen::MatrixXf &vertices);

    MatrixXi &getFaces()
    {
        return this->m_faces;
    }

    MatrixXf computeFaceNormals();
    MatrixXf computeVertexNormals();

private:
    void init(const MatrixXf &vertices, const MatrixXi &faces);
    void center();

    MatrixXf m_vertices;
    MatrixXi m_faces;
};
