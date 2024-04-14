#ifndef SHAPE_H
#define SHAPE_H

#include <GL/glew.h>

#include <Eigen/Dense>

class Shader;

class Shape
{
public:
    Shape();

    inline void init(const Eigen::MatrixXf &vertices, const Eigen::MatrixXf &normals, const Eigen::MatrixXi &triangles) {
        const Eigen::MatrixXf colors = Eigen::MatrixXf::Ones(vertices.rows(), vertices.cols()).array() / 2.f;
        this->init(vertices, normals, triangles, colors);
    };
    void init(const Eigen::MatrixXf &vertices, const Eigen::MatrixXf &normals, const Eigen::MatrixXi &triangles, const Eigen::MatrixXf &colors);

    void setVertices(const Eigen::MatrixXf& vertices, const Eigen::MatrixXf& normals);

    Eigen::Matrix4f& getModelMatrix() {
        return this->m_modelMatrix;
    };
    void setModelMatrix(const Eigen::Affine3f &model) { m_modelMatrix = model.matrix(); };

    void draw(Shader *shader);

private:
    GLuint m_surfaceVao;
    GLuint m_surfaceVbo;
    GLuint m_surfaceIbo;

    Eigen::Matrix4f m_modelMatrix;

    size_t m_nTriangles;
};

#endif // SHAPE_H
