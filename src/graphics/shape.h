#ifndef SHAPE_H
#define SHAPE_H

#include <GL/glew.h>
#include <vector>

#include <Eigen/Dense>

class Shader;

class Shape
{
public:
    Shape();

    void init(const Eigen::MatrixXf &vertices, const Eigen::MatrixXf &normals, const Eigen::MatrixXi &triangles, const Eigen::MatrixXf &colors);

    void init(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3f> &normals, const std::vector<Eigen::Vector3i> &triangles);
    void init(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3i> &triangles);
    void init(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3i> &triangles, const std::vector<Eigen::Vector4i> &tetIndices);

    void setVertices(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3f> &normals);
    void setVertices(const std::vector<Eigen::Vector3f>& vertices);

    std::vector<Eigen::Vector3i>& getFaces() {
        return this->m_faces;
    }

    Eigen::Matrix4f& getModelMatrix() {
        return this->m_modelMatrix;
    };
    void setModelMatrix(const Eigen::Affine3f &model);

    void toggleWireframe();

    void draw(Shader *shader);

private:
    GLuint m_surfaceVao;
    GLuint m_tetVao;
    GLuint m_surfaceVbo;
    GLuint m_tetVbo;
    GLuint m_surfaceIbo;
    GLuint m_tetIbo;

    unsigned int m_numSurfaceVertices;
    unsigned int m_numTetVertices;
    unsigned int m_verticesSize;
    float m_red;
    float m_blue;
    float m_green;
    float m_alpha;

    std::vector<Eigen::Vector3i> m_faces;

    Eigen::Matrix4f m_modelMatrix;
};

#endif // SHAPE_H
