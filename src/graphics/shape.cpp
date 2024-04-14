#include "shape.h"

#include "graphics/shader.h"

using namespace Eigen;
using namespace std;

Shape::Shape()
    : m_modelMatrix(Eigen::Matrix4f::Identity())
{
}


void Shape::init(const Eigen::MatrixXf &vertices, const Eigen::MatrixXf &normals, const Eigen::MatrixXi &triangles, const Eigen::MatrixXf &colors)
{
    glGenBuffers(1, &m_surfaceVbo);
    glGenBuffers(1, &m_surfaceIbo);
    glGenVertexArrays(1, &m_surfaceVao);

    glBindVertexArray(m_surfaceVao);

    // Upload Vertex Buffer Data
    glBindBuffer(GL_ARRAY_BUFFER, m_surfaceVbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (vertices.size() + normals.size() + colors.size()), nullptr, GL_DYNAMIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * vertices.size(), vertices.data());
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), sizeof(float) * normals.size(), normals.data());
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(float) * (vertices.size() + normals.size()), sizeof(float) * colors.size(), colors.data());

    // Set Vertex Attribute Pointers
    glEnableVertexAttribArray(0); // Vertex positions
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(1); // Normal vectors
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)(sizeof(float) * vertices.size()));
    glEnableVertexAttribArray(2); // Colors
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (void*)(sizeof(float) * (vertices.size() + normals.size())));


    // Upload Index Buffer Data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_surfaceIbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * triangles.size(), triangles.data(), GL_STATIC_DRAW);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    m_nTriangles = triangles.cols();
}

void Shape::setVertices(const Eigen::MatrixXf& vertices, const Eigen::MatrixXf& normals)
{
    glBindBuffer(GL_ARRAY_BUFFER, m_surfaceVbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * vertices.size(), vertices.data());
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), sizeof(float) * normals.size(), normals.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void Shape::draw(Shader *shader)
{
    Eigen::Matrix3f m3 = m_modelMatrix.topLeftCorner(3, 3);
    Eigen::Matrix3f inverseTransposeModel = m3.inverse().transpose();


    shader->setUniform("model", m_modelMatrix);
    shader->setUniform("inverseTransposeModel", inverseTransposeModel);
    glBindVertexArray(m_surfaceVao);
    glDrawElements(GL_TRIANGLES, m_nTriangles * 3, GL_UNSIGNED_INT, reinterpret_cast<GLvoid *>(0));
    glBindVertexArray(0);
}

