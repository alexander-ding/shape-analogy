#include "analogy.h"

Analogy::Analogy(Mesh aPrime, Mesh b)
    : m_aPrime(aPrime), m_b(b), m_bPrime(b)
{

}

void Analogy::computeBPrime()
{
    this->m_bTargetNormals = this->computeTargetNormals();
}

MatrixXf Analogy::computeTargetNormals()
{
    // n x 3
    Eigen::MatrixXf aPrimeNormals = this->m_aPrime.computeFaceNormals().transpose();
    // 3 x m
    Eigen::MatrixXf bNormals = this->m_b.computeVertexNormals();
    // n x m
    Eigen::MatrixXf similarities = aPrimeNormals * bNormals;
    // m x n
    for (size_t i = 0; i < bNormals.cols(); i++) {
        size_t maxIndex;
        similarities.col(i).maxCoeff(&maxIndex);
        bNormals.col(i) = aPrimeNormals.row(maxIndex);
    }
    return bNormals;
}
