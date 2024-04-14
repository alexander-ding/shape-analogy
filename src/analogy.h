#ifndef ANALOGY_H
#define ANALOGY_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <vector>
#include <unordered_map>

#include "mesh.h"

class Analogy
{
public:
    Analogy(Mesh aPrime, Mesh b);

    Mesh& getAPrime() {return this->m_aPrime;}
    Mesh& getB() {return this->m_b;}
    Mesh& getBPrime() {return this->m_bPrime;}
    Eigen::MatrixXf& getBTargetNormals() { return this->m_bTargetNormals; }
    void computeBPrime(float lambda = 1);
private:
    Eigen::MatrixXf computeTargetNormals();
    Eigen::VectorXf computeVoronoiAreas();
    std::vector<std::unordered_map<size_t, float>> computeEdgeWeights();
    std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>>> buildL(const std::vector<std::unordered_map<size_t, float> >& edgeWeights);

    Eigen::MatrixXf m_bTargetNormals;
    Mesh m_bPrime;
    Mesh m_aPrime;
    Mesh m_b;
};

#endif // ANALOGY_H
