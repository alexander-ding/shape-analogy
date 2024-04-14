#ifndef ANALOGY_H
#define ANALOGY_H

#include "mesh.h"

class Analogy
{
public:
    Analogy(Mesh aPrime, Mesh b);

    Mesh& getAPrime() {return this->m_aPrime;}
    Mesh& getB() {return this->m_b;}
    Mesh& getBPrime() {return this->m_bPrime;}
    Eigen::MatrixXf& getBTargetNormals() { return this->m_bTargetNormals; }
    void computeBPrime();
private:
    Eigen::MatrixXf computeTargetNormals();

    Eigen::MatrixXf m_bTargetNormals;
    Mesh m_bPrime;
    Mesh m_aPrime;
    Mesh m_b;
};

#endif // ANALOGY_H
