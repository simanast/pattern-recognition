#ifndef LOGIC_H
#define LOGIC_H

#include <vector>
#include <map>
#include <random>
#include <QVector>
#include <qcustomplot.h>

class logic
{
    public:
        logic();
        void Model();

        void SetC(int new_c);
        void SetD(int new_d);
        void SetN(int new_N);
        void ClearData();
        void GenerateClasses();
        void ModelNormalDistributionVector(QVector<double>& result);
        void CalculateRandomVector(int class_number, QVector<double>& result);
        void CalculateCoefMatrix();
        void GetPointsData(int class_number, QVector<double>& x, QVector<double>& y, std::string type);
        void CalculateSigmaInvMatrices();
        double LikehoodFunction(int i, int class_number);
        double CalculateErrorProbability();
        double CalculateMeanRisk();
        double AposterioriProbability(int class_number, int vector_number);
        double g(int class_number, QVector<double>& x);
        double d_function(int class_number, int vector_numer);
        double risk(int class_number, int vector_number);

        int c, d;
        int axis1, axis2;
        int N;
        QVector<double> P;
        QVector<double> classes;
        QVector<QVector<double>> lambdas;
        QVector<QVector<QVector<double>>> sigma_matrices;
        QVector<QVector<QVector<double>>> sigma_inv_matrices;
        QVector<QVector<QVector<double>>> A_matrices;
        QVector<QVector<double>> risk_transformation_matrix;
        QVector<QVector<double>> transformation_matrix;
        QVector<QVector<double>> mu_vectors;
        QVector<QVector<double>> var_vectors;
        QVector<QVector<double>> modelled_data;
        double sum_risk = 0;

        // classif. 1 (P(error)->min)
        // пары типа (номер класса; номера точкек в modelled_data, которые отнесли к этому классу)
        void ClassifyMinError();
        std::map<int, QVector<int>> c1_data_by_class;

        // classif. 2 (E[Risk]->min)
        void ClassifyMinRisk();
        std::map<int, QVector<int>> c2_data_by_class;

        std::mt19937 gen;
        std::normal_distribution<> distr;
};

#endif // LOGIC_H
