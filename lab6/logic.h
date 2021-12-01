#ifndef LOGIC_H
#define LOGIC_H

#include <vector>
#include <map>
#include <random>
#include <QVector>
#include <qcustomplot.h>

class MinErrorClassifier {
    public:
        QVector<double> P;
        QVector<double> train_classes;
        QVector<double> test_classes;
        QVector<QVector<double>> train_data;
        QVector<QVector<double>> test_data;

        int c, d;
        int N_train, N_test;
        double error;

        QVector<double> classes_count;
        QVector<QVector<QVector<double>>> sigma_matrices;
        QVector<QVector<QVector<double>>> sigma_inv_matrices;
        QVector<QVector<QVector<double>>> A_matrices;
        QVector<QVector<double>> transformation_matrix;
        QVector<QVector<double>> mu_vectors;
        std::map<int, QVector<int>> test_data_by_class;

        MinErrorClassifier() = default;
        void Train();
        void Classify();
        void SetC(int new_c);
        void SetD(int new_d);
        void SetN(int new_N);
        void ClearData();

    private:
        void CalculateCoefMatrix();
        void CalculateSigmaInvMatrices();
        double d_function(int class_number, int vector_number);
};

class ParamClassifier {
    public:
        QVector<double> P;
        QVector<double> train_classes;
        QVector<double> test_classes;
        QVector<QVector<double>> train_data;
        QVector<QVector<double>> test_data;

        int c, d;
        int N_train, N_test;

        QVector<double> classes_count;
        QVector<QVector<QVector<double>>> sigma_matrices;
        QVector<QVector<QVector<double>>> sigma_inv_matrices;
        QVector<QVector<QVector<double>>> A_matrices;
        QVector<QVector<double>> transformation_matrix;
        QVector<QVector<double>> mu_vectors;
//        QVector<QVector<double>> var_vectors;

        std::map<int, QVector<int>> test_data_by_class;
        double error;

        ParamClassifier() = default;
//        ParamClassifier(ParamClassisfier& other) = default;
//        ParamClassifier& operator=(ParamClassifier& other);
        ParamClassifier(
                int c, int d,
                int N_train, int N_test,
                QVector<double>& P,
                QVector<double>& train_classes,
                QVector<double>& test_classes,
                QVector<QVector<double>>& train_data,
                QVector<QVector<double>>& test_data
                );
        void Train();
        void Classify();
        void SetC(int new_c);
        void SetD(int new_d);
        void SetN(int new_N);
        void ClearData();

    private:
        void CalculateCoefMatrix();
        void CalculateSigmaInvMatrices();
        double d_function(int class_number, int vector_number);
};

class NeighbourClassifier {
    public:
        QVector<double> P;
        QVector<double> train_classes;
        QVector<double> test_classes;
        QVector<QVector<double>> train_data;
        QVector<QVector<double>> test_data;

        int c, d;
        int N_train, N_test;

        QVector<double> classes_count;
        QVector<QVector<double>> transformation_matrix;
//        QVector<QVector<double>> var_vectors;

        std::map<int, QVector<int>> test_data_by_class;
        double error;

        NeighbourClassifier() = default;

        void Classify();
        void SetC(int new_c);
        void SetD(int new_d);
        void ClearData();

        double d_function(int class_number, int vector_number);
        double f(int vector_number,int class_number);
//        double f(QVector<double>& u_);
};


class ParzenClassifier {
    public:
        QVector<double> P;
        QVector<double> train_classes;
        QVector<double> test_classes;
        QVector<QVector<double>> train_data;
        QVector<QVector<double>> test_data;

        int c, d;
        int N_train, N_test;

        QVector<double> classes_count;
        QVector<QVector<double>> transformation_matrix;
//        QVector<QVector<double>> var_vectors;

        std::map<int, QVector<int>> test_data_by_class;
        double error;

        QVector<double> v;
        ParzenClassifier() = default;

        void Classify();
        void SetC(int new_c);
        void SetD(int new_d);
        void ClearData();

        double d_function(int class_number, int vector_number);
        double Parzen(int vector_number,int class_number);
        double u(QVector<double>& u_);
};



class logic
{
    public:
        logic();
        void Model();

        void SetC(int new_c);
        void SetD(int new_d);
        void SetN(int new_N_train, int new_N_test);
        void ClearData();
        void GenerateClasses();
        void ModelNormalDistributionVector(QVector<double>& result);
        void CalculateRandomVector(int class_number, QVector<double>& result);
        void CalculateCoefMatrix();
        void GetPointsData(int class_number, QVector<double>& x, QVector<double>& y, std::string type);
        void CalculateSigmaInvMatrices();
        double LikehoodFunction(int i, int class_number);
        double CalculateErrorProbability();
        double AposterioriProbability(int class_number, int vector_number);
        double g(int class_number, QVector<double>& x);
        double d_function(int class_number, int vector_numer);
        void CreateClassifiers();

        int c, d;
        int axis1, axis2;
        int N_train, N_test;
        QVector<double> P;
        QVector<double> test_classes;
        QVector<double> train_classes;
        QVector<QVector<QVector<double>>> sigma_matrices;
        QVector<QVector<QVector<double>>> sigma_inv_matrices;
        QVector<QVector<QVector<double>>> A_matrices;
        QVector<QVector<double>> transformation_matrix;
        QVector<QVector<double>> mu_vectors;
        QVector<QVector<double>> var_vectors;
        QVector<QVector<double>> train_data;
        QVector<QVector<double>> test_data;

        ParamClassifier* param_classifier;
        MinErrorClassifier* bayes_classifier;
        ParzenClassifier* parzen_classifier;
        NeighbourClassifier* knn_classifier;

        // classif. 1 (P(error)->min)
        // пары типа (номер класса; номера точкек в modelled_data, которые отнесли к этому классу)
        void ClassifyMinError();
        std::map<int, QVector<int>> c1_data_by_class;

        std::mt19937 gen;
        std::normal_distribution<> distr;
};

#endif // LOGIC_H
