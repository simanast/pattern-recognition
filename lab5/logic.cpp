#include "logic.h"
#include <iostream>
#include <iomanip>

QVector<double> Multiply(QVector<QVector<double>>& lhs, QVector<double>& rhs) {
    // m x m * m x 1 => m x 1
    QVector<double> result(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        for (int j = 0; j < rhs.size(); j++) {
            result[i] += lhs[i][j] * rhs[j];
        }
    }
    return result;
}

QVector<double> Multiply(QVector<double> lhs, QVector<QVector<double>> rhs) {
    // 1 x m * m x m => 1 x m
    QVector<double> result(lhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        for (int j = 0; j < lhs.size(); j++) {
            result[i] += lhs[j] * rhs[j][i];
        }
    }
    return result;
}

double Multiply(QVector<double> lhs, QVector<double> rhs) {
    // 1 x m * m x m => 1 x m
    double result = 0;
    for (int i = 0; i < rhs.size(); i++) {
        result += lhs[i] * rhs[i];
    }
    return result;
}

QVector<QVector<double>> Multiply(QVector<QVector<double>>& lhs, QVector<QVector<double>>& rhs) {
    QVector<QVector<double>> result(lhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        result[i].resize(lhs.size());
        for (int j = 0; j < lhs.size(); j++) {
            result[i][j] += lhs[i][j] * rhs[j][i];
        }
    }
    return result;
}

void Multiply(QVector<QVector<double>>& lhs, double rhs) {
    for (int i = 0; i < lhs.size(); i++) {
        for (int j = 0; j < lhs.size(); j++) {
            lhs[i][j] *= rhs;
        }
    }
}

QVector<QVector<double>> scalar_multiply(QVector<QVector<double>>& lhs, double rhs) {
    QVector<QVector<double>> result(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        result[i].resize(lhs.size());
        for (int j = 0; j < lhs.size(); j++) {
            result[i][j] = lhs[i][j] * rhs;
        }
    }
    return result;
}

QVector<double> Add(QVector<double>& lhs, QVector<double>& rhs) {
    QVector<double> result(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

QVector<double> Subtract(QVector<double>& lhs, QVector<double>& rhs) {
    QVector<double> result(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

void T(QVector<QVector<double>>& A0) {
    //транспонирование
    size_t n = A0.size();
    for (size_t row = 0; row < n - 1; row++) {
        for (size_t col = row + 1; col < n; col++)
            std::swap(A0[col][row], A0[row][col]);
    }
}

double det(QVector<QVector<double>> matrix) {
    int n = matrix.size();
    QVector<QVector<double>> B = matrix;
    //приведение матрицы к верхнетреугольному виду
    for (int step = 0; step < n - 1; step++)
        for (int row = step + 1; row < n; row++) {
            double coeff = -B[row][step] / B[step][step]; //метод Гаусса
            for(int col = step; col < n; col++)
                B[row][col] += B[step][col] * coeff;
        }
    //Рассчитать определитель как произведение элементов главной диагонали
    double Det = 1;
    for (int i = 0; i < n; i++)
        Det *= B[i][i];

    return Det;
}

QVector<QVector<double>> InvMatrix(QVector<QVector<double>>& A) {
    auto matrix = A;
    QVector<QVector<double>> result;
    int size = matrix.size();
    QVector<QVector<double>> E(size, QVector<double>(size));
    //Заполнение единичной матрицы
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if(i == j) E[i][j] = 1.0;
            else E[i][j] = 0.0;
        }
    }
    //Трансформация исходной матрицы в верхнетреугольную
    for (int k = 0; k < size; k++) {
        if(abs(matrix[k][k]) < 1e-8) {
            bool changed = false;
            for (int i = k + 1; i < size; i++) {
                if (abs(matrix[i][k]) > 1e-8) {
                    std::swap(matrix[k], matrix[i]);
                    std::swap(E[k], E[i]);
                    changed = true;
                    break;
                }
            }
            if(!changed)
            return QVector<QVector<double>>();
        }

        double div = matrix[k][k];

        for (int j = 0; j < size; j++) {
            matrix[k][j] /= div;
            E[k][j] /= div;
        }

        for (int i = k + 1; i < size; i++) {
            double multi = matrix[i][k];
            for (int j = 0; j < size; j++) {
                matrix[i][j] -= multi * matrix[k][j];
                E[i][j] -= multi * E[k][j];
            }
        }
    }

    //Формирование единичной матрицы из исходной
    //и обратной из единичной
    for (int k = size - 1; k > 0; k--) {
        for (int i = k - 1; i > -1; i--) {
            double multi = matrix[i][k];

            for (int j = 0; j < size; j++) {
                matrix[i][j] -= multi * matrix[k][j];
                E[i][j] -= multi * E[k][j];
            }
        }
    }
    result = E;
    return result;
}

void logic::SetC(int new_c) {
    c = new_c;
    sigma_matrices.resize(c);
    sigma_inv_matrices.resize(c);
    A_matrices.resize(c);

    mu_vectors.resize(c);
    var_vectors.resize(c);

    P.resize(c);
    transformation_matrix.resize(c);
    risk_transformation_matrix.resize(c);
    lambdas.resize(c);
    for (int i = 0; i < c; i++) {
        transformation_matrix[i].resize(c);
        risk_transformation_matrix[i].resize(c);
        lambdas[i].resize(c);
    }
}

void logic::SetD(int new_d) {
    d = new_d;
    for (int class_number = 0; class_number < c; class_number++) {
        sigma_matrices[class_number].resize(d);
        sigma_inv_matrices[class_number].resize(d);
        A_matrices[class_number].resize(d);
        for (int i = 0; i < d; i++) {
            sigma_matrices[class_number][i].resize(d);
            sigma_inv_matrices[class_number][i].resize(d);
            A_matrices[class_number][i].resize(d);
        }
    }
    for (int class_number = 0; class_number < c; class_number++) {
        mu_vectors[class_number].resize(d);
        var_vectors[class_number].resize(d);
    }
}

void logic::SetN(int new_N) {
    N = new_N;
    modelled_data.resize(N);
    classes.resize(N);
    for (int i = 0; i < c; i++) {
        c1_data_by_class[i];
        c2_data_by_class[i];
    }
}

double logic::g(int class_number, QVector<double>& x) {
    auto& sigma = sigma_matrices[class_number];
    auto& mu = mu_vectors[class_number];
    auto& sigma_inv = sigma_inv_matrices[class_number];
    double c1 = Multiply(Multiply(x, scalar_multiply(sigma_inv, -1/2)), x);
    double c2 = Multiply((Multiply(sigma_inv, mu)), x);
    return c1 + c2 + log(P[class_number]) - log(sqrt(det(sigma)))
            - 1/2 * Multiply(Multiply(mu, sigma_inv), mu);
}

double logic::d_function(int class_number, int vector_number) {
    QVector<double> x = modelled_data[vector_number];
    auto& sigma = sigma_matrices[class_number];
    auto& mu = mu_vectors[class_number];
    auto& sigma_inv = sigma_inv_matrices[class_number];
    double c1 = log(P[class_number]) - log(sqrt(det(sigma)));
    auto t = Subtract(x, mu);
    auto t_vect = Multiply(t, sigma_inv);
    double c2 = Multiply(t_vect, t);
    return c1 - 0.5*c2;
}


double logic::LikehoodFunction(int i, int class_number) {
    auto& mu = mu_vectors[class_number];
    auto& sigma = sigma_matrices[class_number];
    auto& sigma_inv = sigma_inv_matrices[class_number];
    QVector<double> &x = modelled_data[i];

    double result = pow((2*M_PI), (d/2)) * sqrt(det(sigma));

    QVector<double> t = Subtract(x, mu);
    double coef = Multiply(Multiply(t, sigma_inv), t);
    double e = exp(-0.5 * coef);
    result = 1 / result * e;
    return result;
}

logic::logic() {
    std::random_device rd{};
    std::mt19937 g{rd()};
    gen = g;
    std::normal_distribution<> d{0,1};
    distr = d;
}

void logic::ModelNormalDistributionVector(QVector<double>& result) {
    if (result.size() < d)
        result.resize(d);
    for (int i = 0; i < d; i++)
        result[i] = (distr(gen));
}

void logic::GetPointsData(int class_number, QVector<double>& x, QVector<double>& y, std::string type) {
    if (type == "Original")
        for (int i = 0; i < N; i++) {
            if (classes[i] == class_number) {
                x.append(modelled_data[i][axis1]);
                y.append(modelled_data[i][axis2]);
            }
        }
    if (type == "Risk")
        for (int i = 0; i < c2_data_by_class[class_number].size(); i++) {
            x.append(modelled_data[c2_data_by_class[class_number][i]][axis1]);
            y.append(modelled_data[c2_data_by_class[class_number][i]][axis2]);
        }
    if (type == "Error")
        for (int i = 0; i < c1_data_by_class[class_number].size(); i++) {
            x.append(modelled_data[c1_data_by_class[class_number][i]][axis1]);
            y.append(modelled_data[c1_data_by_class[class_number][i]][axis2]);
        }
}

void logic::CalculateRandomVector(int class_number, QVector<double>& result) {
    QVector<QVector<double>>& A = A_matrices[class_number];
    QVector<double>& mean_vector = mu_vectors[class_number];
    QVector<double> random_vector(d);
    ModelNormalDistributionVector(random_vector);
    result = Multiply(A, random_vector);
    result = Add(result, mean_vector);
}

void logic::ClassifyMinRisk() {
    sum_risk = 0;
    for (int i = 0; i < modelled_data.size(); i++) {
        double min = 100000;
        int predict_class_number = -1;
        int true_class = classes[i];        
        double cur_risk = 0;
        for (int class_number = 0; class_number < c; class_number++) {
            cur_risk =  risk(class_number, i);
            if (min > cur_risk) {
                min = cur_risk;
                predict_class_number = class_number;
            }
        }
        sum_risk += min;
        c2_data_by_class[predict_class_number].append(i);
        risk_transformation_matrix[predict_class_number][true_class] += 1;
    }
}

void logic::ClassifyMinError() {
    for (int i = 0; i < modelled_data.size(); i++) {
        double max = -100000;
        int predict_class_number = -1;
        int true_class = classes[i];
        for (int class_number = 0; class_number < c; class_number++) {
            double d_ = d_function(class_number, i);
            if (max < d_) {
                max = d_;
                predict_class_number = class_number;
            }
        }
        c1_data_by_class[predict_class_number].append(i);
        transformation_matrix[predict_class_number][true_class] += 1;
    }
}

void logic::ClearData() {
    classes.clear();
    modelled_data.clear();
    c1_data_by_class.clear();
    c2_data_by_class.clear();
    //    modelled_data.resize(N);
    transformation_matrix.clear();
//    transformation_matrix.resize(c);

    risk_transformation_matrix.clear();
//    risk_transformation_matrix.resize(c);
//    for (int i = 0; i < c; i++) {
//        c1_data_by_class[i];
//        c2_data_by_class[i];
//        transformation_matrix[i].resize(c);
//        risk_transformation_matrix[i].resize(c);
//    }
}

void logic::Model() {
    GenerateClasses(); // pick a class for each one of N vectors
    for (int i = 0; i < modelled_data.size(); i++) {
        int class_numder = classes[i]; // get class number (pregenerated)
        CalculateRandomVector(class_numder, modelled_data[i]);
    }
}

void logic::GenerateClasses() {
    classes.resize(N);
    std::vector<int> classes_count(c);
    for (unsigned i = 0; i < classes.size(); i++) {
        double probability = (double)rand() / (double)RAND_MAX;
        if (probability <= P[0]) {
            classes[i] = 0;
            classes_count[0]++;
            continue;
        }
        double sum = P[1];
        for (unsigned j = 1; j < P.size(); j++) {
            if (P[j - 1] < probability && probability <= sum + P[j]) {
                classes[i] = j;
                classes_count[j]++;
                break;
            }
            sum += P[j];
        }
    }
    for (int j = 0; j < c; j++) {
        std::cout << "class " << j << ": " << classes_count[j] << " ";
    }
    std::cout << std::endl;
}

void logic::CalculateCoefMatrix() {
//    A_matrices.resize(c);
    for (int class_number = 0; class_number < c; class_number++) {
        auto& A_matrix = A_matrices[class_number];
        auto& B_norm = sigma_matrices[class_number];

        auto Var = var_vectors[class_number];
        QVector<QVector<double>> B(d);
        for (unsigned int i = 0; i < B_norm.size(); i++) {
            B[i].resize(B_norm[i].size());
            for (unsigned int j = 0; j < B_norm[i].size(); j++) {
                B[i][j] = B_norm[i][j] * sqrt(Var[i]) * sqrt(Var[j]);
            }
        }
        for (int i = 0; i < d; i++) {
            for (int j = 0; j <= i; j++) {
                double sum1 = 0;
                for (int k = 0; k < j; k++) {
                    sum1 += A_matrix[i][k] * A_matrix[j][k];
                }
                double sum2 = 0;
                for (int k = 0; k < j; k++) {
                    sum2 += pow(A_matrix[j][k], 2);
                }
                A_matrix[i][j] = (B[i][j] - sum1) / sqrt(B[j][j] - sum2);
            }
        }
        sigma_matrices[class_number] = B;
    }
}

double logic::risk(int class_number, int vector_number) {
    double result = 0;
    for (int j = 0; j < c; j++) {
        double a = AposterioriProbability(j, vector_number);
        result += lambdas[class_number][j] * a;
    }
    return result;
}

void logic::CalculateSigmaInvMatrices() {
    for (int i = 0; i < sigma_matrices.size(); i++) {
        sigma_inv_matrices[i] = (InvMatrix(sigma_matrices[i]));
    }
}

double logic::AposterioriProbability(int class_number, int vector_number) {
    double sum = 0;
    for (int i = 0; i < c; i++) {
        double l = LikehoodFunction(vector_number, i);
        sum += P[i] * l;
    }
    double res = P[class_number] * LikehoodFunction(vector_number, class_number);
    return res / sum;
}

double logic::CalculateErrorProbability() {
    double sum = 0;
    for (int i = 0; i < transformation_matrix.size(); i++) {
        for (int j = 0; j < transformation_matrix.size(); j++) {
            sum += (i == j ? 0 : transformation_matrix[i][j]);
        }
    }
    return sum / N;
}

double logic::CalculateMeanRisk() {
    return sum_risk / N;
}
