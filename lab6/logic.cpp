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



QVector<QVector<double>> Multiply(QVector<double> lhs, QVector<double> rhs, bool Matrix) {

    if (!Matrix)
        return QVector<QVector<double>>();
    QVector<QVector<double>> result(lhs.size(), QVector<double>(rhs.size()));
    for (int i = 0; i < rhs.size(); i++)
        for (int j = 0; j < lhs.size(); j++)
            result[i][j] = lhs[i] * rhs[j];
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

QVector<double> scalar_multiply(QVector<double>& lhs, double rhs) {
    QVector<double> result(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        result[i] = lhs[i] * rhs;
    }
    return result;
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

QVector<QVector<double>> Add(QVector<QVector<double>>& lhs, QVector<QVector<double>>& rhs) {
    QVector<QVector<double>> result(lhs.size(), QVector<double>(rhs.size()));
    for (int i = 0; i < lhs.size(); i++) {
        for (int j = 0; j < lhs.size(); j++)
            result[i][j] = lhs[i][j] + rhs[i][j];
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
    //????????????????????????????????
    size_t n = A0.size();
    for (size_t row = 0; row < n - 1; row++) {
        for (size_t col = row + 1; col < n; col++)
            std::swap(A0[col][row], A0[row][col]);
    }
}

double det(QVector<QVector<double>> matrix) {
    int n = matrix.size();
    QVector<QVector<double>> B = matrix;
    //???????????????????? ?????????????? ?? ???????????????????????????????????? ????????
    for (int step = 0; step < n - 1; step++)
        for (int row = step + 1; row < n; row++) {
            double coeff = -B[row][step] / B[step][step]; //?????????? ????????????
            for(int col = step; col < n; col++)
                B[row][col] += B[step][col] * coeff;
        }
    //???????????????????? ???????????????????????? ?????? ???????????????????????? ?????????????????? ?????????????? ??????????????????
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
    //???????????????????? ?????????????????? ??????????????
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if(i == j) E[i][j] = 1.0;
            else E[i][j] = 0.0;
        }
    }
    //?????????????????????????? ???????????????? ?????????????? ?? ??????????????????????????????????
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

    //???????????????????????? ?????????????????? ?????????????? ???? ????????????????
    //?? ???????????????? ???? ??????????????????
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
    for (int i = 0; i < c; i++) {
        transformation_matrix[i].resize(c);
    }
    for (int i = 0; i < c; i++) {
        c1_data_by_class[i];
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

void logic::SetN(int new_N_train, int new_N_test) {
    N_train = new_N_train;
    N_test = new_N_test;
    train_data.resize(N_train);
    train_classes.resize(N_train);
    test_classes.resize(N_test);
    test_data.resize(N_test);
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
    QVector<double> x = train_data[vector_number];
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
    QVector<double> &x = train_data[i];

    double result = pow((2*M_PI), (d/2)) * sqrt(det(sigma));

    QVector<double> t = Subtract(x, mu);
    double coef = Multiply(Multiply(t, sigma_inv), t);
    double e = exp(-0.5 * coef);
    result = 1 / result * e;
//    std::cout << std::setprecision(10) << result << " " << -1/2 << " " << exp(-1/2 * coef) << " " << t[axis1] << "   " << t[axis2] << std::endl;
    return result;
}

logic::logic() {
    std::random_device rd{};
    std::mt19937 g{rd()};
    gen = g;
    std::normal_distribution<> d{0,1};
    distr = d;
    param_classifier = nullptr;
    bayes_classifier = nullptr;
    parzen_classifier = nullptr;
}

void logic::ModelNormalDistributionVector(QVector<double>& result) {
    if (result.size() < d)
        result.resize(d);
    for (int i = 0; i < d; i++)
        result[i] = (distr(gen));
}

void logic::GetPointsData(int class_number, QVector<double>& x, QVector<double>& y, std::string type) {
    if (type == "Original")
        for (int i = 0; i < N_test; i++) {
            if (test_classes[i] == class_number) {
                x.append(test_data[i][axis1]);
                y.append(test_data[i][axis2]);
            }
        }
    if (type == "Bayes") {
        auto& data_by_class = bayes_classifier->test_data_by_class;
        for (int i = 0; i < data_by_class[class_number].size(); i++) {
            x.append(test_data[data_by_class[class_number][i]][axis1]);
            y.append(test_data[data_by_class[class_number][i]][axis2]);
        }
    }
    if (type == "Parametric") {
        auto& data_by_class = param_classifier->test_data_by_class;
        for (int i = 0; i < data_by_class[class_number].size(); i++) {
            x.append(test_data[data_by_class[class_number][i]][axis1]);
            y.append(test_data[data_by_class[class_number][i]][axis2]);
        }
    }
    if (type == "Parzen") {
        auto& data_by_class =
                parzen_classifier->test_data_by_class;
        for (int i = 0; i < data_by_class[class_number].size(); i++) {
            x.append(test_data[data_by_class[class_number][i]][axis1]);
            y.append(test_data[data_by_class[class_number][i]][axis2]);
        }
//        parzen_classifier
    }
    if (type == "kNN") {
        auto& data_by_class =
                knn_classifier->test_data_by_class;
        for (int i = 0; i < data_by_class[class_number].size(); i++) {
            x.append(test_data[data_by_class[class_number][i]][axis1]);
            y.append(test_data[data_by_class[class_number][i]][axis2]);
        }
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

void logic::ClassifyMinError() {
    for (int i = 0; i < train_data.size(); i++) {
        double max = -100000;
        int predict_class_number = -1;
        int true_class = test_classes[i];
        for (int class_number = 0; class_number < c; class_number++) {
//            double d_ = g(class_number, modelled_data[i]);
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
    test_classes.clear();
    train_data.clear();
    test_data.clear();
    train_classes.clear();
    c1_data_by_class.clear();
    //    modelled_data.resize(N);
    transformation_matrix.clear();
//    transformation_matrix.resize(c);

//    risk_transformation_matrix.clear();
//    risk_transformation_matrix.resize(c);
//    for (int i = 0; i < c; i++) {
//        c1_data_by_class[i];
//        c2_data_by_class[i];
//        transformation_matrix[i].resize(c);
//        risk_transformation_matrix[i].resize(c);
//    }
}

void logic::Model() {
    GenerateClasses(); // pick class for each one of N vectors
    for (int i = 0; i < test_data.size(); i++) {
        int class_numder = test_classes[i]; // get class number (pregenerated)
        CalculateRandomVector(class_numder, test_data[i]);
    }
    for (int i = 0; i < train_data.size(); i++) {
        int class_numder = train_classes[i]; // get class number (pregenerated)
        CalculateRandomVector(class_numder, train_data[i]);
    }
}

void logic::GenerateClasses() {
    std::vector<int> classes_count(c);
    for (unsigned i = 0; i < test_classes.size(); i++) {
        double probability = (double)rand() / (double)RAND_MAX;
        if (probability <= P[0]) {
            test_classes[i] = 0;
            classes_count[0]++;
            continue;
        }
        double sum = P[1];
        for (unsigned j = 1; j < P.size(); j++) {
            if (P[j - 1] < probability && probability <= sum + P[j]) {
                test_classes[i] = j;
                classes_count[j]++;
                break;
            }
            sum += P[j];
        }
    }
    classes_count.clear();
    classes_count.resize(c);
    for (unsigned i = 0; i < train_classes.size(); i++) {
        double probability = (double)rand() / (double)RAND_MAX;
        if (probability <= P[0]) {
            train_classes[i] = 0;
            classes_count[0]++;
            continue;
        }
        double sum = P[1];
        for (unsigned j = 1; j < P.size(); j++) {
            if (P[j - 1] < probability && probability <= sum + P[j]) {
                train_classes[i] = j;
                classes_count[j]++;
                break;
            }
            sum += P[j];
        }
    }
//    for (int j = 0; j < c; j++) {
//        std::cout << "class " << j << ": " << classes_count[j] << " ";
//    }
//    std::cout << std::endl;
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
//    std::cout << sum << std::endl;
    return sum / 1; // N;
}


void ParamClassifier::Train() {
    classes_count.resize(c);
    for (int i = 0; i < N_train; i++) {
        int class_number = train_classes[i];
        auto& mu = mu_vectors[class_number];
        auto& sigma = sigma_matrices[class_number];
        for (int j = 0; j < d; j++) {
            mu[j] += train_data[i][j];
        }
        classes_count[class_number]++;
    }
    for (int i = 0; i < c; i++) {
        P[i] /= classes_count[i] / N_train;
    }
    for (int class_number = 0; class_number < c; class_number++) {
        for (int i = 0; i < d; i++) {
            mu_vectors[class_number][i] /= classes_count[class_number];
        }
    }
    for (int i = 0; i < N_train; i++) {
        int class_number = train_classes[i];
        auto& mu = mu_vectors[class_number];
        auto& sigma = sigma_matrices[class_number];
        auto t = Multiply(Subtract(train_data[i], mu), Subtract(train_data[i], mu), true);
        sigma = Add(sigma, t);
    }
    for (int class_number = 0; class_number < c; class_number++) {
        sigma_matrices[class_number] = scalar_multiply(sigma_matrices[class_number], (double)(1 / classes_count[class_number]));
    }


    for (int i = 0; i < sigma_matrices.size(); i++) {
        sigma_inv_matrices[i] = (InvMatrix(sigma_matrices[i]));
    }

}


void ParamClassifier::CalculateCoefMatrix() {
//    A_matrices.resize(c);
    for (int class_number = 0; class_number < c; class_number++) {
        auto& A_matrix = A_matrices[class_number];
//        auto& B_norm = sigma_matrices[class_number];
        auto& B = sigma_matrices[class_number];
//        auto Var = var_vectors[class_number];
//        QVector<QVector<double>> B(d);
//        for (unsigned int i = 0; i < B_norm.size(); i++) {
//            B[i].resize(B_norm[i].size());
//            for (unsigned int j = 0; j < B_norm[i].size(); j++) {
//                B[i][j] = B_norm[i][j] * sqrt(Var[i]) * sqrt(Var[j]);
//            }
//        }
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
//        sigma_matrices[class_number] = B;
    }
}


void ParamClassifier::CalculateSigmaInvMatrices() {
    for (int i = 0; i < sigma_matrices.size(); i++) {
        sigma_inv_matrices[i] = (InvMatrix(sigma_matrices[i]));
    }
}

double ParamClassifier::d_function(int class_number, int vector_number) {
    QVector<double> x = test_data[vector_number];
    auto& sigma = sigma_matrices[class_number];
    auto& mu = mu_vectors[class_number];
    auto& sigma_inv = sigma_inv_matrices[class_number];
    double c1 = log(P[class_number]) - log(sqrt(det(sigma)));
    auto t = Subtract(x, mu);
    auto t_vect = Multiply(t, sigma_inv);
    double c2 = Multiply(t_vect, t);
    return c1 - 0.5*c2;
}


void ParamClassifier::Classify() {
    error = 0;
    for (int i = 0; i < test_data.size(); i++) {
        double max = -100000;
        int predict_class_number = -1;
        int true_class = test_classes[i];
        for (int class_number = 0; class_number < c; class_number++) {
            double d_ = d_function(class_number, i);
            if (max < d_) {
                max = d_;
                predict_class_number = class_number;
            }
        }
        test_data_by_class[predict_class_number].append(i);
        transformation_matrix[predict_class_number][true_class] += 1;
    }
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++)
            error += i == j ? 0 : transformation_matrix[i][j];
    }
    error /= N_test;
}


ParamClassifier::ParamClassifier(
        int c, int d,
        int N_train, int N_test,
        QVector<double>& P,
        QVector<double>& train_classes,
        QVector<double>& test_classes,
        QVector<QVector<double>>& train_data,
        QVector<QVector<double>>& test_data)
    : P(P),
      train_classes(train_classes),
      test_classes(test_classes),
      train_data(train_data),
      test_data(test_data),
//      c(c), d(d),
      N_train(N_train),
      N_test(N_test) { SetC(c); SetD(d); }



void ParamClassifier::SetC(int new_c) {
    c = new_c;
    error = 0;
    sigma_matrices.resize(c);
    sigma_inv_matrices.resize(c);
    A_matrices.resize(c);

    mu_vectors.resize(c);
//    var_vectors.resize(c);
//    test_data_by_class;

//    P.resize(c);
    transformation_matrix.resize(c);
    for (int i = 0; i < c; i++) {
        transformation_matrix[i].resize(c);
        test_data_by_class[i];
    }
}

void ParamClassifier::SetD(int new_d) {
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
//        var_vectors[class_number].resize(d);
    }
}

void ParamClassifier::ClearData() {
    sigma_matrices.clear();
    sigma_inv_matrices.clear();
    A_matrices.clear();
    mu_vectors.clear();
//    var_vectors.clear();
    transformation_matrix.clear();
    test_data_by_class.clear();
    classes_count.clear();
}



void MinErrorClassifier::Train() {

}


void MinErrorClassifier::CalculateCoefMatrix() {
//    A_matrices.resize(c);
    for (int class_number = 0; class_number < c; class_number++) {
        auto& A_matrix = A_matrices[class_number];
//        auto& B_norm = sigma_matrices[class_number];
        auto& B = sigma_matrices[class_number];
//        auto Var = var_vectors[class_number];
//        QVector<QVector<double>> B(d);
//        for (unsigned int i = 0; i < B_norm.size(); i++) {
//            B[i].resize(B_norm[i].size());
//            for (unsigned int j = 0; j < B_norm[i].size(); j++) {
//                B[i][j] = B_norm[i][j] * sqrt(Var[i]) * sqrt(Var[j]);
//            }
//        }
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
//        sigma_matrices[class_number] = B;
    }
}


void MinErrorClassifier::CalculateSigmaInvMatrices() {
    for (int i = 0; i < sigma_matrices.size(); i++) {
        sigma_inv_matrices[i] = (InvMatrix(sigma_matrices[i]));
    }
}

double MinErrorClassifier::d_function(int class_number, int vector_number) {
    QVector<double> x = test_data[vector_number];
    auto& sigma = sigma_matrices[class_number];
    auto& mu = mu_vectors[class_number];
    auto& sigma_inv = sigma_inv_matrices[class_number];
    double c1 = log(P[class_number]) - log(sqrt(det(sigma)));
    auto t = Subtract(x, mu);
    auto t_vect = Multiply(t, sigma_inv);
    double c2 = Multiply(t_vect, t);
    return c1 - 0.5*c2;
}


void MinErrorClassifier::Classify() {
    error = 0;
    for (int i = 0; i < test_data.size(); i++) {
        double max = -100000;
        int predict_class_number = -1;
        int true_class = test_classes[i];
        for (int class_number = 0; class_number < c; class_number++) {
            double d_ = d_function(class_number, i);
            if (max < d_) {
                max = d_;
                predict_class_number = class_number;
            }
        }
        test_data_by_class[predict_class_number].append(i);
        transformation_matrix[predict_class_number][true_class] += 1;
    }
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++)
            error += i == j ? 0 : transformation_matrix[i][j];
    }
    error /= N_test;
}

void MinErrorClassifier::SetC(int new_c) {
    c = new_c;
    sigma_matrices.resize(c);
    sigma_inv_matrices.resize(c);
    A_matrices.resize(c);

    mu_vectors.resize(c);
//    var_vectors.resize(c);
//    test_data_by_class;

//    P.resize(c);
    transformation_matrix.resize(c);
    for (int i = 0; i < c; i++) {
        transformation_matrix[i].resize(c);
        test_data_by_class[i];
    }
}

void MinErrorClassifier::SetD(int new_d) {
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
//        var_vectors[class_number].resize(d);
    }
}

void MinErrorClassifier::ClearData() {
    sigma_matrices.clear();
    sigma_inv_matrices.clear();
    A_matrices.clear();
    mu_vectors.clear();
//    var_vectors.clear();
    transformation_matrix.clear();
    test_data_by_class.clear();
    classes_count.clear();
}


void logic::CreateClassifiers() {
    if (param_classifier)
        delete param_classifier;
    param_classifier = new ParamClassifier;

    param_classifier->P = (P);
    param_classifier->train_classes = train_classes;
    param_classifier->test_classes = test_classes;
    param_classifier->train_data = (train_data);
    param_classifier->test_data = (test_data);
    param_classifier->N_train = N_train;
    param_classifier->N_test = (N_test);
    param_classifier->SetC(c); param_classifier->SetD(d);

    if (bayes_classifier)
        delete bayes_classifier;
    bayes_classifier = new MinErrorClassifier;

    bayes_classifier->P = (P);
    bayes_classifier->train_classes = train_classes;
    bayes_classifier->test_classes = test_classes;
    bayes_classifier->train_data = (train_data);
    bayes_classifier->test_data = (test_data);
    bayes_classifier->N_train = N_train;
    bayes_classifier->N_test = (N_test);
    bayes_classifier->SetC(c); bayes_classifier->SetD(d);
    bayes_classifier->sigma_matrices = sigma_matrices;
    bayes_classifier->sigma_inv_matrices = sigma_inv_matrices;
    bayes_classifier->mu_vectors = mu_vectors;

    if (parzen_classifier)
        delete parzen_classifier;
    parzen_classifier = new ParzenClassifier;
    parzen_classifier->train_classes = train_classes;
    parzen_classifier->test_classes = test_classes;
    parzen_classifier->train_data = (train_data);
    parzen_classifier->test_data = (test_data);
    parzen_classifier->N_train = N_train;
    parzen_classifier->N_test = (N_test);
    parzen_classifier->SetC(c);
    parzen_classifier->SetD(d);
    parzen_classifier->P = P;

    if (knn_classifier)
        delete knn_classifier;
    knn_classifier = new NeighbourClassifier;
    knn_classifier->train_classes = train_classes;
    knn_classifier->test_classes = test_classes;
    knn_classifier->train_data = (train_data);
    knn_classifier->test_data = (test_data);
    knn_classifier->N_train = N_train;
    knn_classifier->N_test = (N_test);
    knn_classifier->SetC(c);
    knn_classifier->SetD(d);
    knn_classifier->P = P;

}



void ParzenClassifier::Classify() {
    for (int i = 0; i < N_train; i++) {
        int class_number = train_classes[i];
        classes_count[class_number]++;
    }
    double v1 = 2;
    for (int i = 0; i < c; i++) {
        v[i] = v1 / (sqrt(classes_count[i]));
    }
//    for (int i = 0; i < c; i++) {
//        P[i] = N_train / classes_count[i];
//    }
    error = 0;
    for (int i = 0; i < test_data.size(); i++) {
        double max = -100000;
        int predict_class_number = -1;
        int true_class = test_classes[i];
        for (int class_number = 0; class_number < c; class_number++) {
            double d_ = d_function(class_number, i);
            if (max < d_) {
                max = d_;
                predict_class_number = class_number;
//                std::cout << " " << test_data[i][0] << " " << d_ ;
            }
        }
        std::cout << "!" << i << std::endl;
        test_data_by_class[predict_class_number].append(i);
        transformation_matrix[predict_class_number][true_class] += 1;
    }
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++)
            error += i == j ? 0 : transformation_matrix[i][j];
    }
    error /= N_test;
}

double ParzenClassifier::Parzen(int vector_number, int class_number) {
    double result = 0;
    for (int i = 0; i < train_data.size(); i++) {
        if (class_number != train_classes[i])
            continue;
        auto t = Subtract(test_data[vector_number], train_data[i]);
        t = scalar_multiply(t, 1 / pow(v[class_number], (double)(d)));
        double c1;
        double sum = 0;
        for (int j = 0; j < t.size(); j++)
            sum += t[j] * t[j];

        c1 = exp(-0.5 * sum) * 1 / pow((2 * M_PI), (double)(d / 2));;

        c1 /= v[class_number];
        result += c1;
    }
    return result / classes_count[class_number];
}

double ParzenClassifier::u(QVector<double>& u_) {
    double sum = 0, result = 0;
    for (int j = 0; j < u_.size(); j++)
        sum += u_[j] * u_[j];

    result = exp(-0.5 * sum) * 1 / pow((2 * M_PI), (double)(d / 2));;
}

void ParzenClassifier::SetD(int new_d) {
    d = new_d;
}

void ParzenClassifier::SetC(int new_c) {
    c = new_c;
    v.resize(c);
    P.resize(c);
    classes_count.resize(c);
    transformation_matrix.resize(c);
    for (int i = 0; i < c; i++) {
        transformation_matrix[i].resize(c);
        test_data_by_class[i];
    }
}
void ParzenClassifier::ClearData() {
    P.clear();
    transformation_matrix.clear();
    test_data_by_class.clear();
    classes_count.clear();
    v.clear();
}
double ParzenClassifier::d_function(int class_number, int vector_number) {
//    std::cout << " " << Parzen(vector_number, class_number) << std::endl;
    return P[class_number] * Parzen(vector_number, class_number); // / classes_count[class_number];
}




void NeighbourClassifier::Classify() {
    for (int i = 0; i < N_train; i++) {
        int class_number = train_classes[i];
        classes_count[class_number]++;
    }
//    for (int i = 0; i < c; i++) {
//        P[i] = N_train / classes_count[i];
//    }
    error = 0;
    for (int i = 0; i < test_data.size(); i++) {
        double max = -100000;
        int predict_class_number = -1;
        int true_class = test_classes[i];
        for (int class_number = 0; class_number < c; class_number++) {
            double d_ = d_function(class_number, i);
            if (max < d_) {
                max = d_;
                predict_class_number = class_number;
//                std::cout << " " << test_data[i][0] << " " << d_ ;
            }
        }
        std::cout << "!!" << i << std::endl;
        test_data_by_class[predict_class_number].append(i);
        transformation_matrix[predict_class_number][true_class] += 1;
    }
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++)
            error += i == j ? 0 : transformation_matrix[i][j];
    }
    error /= N_test;
}

//double ParzenClassifier::Parzen(int vector_number, int class_number) {
//    double result = 0;
//    for (int i = 0; i < train_data.size(); i++) {
//        if (class_number != train_classes[i])
//            continue;
//        auto t = Subtract(test_data[vector_number], train_data[i]);
//        t = scalar_multiply(t, 1 / pow(v[class_number], (double)(d)));
//        double c1;
//        double sum = 0;
//        for (int j = 0; j < t.size(); j++)
//            sum += t[j] * t[j];

//        c1 = exp(-0.5 * sum) * 1 / pow((2 * M_PI), (double)(d / 2));;

//        c1 /= v[class_number];
//        result += c1;
//    }
//    return result / classes_count[class_number];
//}


void NeighbourClassifier::SetD(int new_d) {
    d = new_d;
}

void NeighbourClassifier::SetC(int new_c) {
    c = new_c;
//    v.resize(c);
    P.resize(c);
    classes_count.resize(c);
    transformation_matrix.resize(c);
    for (int i = 0; i < c; i++) {
        transformation_matrix[i].resize(c);
        test_data_by_class[i];
    }
}
void NeighbourClassifier::ClearData() {
    P.clear();
    transformation_matrix.clear();
    test_data_by_class.clear();
    classes_count.clear();
//    v.clear();
}

double NeighbourClassifier::f(int vector_number,int class_number) {
    double k = 0, v = 0;
    while (k < sqrt(classes_count[class_number])) {
        k = 0;
        v += 0.5;
        for (int j = 0; j < train_data.size(); j++){
            if (train_classes[j] != class_number)
                continue;
            auto train_x = train_data[j];
            for (int i = 0; i < d; i++) {
                if ((test_data[vector_number][i] + v > train_x[i]) && test_data[vector_number][i] - v < train_x[i])
                    ;
                else
                    break;
                if (i == d - 1)
                    k++;
            }
        }
        if (k > sqrt(classes_count[class_number]) || k > N_train)
            break;
    }
    return k / (v * sqrt(classes_count[class_number]));
}

double NeighbourClassifier::d_function(int class_number, int vector_number) {
    return P[class_number] * f(vector_number, class_number);
}







