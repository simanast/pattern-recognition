#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <fstream>
#include <iostream>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow) {
        ui->setupUi(this);

        brush_colors.append((Qt::red));
        brush_colors.append((Qt::blue));
        brush_colors.append((Qt::green));
        brush_colors.append((Qt::magenta));
        brush_colors.append((Qt::black));
        brush_colors.append((Qt::yellow));
        brush_colors.append((Qt::cyan));
        Logic = new logic();
        ui->ClassifierType->addItem("Bayes");
        ui->ClassifierType->addItem("Parametric");
        ui->ClassifierType->addItem("Parzen");
        ui->ClassifierType->addItem("kNN");
    }

MainWindow::~MainWindow() {
    delete ui;
    delete Logic;
}

void MainWindow::ReadAxes() {
    Logic->axis1 = ui->Axis1->text().toInt();
    Logic->axis2 = ui->Axis2->text().toInt();
}

void MainWindow::ReadN() {
    Logic->SetN(ui->NTrainLine->text().toInt(), ui->NTestLine->text().toInt());
}

void MainWindow::ReadConfigFile() {
    std::ifstream file("test.txt");
    std::string section;
    int class_number;
    double value;
    int l;
    Logic->ClearData();
    file >> section;
    if (section == "c:")
        file >> l;
    Logic->SetC(l);

    file >> section;
    if (section == "d:")
        file >> l;
    Logic->SetD(l);

    file >> section;
    if (section == "P:") {
        for (int i = 0; i < Logic->c; i++) {
            file >> Logic->P[i];
        }
    }

    while (file >> class_number) {
        file >> section;
        if (section == "Sigma:") {
            for (int i = 0; i < Logic->d; i++) {
                for (int j = 0; j < Logic->d; j++) {
                    file >> value;
                    Logic->sigma_matrices[class_number][i][j] = value;
                }
            }
        }
        file >> section;
        if (section == "var:") {
            for (int i = 0; i < Logic->d; i++) {
                file >> value;
                Logic->var_vectors[class_number][i] = value;
            }
        }
        file >> section;
        if (section == "mu:") {
            for (int i = 0; i < Logic->d; i++) {
                file >> value;
                Logic->mu_vectors[class_number][i] = value;
            }
        }

        continue;
    }
    file.close();
}

void MainWindow::PlotPoints(std::string type) {
    double x_min, y_min, x_max, y_max;
    QCustomPlot* plot;
    if (type == "Original")
        plot = ui->OriginalPoints;
    else
        plot = ui->PointsError;
    plot->clearGraphs();
    QVector<double> x, y;
    for (int i = 0; i < Logic->c; i++) {
        points = new QCPGraph(plot->xAxis, plot->yAxis);
        points->setLineStyle(QCPGraph::lsNone);
        points->setScatterStyle(QCPScatterStyle::ssCircle);
        points->setPen(QPen(brush_colors[i], 1));
        Logic->GetPointsData(i, x, y, type);
        points->addData(x,y);
        x_min = std::min(x_min, *std::min_element(x.begin(), x.end()));
        y_min = std::min(y_min, *std::min_element(y.begin(), y.end()));
        x_max = std::max(x_max, *std::max_element(x.begin(), x.end()));
        y_max = std::max(y_max, *std::max_element(y.begin(), y.end()));
        x.clear(); y.clear();
        if (x_max - x_min > y_max - y_min)
            y_max = y_min + (x_max - x_min);
        else
            x_max = x_min + (y_max - y_min);
    }
    plot->xAxis->setRange(QCPRange(x_min - 1, x_max + 1));
    plot->yAxis->setRange(QCPRange(y_min - 1, y_max + 1));
    plot->replot();
}

void ShowTable(QTableWidget* table, QVector<QVector<double>>& data) {
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            table->setItem(i, j, new QTableWidgetItem(QString::number(data[i][j])));
        }
    }
}

void MainWindow::ShowTransformationMatrix() {
    auto table = ui->TransformationMatrix;
    auto& data = Logic->transformation_matrix;
    ShowTable(table, data);
}

void MainWindow::on_ModelButton_clicked() {
    ReadConfigFile();
    ReadN();
    ReadAxes();
    CreateTables();
    Logic->CalculateSigmaInvMatrices();
    Logic->CalculateCoefMatrix();
    Logic->Model();

//    Logic->ClassifyMinError();

    PlotPoints("Original");
//    ShowTransformationMatrix();
//    ShowMetrics();
}

void MainWindow::ShowMetrics(double error) {
    ui->PError->setText(QString::number(error));
}

void MainWindow::CreateTables() {
    ui->TransformationMatrix->setColumnCount(Logic->c);
    ui->TransformationMatrix->setRowCount(Logic->c);
}

void MainWindow::on_TestButton_clicked() {
    Logic->CreateClassifiers();
    auto type = ui->ClassifierType->currentText();
    if (type == "Parametric") {
        Logic->param_classifier->Train();
        Logic->param_classifier->Classify();
        PlotPoints("Parametric");
        ShowTable(ui->TransformationMatrix, Logic->param_classifier->transformation_matrix);
        ShowMetrics(Logic->param_classifier->error);
//        for (int i = 0; i < Logic->param_classifier->sigma_matrices[1].size(); i++) {
//            for (int j = 0; j < Logic->param_classifier->sigma_matrices[1][i].size(); j++) {
//                std::cout << Logic->param_classifier->sigma_matrices[1][i][j] << " ";
//            }
//            std::cout << std::endl;
//        }
        Logic->param_classifier->ClearData();
    }
    if (type == "Bayes") {
        Logic->bayes_classifier->Classify();
        PlotPoints("Bayes");
        ShowTable(ui->TransformationMatrix, Logic->bayes_classifier->transformation_matrix);
        ShowMetrics(Logic->bayes_classifier->error);
        std::cout << std::endl;

        Logic->bayes_classifier->ClearData();
    }
    if (type == "Parzen") {
        Logic->parzen_classifier->Classify();
        PlotPoints("Parzen");
        ShowTable(ui->TransformationMatrix, Logic->parzen_classifier->transformation_matrix);
        ShowMetrics(Logic->parzen_classifier->error);
        Logic->parzen_classifier->ClearData();
    }
    if (type == "kNN") {
        Logic->knn_classifier->Classify();
        PlotPoints("kNN");
        ShowTable(ui->TransformationMatrix, Logic->knn_classifier->transformation_matrix);
        ShowMetrics(Logic->knn_classifier->error);
        Logic->knn_classifier->ClearData();
    }
    for (int i = 0; i < Logic->sigma_matrices[1].size(); i++) {
        for (int j = 0; j < Logic->sigma_matrices[1][i].size(); j++) {
            std::cout << Logic->sigma_matrices[1][i][j] << " ";
        }
        std::cout << std::endl;
    }
}
