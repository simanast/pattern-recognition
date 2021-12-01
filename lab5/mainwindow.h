#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qcustomplot.h>
#include <logic.h>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
        Q_OBJECT

    public:
        MainWindow(QWidget *parent = nullptr);
        ~MainWindow();

        void ReadConfigFile();
        void ReadAxes();
        void ReadN();
        void PlotPoints();

        void PlotPointsClassified();
        void ShowTransformationMatrix();
        void ShowRiskTransformationMatrix();
        void CreateTables();
        void ShowMetrics();
        QVector<QBrush> brush_colors;

    private slots:
        void on_ModelButton_clicked();

    private:
        Ui::MainWindow *ui;
        logic* Logic;
        QCPCurve* contours;
        QCPGraph* points;
};
#endif // MAINWINDOW_H
