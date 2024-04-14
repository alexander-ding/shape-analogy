#pragma once

#include <QMainWindow>
#include "glwidget.h"
#include "analogy.h"

class MainWindow : public QWidget
{
    Q_OBJECT

public:
    MainWindow(Analogy analogy);
    ~MainWindow();

private:

    GLWidget *glWidget;
};
