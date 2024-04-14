#include "mainwindow.h"
#include <QHBoxLayout>
#include "analogy.h"

MainWindow::MainWindow(Analogy analogy)
{
    glWidget = new GLWidget(analogy);

    QHBoxLayout *container = new QHBoxLayout;
    container->addWidget(glWidget);
    this->setLayout(container);
}

MainWindow::~MainWindow()
{
    delete glWidget;
}
