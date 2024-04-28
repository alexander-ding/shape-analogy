#include "mainwindow.h"
#include <QHBoxLayout>
#include "analogy.h"
#include "mesh.h"
#include "editorwidget.h"

MainWindow::MainWindow(Analogy analogy)
{
    Mesh mesh("meshes/sphere.obj");
    editorWidget = new EditorWidget(mesh);
    analogyWidget = new AnalogyWidget(analogy);
    editorWidget->onUpdate([&](EditorWidget* e) {
        Mesh m = e->getARAP().getMesh();
        analogyWidget->getAnalogy().setAPrime(m);
        analogyWidget->syncAnalogy();
    });

    QHBoxLayout *container = new QHBoxLayout;
    container->addWidget(editorWidget);
    container->addWidget(analogyWidget);
    this->setLayout(container);
}

MainWindow::~MainWindow()
{
    delete editorWidget;
    delete analogyWidget;
}
