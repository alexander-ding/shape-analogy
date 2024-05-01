#include "mainwindow.h"
#include "settings.h"

#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QRadioButton>
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

    QWidget *controlsWidget = new QWidget();
    QVBoxLayout *vLayout = new QVBoxLayout();
    vLayout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    controlsWidget->setLayout(vLayout);
    controlsWidget->setMaximumWidth(150);

    container->addWidget(controlsWidget);

    container->addWidget(editorWidget);
    container->addWidget(analogyWidget);

    this->setLayout(container);

    settings.loadSettingsOrDefaults();

    addHeading(vLayout, "Mode");
    addRadioButton(vLayout, "ARAP", settings.mode == ARAP, [this]{ setMode(ARAP); });
    addRadioButton(vLayout, "Push", settings.mode == PUSH, [this]{ setMode(PUSH); });
    addRadioButton(vLayout, "Hammer", settings.mode == HAMMER, [this]{ setMode(HAMMER); });
}

MainWindow::~MainWindow()
{
    delete editorWidget;
    delete analogyWidget;
}

void MainWindow::addHeading(QBoxLayout *layout, QString text) {
    QFont font;
    font.setPointSize(16);
    font.setBold(true);

    QLabel *label = new QLabel(text);
    label->setFont(font);
    layout->addWidget(label);
}

void MainWindow::addRadioButton(QBoxLayout *layout, QString text, bool value, auto function) {
    QRadioButton *button = new QRadioButton(text);
    button->setChecked(value);
    layout->addWidget(button);
    connect(button, &QRadioButton::clicked, this, function);
}

void MainWindow::setMode(int type) {
    settings.mode = type;
    // m_canvas->settingsChanged();
}


