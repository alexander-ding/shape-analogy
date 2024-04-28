#pragma once

#include <QMainWindow>
#include "analogywidget.h"
#include "editorwidget.h"
#include "analogy.h"
#include "arap.h"

class MainWindow : public QWidget
{
    Q_OBJECT

public:
    MainWindow(Analogy analogy);
    ~MainWindow();

private:
    EditorWidget *editorWidget;
    AnalogyWidget *analogyWidget;
};
