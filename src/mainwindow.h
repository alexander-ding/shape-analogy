#pragma once

#include <QMainWindow>
#include "analogywidget.h"
#include "editorwidget.h"
#include "analogy.h"
#include "arap.h"
#include <QBoxLayout>

class MainWindow : public QWidget
{
    Q_OBJECT

public:
    MainWindow(Analogy analogy);
    ~MainWindow();

private:
    EditorWidget *editorWidget;
    AnalogyWidget *analogyWidget;

    void addHeading(QBoxLayout *layout, QString text);
    void addRadioButton(QBoxLayout *layout, QString text, bool value, auto function);

    void setMode(int type);
};
