#include <QCommandLineParser>
#include <QtCore>
#include <QString>
#include <QApplication>
#include <QSurfaceFormat>
#include <QScreen>
#include <iostream>

#include "mesh.h"
#include "mainwindow.h"

using namespace std;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QCommandLineParser parser;
    parser.addOptions({
        {"a", "Path to a' .obj file", "path"},
        {"b", "Path to b .obj file", "path"},
        {{"o", "out"}, "Path to output the generated b' .obj file", "path"},
        {{"l", "lambda"}, "Value for lambda (float)", "lambda", "10"},
        {"s", "Sphere deformation", "true"}
    });
    parser.addHelpOption();
    parser.process(a);

    if (!parser.isSet("a")) {
        cerr << "flag -a not set" << std::endl;
        return -1;
    }
    if (!parser.isSet("b")) {
        cerr << "flag -b not set" << std::endl;
        return -1;
    }
    if (!parser.isSet("o")) {
        cerr << "flag -o not set" << std::endl;
        return -1;
    }
    if (!parser.isSet("s")) {
        cerr << "flag -s not set" << std::endl;
        return -1;
    }

    QString aPath = parser.value("a");
    QString bPath = parser.value("b");
    QString oPath = parser.value("o");
    float lambdaVal = parser.value("l").toFloat();
    bool sphereDeform = (parser.value("s")).toLower() == "true";

    Mesh aPrimeMesh(aPath);
    Mesh bMesh(bPath);
    Analogy analogy(aPrimeMesh, bMesh, sphereDeform);
    analogy.computeBPrime(lambdaVal);
    analogy.getBPrime().saveToFile(oPath);

    QApplication::setApplicationName("Shape Shifters");
    QApplication::setOrganizationName("CS 2240");
    QApplication::setApplicationVersion(QT_VERSION_STR);

    // Set OpenGL version to 4.1 and context to Core
    QSurfaceFormat fmt;
    fmt.setVersion(4, 1);
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(fmt);


    // Create a GUI window
    MainWindow w(analogy, oPath, aPath);
    w.resize(1000, 500);
    int desktopArea = QGuiApplication::primaryScreen()->size().width() *
                      QGuiApplication::primaryScreen()->size().height();
    int widgetArea = w.width() * w.height();
    if (((float)widgetArea / (float)desktopArea) < 0.75f)
        w.show();
    else
        w.showMaximized();


    return a.exec();
}
