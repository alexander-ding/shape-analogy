#include <QCoreApplication>
#include <QCommandLineParser>
#include <QtCore>
#include <QString>
#include <iostream>

#include "src/mesh.h"

using namespace std;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    QCommandLineParser parser;
    parser.addOptions({
        {"a", "Path to a' .obj file", "path"},
        {"b", "Path to b .obj file", "path"},
        {{"o", "out"}, "Path to output the generated b' .obj file", "path"},
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

    QString aPath = parser.value("a");
    QString bPath = parser.value("b");
    QString oPath = parser.value("o");

    Mesh aPrimeMesh(aPath);
    Mesh bMesh(bPath);
    Mesh bPrimeMesh(bMesh);
    // TODO: run algorithm

    bPrimeMesh.saveToFile(oPath);

    a.exit();
}
