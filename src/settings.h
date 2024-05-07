#ifndef SETTINGS_H
#define SETTINGS_H

#include <QObject>

enum Mode {
    ARAP,
    PUSH,
    HAMMER
};

struct Settings {
    int mode;
    float radius;
    bool mirrorX;
    bool mirrorY;
    bool mirrorZ;
    void loadSettingsOrDefaults();
};

extern Settings settings;

#endif // SETTINGS_H
