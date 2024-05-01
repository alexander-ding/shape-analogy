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
    void loadSettingsOrDefaults();
};

extern Settings settings;

#endif // SETTINGS_H
