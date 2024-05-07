#include "settings.h"

Settings settings;

void Settings::loadSettingsOrDefaults() {
    mode = ARAP;
    radius = 50.0;
    mirrorX = false;
    mirrorY = false;
    mirrorZ = false;
}
