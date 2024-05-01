#include "settings.h"

Settings settings;

void Settings::loadSettingsOrDefaults() {
    mode = ARAP;
    radius = 50.0;
}
