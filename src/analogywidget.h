#pragma once

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include "graphics/camera.h"
#include "graphics/shader.h"

#include <QOpenGLWidget>
#include <QElapsedTimer>
#include <QTimer>
#include <memory>
#include "analogy.h"
#include "graphics/shape.h"

class AnalogyWidget : public QOpenGLWidget
{
    Q_OBJECT

public:
    AnalogyWidget(Analogy analogy, QWidget *parent = nullptr);
    ~AnalogyWidget();

    Analogy& getAnalogy() { return this->m_analogy; }
    void syncAnalogy();

private:
    static const int FRAMES_TO_AVERAGE = 60;

    // private:
    // Basic OpenGL Overrides
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int w, int h) override;

    // Event Listeners
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;

private:
    const double m_timestep = 0.001;
    QElapsedTimer m_deltaTimeProvider; // For measuring elapsed time
    QTimer m_intervalTimer;            // For triggering timed events

    Camera m_camera;
    Shape m_b;
    Shape m_bPrime;
    Shader *m_shader;
    Analogy m_analogy;

    float m_forward;
    float m_sideways;
    float m_vertical;

    int m_lastX;
    int m_lastY;

    bool m_capture;

private slots:

    // Physics Tick
    void tick();
};