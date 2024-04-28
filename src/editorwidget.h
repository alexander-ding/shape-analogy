#pragma once

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include "graphics/camera.h"
#include "graphics/shader.h"
#include "graphics/shape.h"
#include "arap.h"

#include <QOpenGLWidget>
#include <QElapsedTimer>
#include <QTimer>
#include <functional>

class EditorWidget : public QOpenGLWidget
{
    Q_OBJECT

public:
    EditorWidget(Mesh& mesh, QWidget *parent = nullptr);
    ~EditorWidget();

    void onUpdate(std::function<void(EditorWidget*)> f);
    ARAP& getARAP() { return this->m_arap; }

private:
    static const int FRAMES_TO_AVERAGE = 30;

    Eigen::Vector3f transformToWorldRay(int x, int y);

    void syncShape();

    // Basic OpenGL Overrides
    void initializeGL()         override;
    void paintGL()              override;
    void resizeGL(int w, int h) override;

    // Event Listeners
    void mousePressEvent  (QMouseEvent *event) override;
    void mouseMoveEvent   (QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent       (QWheelEvent *event) override;
    void keyPressEvent    (QKeyEvent   *event) override;
    void keyReleaseEvent  (QKeyEvent   *event) override;

private slots:
    // Physics Tick
    void tick();

private:
    ARAP    m_arap;
    Shape   m_shape;
    Camera  m_camera;
    Shader *m_defaultShader;
    Shader *m_pointShader;
    std::function<void(EditorWidget*)> m_onUpdate;

    float m_movementScaling;
    float m_vertexSelectionThreshold;
    float m_vSize;

    // Timing
    QElapsedTimer m_deltaTimeProvider; // For measuring elapsed time
    QTimer        m_intervalTimer;     // For triggering timed events

    // Movement
    float m_forward;
    float m_sideways;
    float m_vertical;

    // Mouse handler stuff
    int m_lastX;
    int m_lastY;
    bool m_leftCapture;
    bool m_rightCapture;
    SelectMode m_rightClickSelectMode;
    int m_lastSelectedVertex = -1;
};
