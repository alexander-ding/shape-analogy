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

    void reset();
    void undo();

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
    Mesh m_mesh;

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

    // Pushing
    float m_push;
    Vector3f m_push_direction;
    Vector3f m_push_position;

    // Hammer
    Vector3f m_hammer_position;

    // Mouse handler stuff
    int m_lastX;
    int m_lastY;
    bool m_leftCapture;
    bool m_rightCapture;
    SelectMode m_rightClickSelectMode;
    int m_lastSelectedVertex = -1;
};
