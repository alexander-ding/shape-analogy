#include "analogywidget.h"

#include <QApplication>
#include <QKeyEvent>
#include <QTime>
#include <iostream>

#include "mesh.h"
#include "analogy.h"

#define SPEED 1.5
#define ROTATE_SPEED 0.0025

using namespace std;

AnalogyWidget::AnalogyWidget(Analogy analogy, QWidget *parent) : m_analogy(analogy),
                                                       QOpenGLWidget(parent),
                                                       m_deltaTimeProvider(),
                                                       m_intervalTimer(),
                                                       m_camera(),
                                                       m_shader(),
                                                       m_forward(),
                                                       m_b(),
                                                       m_bPrime(),
                                                       m_sideways(),
                                                       m_vertical(),
                                                       m_lastX(),
                                                       m_lastY(),
                                                       m_capture(false)
{
    // AnalogyWidget needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    QApplication::setOverrideCursor(Qt::ArrowCursor);

    // AnalogyWidget needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // Function tick() will be called once per interva
    connect(&m_intervalTimer, SIGNAL(timeout()), this, SLOT(tick()));
}

AnalogyWidget::~AnalogyWidget()
{
    if (m_shader != nullptr)
        delete m_shader;
}

void AnalogyWidget::syncAnalogy()
{
    this->m_analogy.setPrevFrameBPrimeVerts(this->m_analogy.getBPrime().getVertices());
    this->m_analogy.computeBPrime(1000);
    Eigen::MatrixXf colors = (m_analogy.getBTargetNormals().array() + 1.f) / 2.f;
    Mesh &b = this->m_analogy.getB();
    m_b.init(b.getVertices(), b.computeVertexNormals(), b.getFaces(), colors);
    Eigen::Affine3f bTranslation = Eigen::Affine3f::Identity() * Eigen::Translation3f(-2, 0, 0);
    m_b.setModelMatrix(bTranslation);

    Mesh &bPrime = this->m_analogy.getBPrime();
    m_bPrime.init(bPrime.getVertices(), bPrime.computeVertexNormals(), bPrime.getFaces());
}


void AnalogyWidget::reset() {
    this->m_analogy.setAPrime(this->m_analogy.getAPrimeCache());
    this->m_analogy.setBPrime(this->m_analogy.getB());
    syncAnalogy();
}

void AnalogyWidget::undo() {
    this->m_analogy.undo();
    syncAnalogy();
}
// ================== Basic OpenGL Overrides

void AnalogyWidget::initializeGL()
{
    // Initialize GL extension wrangler
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK)
        fprintf(stderr, "Error while initializing GLEW: %s\n", glewGetErrorString(err));
    fprintf(stdout, "Successfully initialized GLEW %s\n", glewGetString(GLEW_VERSION));

    // Set clear color to white
    glClearColor(1, 1, 1, 1);

    // Enable depth-testing and backface culling
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Initialize the shader and simulation
    m_shader = new Shader(":/resources/shaders/shader.vert", ":/resources/shaders/shader.frag");

    this->m_analogy.computeBPrime(1000);
    Eigen::MatrixXf colors = (m_analogy.getBTargetNormals().array() + 1.f) / 2.f;
    Mesh &b = this->m_analogy.getB();
    m_b.init(b.getVertices(), b.computeVertexNormals(), b.getFaces(), colors);
    Eigen::Affine3f bTranslation = Eigen::Affine3f::Identity() * Eigen::Translation3f(-2, 0, 0);
    m_b.setModelMatrix(bTranslation);

    Mesh &bPrime = this->m_analogy.getBPrime();
    m_bPrime.init(bPrime.getVertices(), bPrime.computeVertexNormals(), bPrime.getFaces());

    // --------
    Vector3f coeffMin, coeffMax;
    MatrixXf& vertices = bPrime.getVertices();
    coeffMin = vertices.rowwise().minCoeff();
    coeffMax = vertices.rowwise().maxCoeff();


    Vector3f center = (coeffMax + coeffMin) / 2.0;
    center = center - Vector3f(coeffMax(0) * 2, 0, 0);
    float extentLength = (coeffMax - coeffMin).norm();


    // Note for maintainers: Z-up
    float fovY = 120;
    float nearPlane = 0.001f;
    float farPlane = 4 * extentLength;

    // Initialize camera with a reasonable transform
    Eigen::Vector3f eye = center - Eigen::Vector3f::UnitZ() * extentLength * 2.5 + Eigen::Vector3f::UnitY() * 1.5;
    Eigen::Vector3f target = center;
    m_camera.lookAt(eye, target);
    m_camera.setOrbitPoint(target);
    m_camera.setPerspective(120, width() / static_cast<float>(height()), nearPlane, farPlane);

    // ------

    m_deltaTimeProvider.start();
    m_intervalTimer.start(1000. / FRAMES_TO_AVERAGE);
}

void AnalogyWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    m_shader->bind();
    m_shader->setUniform("proj", m_camera.getProjection());
    m_shader->setUniform("view", m_camera.getView());
    this->m_bPrime.draw(m_shader);
    this->m_b.draw(m_shader);

    m_shader->unbind();
}

void AnalogyWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    m_camera.setAspect(static_cast<float>(w) / h);
}

// ================== Event Listeners

void AnalogyWidget::mousePressEvent(QMouseEvent *event)
{
    m_capture = true;
    m_lastX = event->position().x();
    m_lastY = event->position().y();
}

void AnalogyWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (!m_capture)
        return;

    int currX = event->position().x();
    int currY = event->position().y();

    int deltaX = currX - m_lastX;
    int deltaY = currY - m_lastY;

    if (deltaX == 0 && deltaY == 0)
        return;

    m_camera.rotate(deltaY * ROTATE_SPEED,
                    -deltaX * ROTATE_SPEED);

    m_lastX = currX;
    m_lastY = currY;
}

void AnalogyWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_capture = false;
}

void AnalogyWidget::wheelEvent(QWheelEvent *event)
{
    float zoom = 1 - event->pixelDelta().y() * 0.1f / 120.f;
    m_camera.zoom(zoom);
}

void AnalogyWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat())
        return;

    switch (event->key())
    {
    case Qt::Key_W:
        m_forward = SPEED;
        break;
    case Qt::Key_S:
        m_forward = -SPEED;
        break;
    case Qt::Key_A:
        m_sideways = -SPEED;
        break;
    case Qt::Key_D:
        m_sideways = SPEED;
        break;
    case Qt::Key_F:
        m_vertical = -SPEED;
        break;
    case Qt::Key_R:
        m_vertical = SPEED;
        break;
    case Qt::Key_C:
        m_camera.toggleIsOrbiting();
        break;
    case Qt::Key_Escape:
        QApplication::quit();
    }
}

void AnalogyWidget::keyReleaseEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat())
        return;

    switch (event->key())
    {
    case Qt::Key_W:
        m_forward = 0;
        break;
    case Qt::Key_S:
        m_forward = 0;
        break;
    case Qt::Key_A:
        m_sideways = 0;
        break;
    case Qt::Key_D:
        m_sideways = 0;
        break;
    case Qt::Key_F:
        m_vertical = 0;
        break;
    case Qt::Key_R:
        m_vertical = 0;
        break;
    }
}

// ================== Physics Tick

void AnalogyWidget::tick()
{
    float deltaSeconds = m_deltaTimeProvider.restart() / 1000.f;

    // Move camera
    auto look = m_camera.getLook();
    look.y() = 0;
    look.normalize();
    Eigen::Vector3f perp(-look.z(), 0, look.x());
    Eigen::Vector3f moveVec = m_forward * look.normalized() + m_sideways * perp.normalized() + m_vertical * Eigen::Vector3f::UnitY();
    moveVec *= deltaSeconds;
    m_camera.move(moveVec);

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}


