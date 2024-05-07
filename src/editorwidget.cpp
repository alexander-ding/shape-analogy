#include "editorwidget.h"
#include "arap.h"
#include "settings.h"

#include <QApplication>
#include <QKeyEvent>
#include <iostream>

#define SPEED 1.5
#define ROTATE_SPEED 0.0025
#define PUSH_SPEED 0.0025
#define MAX_RADIUS 0.5f

using namespace std;
using namespace Eigen;

inline MatrixXf vectorsToMatrix(const std::vector<Eigen::Vector3f>& vertices) {
    MatrixXf m(3, vertices.size());
    for (size_t i = 0; i < vertices.size(); i++) {
        m.col(i) = vertices[i];
    }
    return m;
}

inline std::vector<Eigen::Vector3f> matrixToVectors(const MatrixXf& matrix) {
    std::vector<Vector3f> v(matrix.cols());
    for (size_t i = 0; i < matrix.cols(); i++) {
        v[i] = matrix.col(i);
    }
    return v;
}


EditorWidget::EditorWidget(Mesh& mesh, QWidget *parent) :
    QOpenGLWidget(parent),
    m_arap(mesh),
    m_onUpdate(nullptr),
    m_camera(),
    m_defaultShader(),
    m_pointShader(),
    m_vSize(),
    m_movementScaling(),
    m_vertexSelectionThreshold(),
    // Movement
    m_deltaTimeProvider(),
    m_intervalTimer(),
    // Timing
    m_forward(),
    m_sideways(),
    m_vertical(),
    // Mouse handler stuff
    m_lastX(),
    m_lastY(),
    m_leftCapture(false),
    m_rightCapture(false),
    m_rightClickSelectMode(SelectMode::None),
    m_lastSelectedVertex(-1),
    m_mesh(mesh)
{
    // EditorWidget needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    QApplication::setOverrideCursor(Qt::ArrowCursor);

    // EditorWidget needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // Function tick() will be called once per interva
    connect(&m_intervalTimer, SIGNAL(timeout()), this, SLOT(tick()));
}

EditorWidget::~EditorWidget()
{
    if (m_defaultShader != nullptr)
        delete m_defaultShader;
    if (m_pointShader != nullptr)
        delete m_pointShader;
}

void EditorWidget::onUpdate(std::function<void(EditorWidget*)> f) {
    this->m_onUpdate = f;
}

// ================== Basic OpenGL Overrides

void EditorWidget::initializeGL()
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

    // Initialize shaders
    m_defaultShader = new Shader(":resources/shaders/arapShader.vert", ":resources/shaders/arapShader.frag");
    m_pointShader = new Shader(":resources/shaders/anchorPoint.vert", ":resources/shaders/anchorPoint.geom", ":resources/shaders/anchorPoint.frag");

    // Initialize ARAP, and get parameters needed to decide the camera position, etc
    Vector3f coeffMin, coeffMax;
    m_arap.init(coeffMin, coeffMax);
    Mesh& mesh = this->m_arap.getMesh();
    m_shape.init(mesh.getVertices(), mesh.computeVertexNormals(), mesh.getFaces());
    syncShape();

    Vector3f center = (coeffMax + coeffMin) / 2.0f;
    float extentLength = (coeffMax - coeffMin).norm();

    // Screen-space size of vertex points
    m_vSize = 0.005 * extentLength;

    // Scale all movement by this amount
    m_movementScaling = extentLength * 0.5;

    // When raycasting, select closest vertex within this distance
    m_vertexSelectionThreshold = extentLength * 0.025;

    // Note for maintainers: Z-up
    float fovY = 120;
    float nearPlane = 0.001f;
    float farPlane = 4 * extentLength;

    // Initialize camera with a reasonable transform
    Eigen::Vector3f eye = center - Eigen::Vector3f::UnitZ() * extentLength * 2;
    Eigen::Vector3f target = center;
    m_camera.lookAt(eye, target);
    m_camera.setOrbitPoint(target);
    m_camera.setPerspective(120, width() / static_cast<float>(height()), nearPlane, farPlane);

    m_deltaTimeProvider.start();
    m_intervalTimer.start(1000 / 60);
}

void EditorWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_defaultShader->bind();
    m_defaultShader->setUniform("proj", m_camera.getProjection());
    m_defaultShader->setUniform("view", m_camera.getView());
    m_shape.draw(m_defaultShader);
    m_defaultShader->unbind();

    glClear(GL_DEPTH_BUFFER_BIT);

    m_pointShader->bind();
    m_pointShader->setUniform("proj", m_camera.getProjection());
    m_pointShader->setUniform("view", m_camera.getView());
    m_pointShader->setUniform("vSize", m_vSize);
    m_pointShader->setUniform("width", width());
    m_pointShader->setUniform("height", height());
    m_shape.draw(m_pointShader, GL_POINTS);
    m_pointShader->unbind();
}

void EditorWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    m_camera.setAspect(static_cast<float>(w) / h);
}

// ================== Event Listeners

Eigen::Vector3f EditorWidget::transformToWorldRay(int x, int y)
{
    Eigen::Vector4f clipCoords = Eigen::Vector4f(
        (float(x) / width()) * 2.f - 1.f,
        1.f - (float(y) / height()) * 2.f,
        -1.f,
        1.f);

    Eigen::Vector4f transformed_coords = m_camera.getProjection().inverse() * clipCoords;
    transformed_coords = Eigen::Vector4f(transformed_coords.x(), transformed_coords.y(), -1.f, 0.f);
    transformed_coords = m_camera.getView().inverse() * transformed_coords;

    return Eigen::Vector3f(transformed_coords.x(), transformed_coords.y(), transformed_coords.z()).normalized();
}

void EditorWidget::mousePressEvent(QMouseEvent *event)
{
    // Get current mouse coordinates
    const int currX = event->position().x();
    const int currY = event->position().y();

    // Get closest vertex to ray
    const Vector3f ray = transformToWorldRay(currX, currY);

    const Vector3f start = m_camera.getPosition();

    const int closest_vertex = m_arap.getClosestVertex(start, ray, m_vertexSelectionThreshold);

    // Switch on button
    switch (event->button())
    {
    case Qt::MouseButton::RightButton:
    {
        // Capture
        m_rightCapture = true;

        switch (settings.mode) {
        case ARAP:
            // Anchor/un-anchor the vertex
            m_rightClickSelectMode = m_arap.select(closest_vertex);
            break;
        case PUSH: {
            m_arap.setPrevFrameVerts(m_arap.getVerts());
            Intersection shapeIntersect = m_arap.intersectMesh(start, ray);
            if (!shapeIntersect.hit) {
                break;
            }

            // std::cout << "hit at: " << shapeIntersect.position << std::endl;

            m_push_position = shapeIntersect.position;
            m_push_direction = -shapeIntersect.normal;
            m_push = PUSH_SPEED;
            break;
        }
        case HAMMER: {
            m_arap.setPrevFrameVerts(m_arap.getVerts());
            Intersection shapeIntersect = m_arap.intersectMesh(start, ray);
            if (!shapeIntersect.hit) {
                break;
            }

            float normalizedRadius = settings.radius / 100.f;
            m_arap.hammer(shapeIntersect.position, normalizedRadius);
            break;
        }
        default:
            break;
        }

        syncShape();

        break;
    }
    case Qt::MouseButton::LeftButton:
    {
        // Capture
        m_leftCapture = true;
        switch (settings.mode) {
        case ARAP:
            m_arap.setPrevFrameVerts(m_arap.getVerts());
            // Select this vertex
            m_lastSelectedVertex = closest_vertex;
        default:
            break;
        }
        break;
    }
    default:
        break;
    }

    // Set last mouse coordinates
    m_lastX = currX;
    m_lastY = currY;
}

void EditorWidget::mouseMoveEvent(QMouseEvent *event)
{
    // Return if neither mouse button is currently held down
    if (!(m_leftCapture || m_rightCapture))
    {
        return;
    }

    // Get current mouse coordinates
    const int currX = event->position().x();
    const int currY = event->position().y();

    // Find ray
    const Vector3f ray = transformToWorldRay(event->position().x(), event->position().y());
    Vector3f pos;

    // If right is held down
    if (settings.mode == ARAP && m_rightCapture)
    {
        // Get closest vertex to ray
        const int closest_vertex = m_arap.getClosestVertex(m_camera.getPosition(), ray, m_vertexSelectionThreshold);

        // Anchor/un-anchor the vertex
        if (m_rightClickSelectMode == SelectMode::None)
        {
            m_rightClickSelectMode = m_arap.select(closest_vertex);
        }
        else
        {
            m_arap.selectWithSpecifiedMode(closest_vertex, m_rightClickSelectMode);
        }

        this->syncShape();

        return;
    }

    // If the selected point is an anchor point
    if (settings.mode == ARAP && m_lastSelectedVertex != -1 && m_arap.getAnchorPos(m_lastSelectedVertex, pos, ray, m_camera.getPosition()))
    {
        // Move it
        m_arap.move(m_lastSelectedVertex, pos);
        this->syncShape();
    }
    else
    {
        if (settings.mode == PUSH && m_rightCapture) {
            m_lastX = currX;
            m_lastY = currY;
            return;
        }
        // Rotate the camera
        const int deltaX = currX - m_lastX;
        const int deltaY = currY - m_lastY;
        if (deltaX != 0 || deltaY != 0)
        {
            m_camera.rotate(deltaY * ROTATE_SPEED, -deltaX * ROTATE_SPEED);
        }
    }

    // Set last mouse coordinates
    m_lastX = currX;
    m_lastY = currY;
}

void EditorWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_leftCapture = false;
    m_lastSelectedVertex = -1;
    bool wasPush = m_push;
    m_push = 0;

    m_rightCapture = false;
    m_rightClickSelectMode = SelectMode::None;
    bool wasARAP = this->m_arap.invalidate();
    if (wasPush || wasARAP) {
        // std::cout << "Committing update" << std::endl;
        this->m_arap.commitUpdate(this->m_arap.getVerts());
    }
    if (this->m_arap.getIsUnsynced() && this->m_onUpdate) {
        // cout << "updating" << endl;
        this->m_arap.setIsUnsynced(false);
        syncShape();
        this->m_onUpdate(this);
    }
}

void EditorWidget::wheelEvent(QWheelEvent *event)
{
    float zoom = 1 - event->pixelDelta().y() * 0.1f / 120.f;
    m_camera.zoom(zoom);
}

void EditorWidget::keyPressEvent(QKeyEvent *event)
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
    case Qt::Key_Equal:
        m_vSize *= 11.0f / 10.0f;
        break;
    case Qt::Key_Minus:
        m_vSize *= 10.0f / 11.0f;
        break;
    case Qt::Key_Escape:
        QApplication::quit();

    }
}

void EditorWidget::keyReleaseEvent(QKeyEvent *event)
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

void EditorWidget::syncShape() {
    Mesh& mesh = this->m_arap.getMesh();
    Eigen::MatrixXf colors(3, mesh.getVertices().cols());
    colors.setConstant(0.3);
    auto& anchors = this->m_arap.getAnchors();
    for (size_t i = 0; i < colors.cols(); i++) {
        if (anchors.find(i) == anchors.end()) {
            colors.col(i) = Vector3f(0, 0.8, 0.8);
        } else {
            colors.col(i) = Vector3f(1, 1, 0);
        }
    }
    this->m_shape.setVertices(mesh.getVertices(), mesh.computeVertexNormals(), colors);
}

void EditorWidget::reset() {
    // return mesh to sphere
    this->m_arap.setMesh(m_mesh);
    // clear anchor
    this->m_arap.clearAnchors();
    // update visualization
    syncShape();
}

void EditorWidget::undo() {
    this->m_arap.undo();
    syncShape();
}

// ================== Physics Tick

void EditorWidget::tick()
{
    float deltaSeconds = m_deltaTimeProvider.restart() / 1000.f;

    // Move camera
    auto look = m_camera.getLook();
    look.y() = 0;
    look.normalize();
    Eigen::Vector3f perp(-look.z(), 0, look.x());
    Eigen::Vector3f moveVec = m_forward * look.normalized() + m_sideways * perp.normalized() + m_vertical * Eigen::Vector3f::UnitY();
    moveVec *= m_movementScaling;
    moveVec *= deltaSeconds;
    m_camera.move(moveVec);

    if (settings.mode == PUSH && m_push > 0) {

        m_push_position += m_push * m_push_direction;
        m_arap.push(m_push_position, m_push_direction);
        syncShape();
    }

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
