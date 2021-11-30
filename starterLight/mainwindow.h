/*
 * MASTER 2 GIG - Modèles Géométriques - 2021-2022
 * Auteurs : Gabriel PAREL & Samantha THIEL
 * Merci à A. POLETTE et A. BAC pour la base de code.
 *
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QVector3D>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness; float value; Color faceShadingColor;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;


enum DisplayMode {Normal, TemperatureMap, ColorShading, VertexColorShading};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    /**
     * @brief MainWindow::faceArea : aire d'une face du maillage
     * @param _mesh : le maillage ouvert dans la fenêtre
     * @param faceID : face dont on calcule l'aire
     * @return l'aire de la face courante
     */
    float faceArea(MyMesh* _mesh, int faceID);

    ///
    /// \brief MainWindow::computeLaplace : lisse le maillage courant par
    /// application d'un opérateur de Laplace cotangentiel ou uniforme
    /// \param _mesh : le maillage courant
    /// \param typeIdx : 0 pour approche cotangentielle,
    /// 1 pour approche uniforme
    /// \return un vector des coordonnées après applications du laplacien
    ///
    std::vector<QVector3D> computeLaplace(MyMesh *_mesh, int typeIdx);

    ///
    /// \brief are1Neighbors : retourne vrai si les sommets vh1 et vh2
    /// sont voisins (1-voisinage), faux sinon
    /// \param vh1 : sommet 1
    /// \param vh2 : sommet 2
    /// \return vrai si vh2 appartient au 1-voisinage de vh2, faux sinon
    ///
    bool are1Neighbors(const MyMesh::VertexHandle vh1, const MyMesh::VertexHandle vh2);

    void displayMesh(MyMesh *_mesh, DisplayMode mode = DisplayMode::Normal);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();

    void on_pushButton_courbures_clicked(); //lissage

private:
    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
