/*
 * MASTER 2 GIG - Modèles Géométriques - 2021-2022
 * Auteurs : Gabriel PAREL & Samantha THIEL
 * Merci à A. POLETTE et A. BAC pour la base de code.
 *
 */

#include "mainwindow.h"
#include "ui_mainwindow.h"

bool MainWindow::are1Neighbors(const MyMesh::VertexHandle vh1, const MyMesh::VertexHandle vh2)
{
    bool areNeighbors = false;
    //circulateur sur le 1-voisinage
    for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh1);vv_it.is_valid();++vv_it)
    {
        if (*vv_it == vh2) areNeighbors = true;
    }
    return areNeighbors;
}

float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    float aire = 0.;
    FaceHandle fh = _mesh->face_handle(faceID);  //face courante

    //iterator sur les sommets de la face courante
    MyMesh::FaceVertexIter fv_it = _mesh->fv_iter(fh);

    //points du triangle courant ABC
    const MyMesh::Point& A = _mesh->point(*fv_it);  ++fv_it;
    const MyMesh::Point& B = _mesh->point(*fv_it);  ++fv_it;
    const MyMesh::Point& C = _mesh->point(*fv_it);

    //aire du triangle
    aire += (((B - A)%(C - A)).norm())/2;

    return aire;
}

std::vector<QVector3D> MainWindow::computeLaplace(MyMesh *_mesh, int typeIdx)
{
    float lambda = 0., h = 0.; //pour flou de diffusion
    std::vector<QVector3D> matrixOp;
    std::vector<int> idx; //stockage des indices
    QVector<float> areas; //D

    // parcours des sommets
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert !=
         _mesh->vertices_end(); curVert++)
    {
        auto point = _mesh->point(*curVert);
        QVector3D v = QVector3D(point[0], point[1], point[2]);

        //calcul de la composante angulaire
        //chercher le point suivant (next) et précédent (previous)
        //pour le calcul d'angles
        float cot = 0., area = 0.;
        std::vector<float> result(3);
        for (MyMesh::VertexOHalfedgeCWIter vheh_cwit = _mesh->voh_cwiter(*curVert);
             vheh_cwit.is_valid(); vheh_cwit ++)
        {
            //on récupère l'arête courante, suivante, précédente
             HalfedgeHandle oheh = _mesh->opposite_halfedge_handle(*vheh_cwit);
             HalfedgeHandle heh1 = _mesh->next_halfedge_handle(oheh);
             HalfedgeHandle heh2 = _mesh->next_halfedge_handle(*vheh_cwit);

             FaceHandle fh = mesh.face_handle(heh1);
             area += faceArea(&mesh, fh.idx());

             //on récupère le point courant, suivant, précédent
             auto point_i = _mesh->point( _mesh->to_vertex_handle(*vheh_cwit));
             auto next_point_i = _mesh->point( _mesh->to_vertex_handle(heh1));
             auto prev_point_i = _mesh->point( _mesh->to_vertex_handle(heh2));

             //si approche cotangentielle, on calcule les angles alpha et beta afin
             //de déterminer les cotangentes de ceux-ci
             if (typeIdx == 0)
             {
                 lambda = 0.01, h = 0.01; //bons résultats avec ces valeurs
                 QVector3D vi = QVector3D(point_i[0], point_i[1], point_i[2]);
                 QVector3D nvi = QVector3D(next_point_i[0], next_point_i[1], next_point_i[2]);
                 QVector3D pvi = QVector3D(prev_point_i[0], prev_point_i[1], prev_point_i[2]);

                 //qDebug() << "vi, nvi, pvi" << vi << nvi << pvi; //pour vérification

                 //calcul des vecteurs d'intérêt; attention à l'orientation !
                 QVector3D nviv = v-nvi;
                 QVector3D pviv = v-pvi;
                 QVector3D nvivi = vi - nvi;
                 QVector3D pvivi = vi - pvi;

                 float alpha = 0., beta = 0.;
                 if (QVector3D::dotProduct(nviv, nvivi) != 0)
                     alpha = acos(QVector3D::dotProduct(nviv, nvivi)/
                                  (nviv.length()*nvivi.length()));
                 if (QVector3D::dotProduct(pviv, pvivi) != 0)
                     beta = acos(QVector3D::dotProduct(pviv, pvivi)/
                                 (pviv.length()*pvivi.length()));

                 //qDebug() << "alpha, beta" << alpha << beta; //pour vérification

                 cot = ((float)(1/tan(alpha))+(float)(1/tan(beta)));

                 //result correspond à la composante cotangentielle du calcul
                 result[0] += (float) cot*(point_i[0] - point[0]);
                 result[1] += (float) cot*(point_i[1] - point[1]);
                 result[2] += (float) cot*(point_i[2] - point[2]);
             }
             //si approche uniforme, on ne calcule pas d'angles
             else if (typeIdx == 1)
             {
                 lambda = 0.15, h = 0.5; //bons résultats avec ces valeurs
                 result[0] += (point_i[0] - point[0]);
                 result[1] += (point_i[1] - point[1]);
                 result[2] += (point_i[2] - point[2]);
             }
        }

        std::vector<float> laplaceOp(3);
        //si approche cotangentielle, on calcule l'aire barycentrique
        if (typeIdx == 0 || typeIdx == 2)
        {
            // calcul de l'aire barycentrique
            float barycentricArea = (float) 1/3*area;
            areas.push_back(1/(2*barycentricArea));

            //calcul du laplacien cotangentiel
            laplaceOp[0] = (1/(2*barycentricArea))*result[0];
            laplaceOp[1] = (1/(2*barycentricArea))*result[1];
            laplaceOp[2] = (1/(2*barycentricArea))*result[2];
        }

        else if (typeIdx == 1)
        {
            //calcul du laplacien uniforme
            laplaceOp[0] = result[0];
            laplaceOp[1] = result[1];
            laplaceOp[2] = result[2];
        }

        matrixOp.push_back(QVector3D(laplaceOp[0], laplaceOp[1], laplaceOp[2]));
        idx.push_back(curVert->idx());

        /// APPROCHE MATRICIELLE; non terminée
        /// voir commentaires associés pour plus d'infos

        /*//si approche matricielle, on calcule les matrices M et D
        if (typeIdx == 2)
        {
            QVector<QVector<float>> L, D, M;
            for (int i = 0; i<(int)matrixOp.size(); i++)
            {
                auto vhI = _mesh->vertex_handle((idx[i])); //vi
                QVector<float> rowL, rowD;
                for(int j = 0; j<(int)matrixOp.size(); j++)
                {
                    auto vhJ = _mesh->vertex_handle((idx[j])); //vj
                    if (i==j) //diagonale
                    {
                        //rowL.push_back(); //-somme des wi,k
                        rowD.push_back(areas[idx[i]]); //1/2*A(vi)
                    }
                    else
                    {
                        //si vj est dans le 1-voisinage de vi
                        if(are1Neighbors(vhI, vhJ))
                        {
                            //rowL.push_back(); //wi,j
                            break;
                        }
                        rowL.push_back(0);
                        rowD.push_back(0);
                    }
                }
                L.push_back(rowL);
                D.push_back(rowD);
            }

            //L = DM
            //besoin de multiplier les matrices entre elles...
            //la structure de données QVector<QVector> n'est pas très adéquate
            //il faudrait utiliser la bibliothèque Eigen ou OpenCV
        } */
    }

    //flou de diffusion; on applique directement à chaque point du
    //maillage courant l'opérateur de Laplace
    for (int i = 0; i < (int) idx.size(); i++)
    {
        auto point = _mesh->point(_mesh->vertex_handle((idx[i])));
        point[0] += (float) h*lambda*matrixOp[i][0];
        point[1] += (float) h*lambda*matrixOp[i][1];
        point[2] += (float) h*lambda*matrixOp[i][2];

        _mesh->set_point(_mesh->vertex_handle((idx[i])), point);
    }

    displayMesh(_mesh); //affichage du nouveau maillage
    return matrixOp;
}

/* **** début de la partie boutons et IHM **** */


void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "",
                                                    tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

void MainWindow::on_pushButton_courbures_clicked()
{
    mesh.update_normals();
    //mesh.request_vertex_colors() ;

    int typeIdx = ui->idxComboBox->currentIndex();
    computeLaplace(&mesh, typeIdx);

    displayMesh(&mesh/*,DisplayMode::VertexColorShading*/);
}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */

// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert !=
         _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace !=
         _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge !=
         _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, DisplayMode mode)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(mode == DisplayMode::TemperatureMap)
    {
        QVector<float> values;
        for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert !=
             _mesh->vertices_end(); curVert++)
            values.append(fabs(_mesh->data(*curVert).value));
        qSort(values);

        float range = values.at(values.size()*0.8);

        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::Normal)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::ColorShading)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if (mode == DisplayMode::VertexColorShading)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

