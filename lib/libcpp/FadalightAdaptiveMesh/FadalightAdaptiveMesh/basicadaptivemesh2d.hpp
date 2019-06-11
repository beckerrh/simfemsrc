#ifndef __FadalightAdaptiveMesh_BasicAdaptiveMesh2d_h
#define __FadalightAdaptiveMesh_BasicAdaptiveMesh2d_h

// //#include  "Alat/vector.hpp"
#include  "FadalightMesh/curvedboundaryinformation.hpp"
#include  "FadalightMesh/fadalightmeshbase2d.hpp"
// #include  "FadalightMesh/quadtotri.hpp"
#include  "adaptivemeshinterface.hpp"
#include  "edge.hpp"
#include  "faceinterface.hpp"
#include  "volumeinterface.hpp"
#include  "Alat/node.hpp"
#include  "quad.hpp"
#include  "setdefs.hpp"
#include  "tri.hpp"
#include  "typedefs.hpp"
#include  <map>

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  template<int NODESPERCELL>
  class BasicAdaptiveMesh2d : virtual public AdaptiveMeshInterface, virtual public FadalightMesh::FadalightMeshBase2d<NODESPERCELL>
  {
protected:
    /// arbre des Volumes
    tree<VolumeInterface*> _Volumes;
    /// arbre des faces
    tree<FaceInterface*> _Faces;
    /// arbre des Edges
    tree<Edge*> _Edges;
    /// ensemble des noeuds
    NodeSet _Nodes;
    /// information pour les bords courbes
    // FadalightMesh::CurvedBoundaryInformation _curvedboundaries;
    /// pointeur sur un quadtotri pour un redécoupage éventuel
    // FadalightMesh::QuadToTri* _quad_to_tri_pointer;
    /// indicateur de redécoupage en quadtotri
    // bool _quad_to_tri;
    /// nom du maillage triangluaire lors d'un redécoupage en quadtotri
    // std::string _quad_to_tri_name;

    /// map affectant à chaque feuille de l'arbre des faces un numéro de feuille
    alat::IntMap _faceid2id;
    /// map affectant à chaque feuille de l'arbre des edges un numéro de feuille
    alat::IntMap _edgeid2id;
    /// map permettant de renuméroter de manière continue, les nodes du maillages
    /// important lorsque des noeuds ont disparu lors d'une procédure de déraffinement
    alat::IntMap _nodeid2id;
    /// map donnant un pointer sur une face en fonction du numéro de cell dans
    /// le maillage FadalightMesh
    std::map<int, face_pointer> _cellmap_fadalightmesh;

    /// le maillage FadalightMesh
    std::map<int, volume_pointer> _dummycellmap;

    /// dernière identité de node utilisée
    int _last_node_id;
    /// dernière identité de edge utilisée
    int _last_edge_id;
    /// dernière identité de face utilisée
    int _last_face_id;
    /// dernière identité de face utilisée
    int _last_volume_id;
    /// booleéns indiquant si les map de numérotation et le cell_map sont à jour
    bool _cell_map_ok, _numbering_ok;

    /// Ensemble des hanging edges
    EdgeSet _copy_edge;

protected:
    /// fonction renvoyant la cellule voisine de F par rapport à l'edge E
    const face_pointer _neighbour(const face_pointer F, const edge_pointer E)
    {
      if( ( *E )->face(0) != F )
      {
        return ( *E )->face(0);
      }
      else
      {
        return ( *E )->face(1);
      }
    }

    /// fonction de création d'une face
    /// utile dans le cas de raffinement mais aussi lors de la création de l'arbre des faces
    FaceInterface* _newFace(int iK)
    {
      if(getClassName() == "FadalightMesh::QuadrilateralMesh")
      {
        return new Quad(iK);
      }
      else if(getClassName() == "FadalightMesh::TriangleMesh")
      {
        return new Tri(iK);
      }
      assert(0);
    }

public:
    /// constructeur sans arguments
    BasicAdaptiveMesh2d();
    /// destructeur
    ~BasicAdaptiveMesh2d(){}

    /// Initialisation
    void BasicInit(const std::string& t);
    /// Nom du maillage en fonction du paramètre NODESPERCELLS
    std::string getClassName() const
    {
      if(NODESPERCELL == 4)
      {
        return "FadalightMesh::QuadrilateralMesh";
      }
      else if(NODESPERCELL == 3)
      {
        return "FadalightMesh::TriangleMesh";
      }
      return "BasicAdaptiveMesh2d";
    }

    /// Procédure de lecture du maillage adaptatif (arbre et numérotations)
    void readAdaptiveMesh(std::string name, std::string last_in);
    /// Procédure d'écriture du maillage adaptatif (arbre et numérotations)
    void writeAdaptiveMesh(std::string name, arma::file_type datatype = arma::arma_binary);
    /// Lecture d'un maillage de type FadalightMesh et initialisation des arbres
    void readFadalightMeshAndInitTrees(std::string name);
    /// Lecture écriture du maillage au format FadalightMesh
    void writeFadalightMesh(std::string name, arma::file_type datatype = arma::arma_binary);
    /// Mise à jour des map de numérotation
    void constructNumbering();
    /// Mise à jour du cell_map
    void constructCellMap();
    /// Mise à jour des informations sur les hanging nodes et sides
    void updateHangingInfo();
    /// Mise à jour du maillage FadalightMesh
    void reInitFadalightMesh();
    /// Nombre de nodes par cell
    int getNNodesPerCell()
    {
      return NODESPERCELL;
    }

    /// Nombre de sides par cells
    int getNSidesPerCell()
    {
      return NODESPERCELL;
    }

    /// Nombre de nodes par side
    int getNNodesPerSide()
    {
      return 2;
    }

    /// Renvoie l'ensemble des nodes
    NodeSet& getNodes()
    {
      return _Nodes;
    }

    /// Renvoie l'arbre des edges
    tree<Edge*>& getEdges()
    {
      return _Edges;
    }

    /// Renvoie l'arbre des faces
    tree<FaceInterface*>& getFaces()
    {
      return _Faces;
    }

    /// Renvoie l'arbre des faces
    tree<VolumeInterface*>& getVolumes()
    {
      return _Volumes;
    }

    /// Renvoie les informations sur les bords courbes
    // const FadalightMesh::CurvedBoundaryInformation* getCurvedBoundaries()
    // {
    //   return FadalightMesh::FadalightMeshBase2d<NODESPERCELL>::getCurvedBoundaryInformation();
    //   // return _curvedboundaries;
    // }

    /// renvoie la dernière identité de node utilisée
    int& getLastNodeId()
    {
      return _last_node_id;
    }

    /// renvoie la dernière identité d'edge utilisée
    int& getLastEdgeId()
    {
      return _last_edge_id;
    }

    /// renvoie la dernière identité de face utilisée
    int& getLastFaceId()
    {
      return _last_face_id;
    }
    /// renvoie la dernière identité de volume utilisée (0 en 2D)
    int& getLastVolumeId()
    {
      return _last_volume_id;
    }


    /// renvoie le map de renumérotation des nodes
    alat::IntMap& getNodeId2Id()
    {
      return _nodeid2id;
    }

    /// renvoie le map de renumérotation des edges actifs (feuilles)
    alat::IntMap& getEdgeId2Id()
    {
      return _edgeid2id;
    }

    /// renvoie le map de renumérotation des edges actifs (feuilles)
    alat::IntMap& getFaceId2Id()
    {
      return _faceid2id;
    }
    alat::IntMap& getVolumeId2Id()
    {
      assert(0);
      return _faceid2id;
    }

    /// renvoie le cell_map
    std::map<int, face_pointer>&  getCellMap2d()
    {
      return _cellmap_fadalightmesh;
    }
    /// renvoie le dummy cell_map pour compatibilité
    std::map<int, volume_pointer>&  getCellMap3d()
    {
      return _dummycellmap;
    }

    /// écriture au format vtk pour visualisation
    void writeVtk(std::string name)
    {
      FadalightMesh::FadalightMeshBase2d<NODESPERCELL>::writeVtk(name);
    }

    // /// écriture au format Ensight pour visualisation
    // void writeEnsightGeometry(std::string name)
    // {
    //   FadalightMesh::FadalightMeshBase2d<NODESPERCELL>::writeEnsightGeometry(name);
    // }

    /// renvoie l'indicateur de mise à jour du cell_map
    bool&  getCellMapOk()
    {
      return _cell_map_ok;
    }

    /// renvoie l'indicateur de mise à jour des numérotations
    bool&  getNumberingOk()
    {
      return _numbering_ok;
    }

    /// Récriture de la fonction pour permettre son appel dans les outils d'adaptation
    FadalightMesh::GeometryObject* getGeometryObject(std::string name)
    {
      return FadalightMesh::FadalightMeshBase2d<NODESPERCELL>::getGeometryObject(name);
    }

    /// Récriture de la fonction pour permettre son appel dans les outils d'adaptation
    bool geometryObjectExists(std::string name) const
    {
      return FadalightMesh::FadalightMeshBase2d<NODESPERCELL>::geometryObjectExists(name);
    }

    /// Récriture de la fonction pour permettre son appel dans les outils d'adaptation
    void createGeometryObject(std::string name)
    {
      FadalightMesh::FadalightMeshBase2d<NODESPERCELL>::createGeometryObject(name);
    }

    /// affecte l'indicateur pour un découpage en quadtotri
    // void setQuadToTri(bool quadtotri)
    // {
    //   if(getClassName() == "FadalightMesh::QuadrilateralMesh")
    //   {
    //     _quad_to_tri = quadtotri;
    //     _quad_to_tri_pointer = new FadalightMesh::QuadToTri;
    //   }
    //   else
    //   {
    //     _quad_to_tri = false;
    //   }
    // }

    // void setQuadToTriFileName(std::string quadtotriname)
    // {
    //   _quad_to_tri_name = quadtotriname;
    // }
    //
    // bool quadToTri()
    // {
    //   return _quad_to_tri;
    // }

    // FadalightMesh::QuadToTri* getQuadToTri()
    // {
    //   return _quad_to_tri_pointer;
    // }
    //
    // std::string getQuadToTriFileName()
    // {
    //   return _quad_to_tri_name;
    // }
  };
  typedef BasicAdaptiveMesh2d<3> TriangleAdaptiveMesh;
  typedef BasicAdaptiveMesh2d<4> QuadrilateralAdaptiveMesh;
}

/*--------------------------------------------------------------------------*/

#endif
