
#if !defined(_pcs_pof_bspfuncs_h_)
#define _pcs_pof_bspfuncs_h_

#include <memory>
#include <unordered_map>

#include "pcs_file.h"
#include "pcs_file_dstructs.h"
#include "POFHandler.h"
#include "BSPDataStructs.h"

vector3d POFTranslate(vector3d v);

struct bsp_vert
{
	vector3d point;
	std::vector<vector3d> norms;
};

bool operator==(const bsp_vert &one, const bsp_vert &two);

bool Neighbor(pcs_shield_triangle &face1, pcs_shield_triangle &face2);

pcs_polygon RebuildCleanPolygon(pcs_polygon &src);

//-----------------------------------------------
// mesh and math functions
//-----------------------------------------------
void TriangulateMesh(std::vector<pcs_polygon> &polygons);
void SplitPolygon(std::vector<pcs_polygon> &polygons, int polynum, vector3d plane_point, vector3d plane_normal, std::vector<pcs_polygon> &newpolys);
void SplitIntersecting(std::vector<pcs_polygon> &polygons, vector3d plane_point, vector3d plane_normal);
float FindIntersection(vector3d &intersect, vector3d p1, vector3d p2, vector3d plane_point, vector3d plane_normal);
float DistanceToPlane(vector3d point, vector3d plane_point, vector3d plane_normal);
bool InFrontofPlane(vector3d point, vector3d plane_point, vector3d plane_normal);


//-----------------------------------------------
// code for writing a BSP from PMF -- new style
//-----------------------------------------------

enum node_type { SPLIT, POLY, INVALID };

struct bsp_tree_node
{
	bsp_tree_node()
		:Type(INVALID), used(false), counted(false)
	{};
	node_type Type;

	// used by SPLIT only
	vector3d normal;
	vector3d point;

	// used by POLY only
	std::vector<int> poly_num;

	// used by BOTH
	vector3d bound_min;
	vector3d bound_max;

	// tree
	std::unique_ptr<bsp_tree_node> front; // front recruse
	std::unique_ptr<bsp_tree_node> back; // back recruse

	// safety variables - counted is true when calculate size has hit a node, used is true when writing has hit a node
	bool used;
	bool counted;
};


//-----------------------------------------------
// interfacable functions for generating the BSP tree/SLDC tree
//-----------------------------------------------

//shared
std::unique_ptr<bsp_tree_node> MakeTree(std::vector<pcs_polygon> &polygons, vector3d &Max, vector3d &Min);
void DebugPrintTree(bsp_tree_node* root, std::ostream &out);

// BSP tree functions
#define BSP_NOERRORS				0

// --------- errors raised while packing --------- 
// raised if it detects less than 7 bytes free in buffer at start of a call
#define BSP_PACK_PREOVERFLOW		0x00000001
// raised if it detects a double usage of a node
#define BSP_PACK_DOUBLEUSE			0x00000002
// raised if it detects an attempt to use an uncounted node
#define BSP_PACK_UNCOUNTED			0x00000004
// raised if it detects an overflow in a polygon after-writing
#define BSP_PACK_POLYOVERFLOW		0x00000008
// raised if it detects an overflow in a split after-writing
#define BSP_PACK_SPLITOVERFLOW		0x00000010
// raised if it detects an overflow in a polygon before writing
#define BSP_PACK_PREPOLYOVERFLOW	0x00000020
// raised if it detects an overflow in a split before writing
#define BSP_PACK_PRESPLITOVERFLOW	0x00000040

int CalculateTreeSize(bsp_tree_node* root, std::vector<pcs_polygon> &polygons);
int PackTreeInBSP(bsp_tree_node* root, int offset, char *buffer, std::vector<pcs_polygon> &polygons,
	std::unordered_map<vector3d, int> &norms, std::unordered_map<vector3d, int> &verts, BSP_DefPoints &dpnts, vector3d geo_center, int buffsize, int &error_flags);

// closely related functions for SLDC meshes

int CalcSLDCTreeSize(bsp_tree_node* root);
int PackTreeInSLDC(bsp_tree_node* root, int offset, char *buffer, int bufsz);


//-----------------------------------------------
// aux functions for generating the BSP tree
//-----------------------------------------------
std::unique_ptr<bsp_tree_node> GenerateTreeRecursion(std::vector<pcs_polygon> &polygons, std::vector<int>&);

bool Bisect(const vector3d& cmax, const vector3d& cmin,
	vector3d &p_point, vector3d &p_norm,
	const std::vector<pcs_polygon>& polygons,
	std::vector<int>& contained,
	std::vector<int>& front,
	std::vector<int>& back,
	vector3d *centera = NULL, vector3d *centerb = NULL);
vector3d PolygonCenter(pcs_polygon &polygon);
void BoundPolygon(vector3d &Max, vector3d &Min, int polygon, std::vector<pcs_polygon> &polygons);
void MakeBound(vector3d &Max, vector3d &Min, std::vector<int> &polylist, std::vector<pcs_polygon> &polygons);

void AddIfNotInList(std::vector<pcs_vertex> &list, pcs_vertex &point);


//-----------------------------------------------
// functions for writing a BSP from PMF -- old style
//-----------------------------------------------

void MakeTmapPoly(BSP_TmapPoly &dst, pcs_polygon &src, std::unordered_map<vector3d, int> &norms, std::unordered_map<vector3d, int> &verts, BSP_DefPoints &dpnts);
void MakeFlatPoly(BSP_FlatPoly &dst, pcs_polygon &src, std::unordered_map<vector3d, int> &norms, std::unordered_map<vector3d, int> &verts, BSP_DefPoints &dpnts);
void MakeDefPoints(BSP_DefPoints& dpnts, std::vector<bsp_vert> &pntslist);


//-----------------------------------------------
// functions for reading a BSP into PMF
//-----------------------------------------------

void BSPTransPMF(unsigned int offset, unsigned char *data,
	BSP_DefPoints &points, std::vector<pcs_polygon> &polygons,
	unsigned int &upolys);

void TranslateFPoly(unsigned int offset, unsigned char *data,
	BSP_DefPoints &points, std::vector<pcs_polygon> &polygons,
	unsigned int &upolys);

void TranslateTPoly(unsigned int offset, unsigned char *data,
	BSP_DefPoints &points, std::vector<pcs_polygon> &polygons,
	unsigned int &upolys);

void InterpretSortNorm(unsigned int offset, unsigned char *data,
	BSP_DefPoints &points, std::vector<pcs_polygon> &polygons,
	unsigned int &upolys);

//-----------------------------------------------
// functions for rendering a BSP's sortnorms and BBoxes
//-----------------------------------------------
void RenderBSP(unsigned int offset, unsigned char *data, vector3d obj_center);
void RenderSortnorm(unsigned int offset, unsigned char *data, vector3d obj_center);
void RenderBBox(unsigned int offset, unsigned char *data, vector3d obj_center);
void RenderUntextured(unsigned int offset, unsigned char *data, vector3d obj_center);
void RenderTextured(unsigned int offset, unsigned char *data, vector3d obj_center);
void OpenGL_RenderBox(vector3d min, vector3d max);


#endif // _pcs_pof_bspfuncs_h_

