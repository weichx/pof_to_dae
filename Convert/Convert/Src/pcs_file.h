#if !defined(_pcs_file_h_)
#define _pcs_file_h_


#include "vector3d.h"
#include "POFHandler.h"

#include "pcs_file_dstructs.h"
#include <string>
#include <vector>
#include <math.h>

#include "async_progress.h"

// forward dec
struct omnipoints;

#define PMF_VERSION 103
#define PMF_MIN_VERSION 100
#define PMF_MAX_VERSION PMF_VERSION

enum CHUNK {
	ROOT, HDR2, HDR2_SUB_DTL, HDR2_SUB_DEBRIS, ACEN, TXTR, SOBJ, PINF, EYE,
	SPCL, WEAP, GPNT, MPNT, WEAP_SUB, TGUN, TGUN_SUB,
	DOCK, FUEL, FUEL_SUB, SHLD, INSG, INSG_SUB,
	PATH, PATH_SUB, GLOW, GLOW_SUB
};

class PCS_Model
{
private:

	//********************************************************************
	// Data saved to disk
	//********************************************************************
	// ---- -------------------- ----
	// ----  Model Editing Data  ----
	// ---- -------------------- ----

	header_data header;

	vector3d autocentering; //ACEN

	// ---- -------------------- ----
	std::vector<std::string> textures; //TXTR
	std::vector<pcs_sobj> subobjects; // OBJ2
	std::vector<std::string> model_info; // PINF
	std::vector<pcs_eye_pos> eyes; // ' EYE'
	std::vector<pcs_special> special; // SPCL
	std::vector<pcs_slot> weapons; // GPNT and MPNT
	std::vector<pcs_turret> turrets; // TGUN and TMIS
	std::vector<pcs_dock_point> docking; // DOCK
	std::vector<pcs_thruster> thrusters; // FUEL
	std::vector<pcs_shield_triangle> shield_mesh; // SHLD
	std::vector<pcs_insig> insignia; // INSG
	std::vector<pcs_path> ai_paths; // PATH
	std::vector<pcs_glow_array> light_arrays; // GLOW
	// ---- -------------------- ----
	// ----       BSP Cache      ----
	// ---- -------------------- ----
	// #Caveat: BSP Cache/can_cache saved to PMF files starting with PMF version 102
	std::vector<pmf_bsp_cache> bsp_cache;
	bool can_bsp_cache;

	bool has_fullsmoothing_data;

	//********************************************************************
	// Data NOT saved to disk
	//********************************************************************

	// ---- -------------------- ----
	// ----    Rendering Data    ----
	// ---- -------------------- ----
	int active_submodel;
	int active_texture;
	bool Wireframe;
	bool Textureless;
	bool highlight_active_model;
	bool vbos_enabled;




	//// priv funcs
	//void RenderGeometryRecursive(int sobj, TextureControl &tc, bool use_vbos);
	//void RenderGeometry_vertex_buffers(int sobj, TextureControl &tc);
	int FindTexture(std::string name);


	bool PMFObj_to_POFObj2(int src_num, OBJ2 &dst, bool &bsp_compiled, float& model_radius);

	inline void POFTranslateBoundingBoxes(vector3d& min, vector3d& max) {
		float temp = -min.x;
		min.x = -max.x;
		max.x = temp;
	}

	bobboau::matrix moi_recalculate(float X, float Z, float xw, float zw);
	bobboau::matrix moi_recalculate(float X, float Z, float xw, float zw, float Y, float dy);
	bool moi_colide(std::vector<vector3d>&cpoints, float x, float z);

public:
	bobboau::matrix moi_recalculate(int Xres, int Yres);
	vector3d OffsetFromParent(int ObjNum);
	void Transform(const matrix& transform, const vector3d& translation);

	PCS_Model() : header(), can_bsp_cache(false), has_fullsmoothing_data(false), active_submodel(0), Wireframe(false), Textureless(false), vbos_enabled(false), draw_bsp(false)
	{

	}
	~PCS_Model(){
		for (unsigned int i = 0; i<subobjects.size(); i++){
			subobjects[i].destroy_vertex_buffer();
		}
	}

	// Renderer commands
	void Rcall_Wireframe(bool tf) { Wireframe = tf; }
	void Rcall_Textureless(bool tf) { Textureless = tf; }

	bool draw_bsp;
	//

	//statistic trackers
	static unsigned long BSP_TREE_TIME;
	static unsigned int BSP_MAX_DEPTH;
	static unsigned int BSP_CUR_DEPTH;
	static unsigned int BSP_NODE_POLYS;
	static bool BSP_COMPILE_ERROR;


	void Reset();
	void PurgeBSPcache() { bsp_cache.resize(0); can_bsp_cache = false; }
	// Loaders
	int LoadFromPMF(std::string filename, AsyncProgress* progress); //PMF = PCS Model File
	int LoadFromPOF(std::string filename, AsyncProgress* progress);

	int LoadFromCOB(std::string filename, AsyncProgress* progress, float scaler, bool useFilter);
	int LoadFromDAE(std::string filename, AsyncProgress* progress, bool mirror_x, bool mirror_y, bool mirror_z);


	//geometery manipulation

	//splits the poly at position I into two seperate polyogns, along the ith and jth vert
	static bool split_poly(std::vector<pcs_polygon>&polys, int I, int i, int j);

	//poly fixing

	//fixes a convex polygon within a list
	//also fixes polygons whith more than 20 points
	static void filter_polygon(std::vector<pcs_polygon>&polys, int i);

	//calls the above on all polys in a list
	static void filter_geometry(std::vector<pcs_polygon>&polys);


	// brute-force recovers the smoothing angle numbers for each polygon
	void Calculate_Smoothing_Data(int &sobjs_comp); // will IMMEDIATELY TERMINATE if has_fullsmoothing_data == true

	// Savers
	int SaveToPMF(std::string filename, AsyncProgress* progress); // PMF = PCS Model File
	int SaveToPOF(std::string filename, AsyncProgress* progress);

	int SaveToCOB(std::string filename, AsyncProgress* progress, float scaler);

	// These functions don't actually work, and they're not likely to - just kept here for historical purposes
	// the ASE format was found to be unreliable
	//int SaveToASE(std::string filename, AsyncProgress* progress, float scaler);
	//int LoadFromASE(std::string filename, AsyncProgress* progress, float scaler);
	int SaveToDAE(std::string filename, AsyncProgress* progress, int helpers, int props_as_helpers);


	// Rendering
	/*void Render(TextureControl &tc, bool use_vbos, bool highlight = false);
	void draw_shields();
	void draw_insignia(int lod, const omnipoints& omni);*/


	// Accessors
	const header_data&get_header(){ return header; }
	void set_header(const header_data&h){ header = h; }
	float GetMaxRadius();
	float GetMass();
	void GetMOI(std::vector<float>& tensor);

	vector3d GetMinBounding();
	vector3d GetMaxBounding();
	vector3d GetCenterOfMass();
	vector3d GetAutoCenter();


	size_t GetModelInfoCount() { return model_info.size(); }
	int GetLODCount();
	int GetDebrisCount();
	int GetCrossSectCount();
	int GetTexturesCount();
	int GetSOBJCount();
	int GetEyeCount();
	int GetSpecialCount();
	int GetWeaponCount();
	int GetTurretCount();
	int GetDockingCount();
	int GetThrusterCount();
	int GetShldTriCount();
	int GetInsigniaCount();
	int GetPathCount();
	int GetLightCount();

	// Referencers (both Acc/Mod)

	int&					LOD(unsigned int idx);
	int&					Debris(unsigned int idx);
	pcs_crs_sect&			CrossSect(unsigned int idx);
	std::string&				Texture(unsigned int idx);
	pcs_sobj&				SOBJ(unsigned int idx);
	pcs_eye_pos&			Eye(unsigned int idx);
	pcs_special&			Special(unsigned int idx);
	pcs_slot&				Weapon(unsigned int idx);
	pcs_turret&				Turret(unsigned int idx);
	pcs_dock_point&			Dock(unsigned int idx);
	pcs_thruster&			Thruster(unsigned int idx);
	pcs_shield_triangle&	ShldTri(unsigned int idx);
	pcs_insig&				Insignia(unsigned int idx);
	pcs_path&				Path(unsigned int idx);
	pcs_glow_array&			Light(unsigned int idx);
	std::string&				ModelInfo(unsigned int idx);

	// Modifiers
	void SetMaxRadius(float rad);
	void SetMass(float mass);
	void SetMOI(float tensor[3][3]); // float[3][3]

	void SetMinBounding(const vector3d &bnd);
	void SetMaxBounding(const vector3d &bnd);
	void SetCenterOfMass(const vector3d &cnt);
	void SetAutoCenter(const vector3d &cnt);

	void AddModelInfo(std::string info = "");
	void SetModelInfo(unsigned int idx, std::string info);

	void SetNumLODs(int num) { header.detail_levels.resize(num); }
	void AddLOD(int sobj = -1); // we can add an emtpy lod, or initialize it upon creation
	void DelLOD(unsigned int idx);

	void SetNumDebris(int num) { header.debris_pieces.resize(num); }
	void AddDebris(int sobj = -1); // we can add an emtpy lod, or initialize it upon creation
	void DelDebris(unsigned int idx);

	void SetNumCrossSects(int num) { header.cross_sections.resize(num); }
	void AddCrossSect(pcs_crs_sect *cs = NULL); // we can add an emtpy lod, or initialize it upon creation
	void DelCrossSect(unsigned int idx);

	int maybeAddTexture(std::string txt);
	void AddTextures(std::string txt = "");
	void DelTextures(unsigned int idx);

	void AddSOBJ(pcs_sobj *obj = NULL);
	void DelSOBJ(int idx);
	void SetObjectChanged(unsigned int idx);

	void AddEye(pcs_eye_pos *eye = NULL);
	void DelEye(unsigned int idx);

	void AddSpecial(pcs_special *spcl = NULL);
	void DelSpecial(unsigned int idx);

	void AddWeapon(pcs_slot *weap = NULL);
	void DelWeapon(unsigned int idx);

	void AddTurret(pcs_turret *trrt = NULL);
	void DelTurret(unsigned int idx);

	void AddDocking(pcs_dock_point *dock = NULL);
	void DelDocking(unsigned int idx);

	void AddThruster(pcs_thruster *thrust = NULL);
	void DelThruster(unsigned int idx);

	void AddShldTri(pcs_shield_triangle *stri = NULL);
	void DelShldTri(unsigned int idx);

	void AddInsignia(pcs_insig *insig = NULL);
	void DelInsignia(unsigned int idx);

	void AddPath(pcs_path *path = NULL);
	void DelPath(unsigned int idx);

	void AddLight(pcs_glow_array *lights = NULL);
	void DelLight(unsigned int idx);

	//interface accessors
	std::vector<std::string> get_textures(){ return textures; }
	void set_textures(const std::vector<std::string> &t){
		if (t.size() != textures.size()) {
			textures = t;
			make_vertex_buffers();
		}
		else {
			textures = t;
		}
	}

	std::vector<pcs_sobj> get_subobjects(){ return subobjects; }
	void set_subobjects(const std::vector<pcs_sobj> &t){ subobjects = t; }

	std::vector<std::string> get_model_info(){ return model_info; }
	void set_model_info(const std::vector<std::string> &t){ model_info = t; }

	std::vector<pcs_eye_pos> get_eyes(){ return eyes; }
	void set_eyes(const std::vector<pcs_eye_pos> &t){ eyes = t; }

	std::vector<pcs_special> get_special(){ return special; }
	void set_special(const std::vector<pcs_special> &t){ special = t; }

	std::vector<pcs_slot> get_weapons(){ return weapons; }
	void set_weapons(const std::vector<pcs_slot> &t){ weapons = t; }

	std::vector<pcs_turret> get_turrets(){ return turrets; }
	void set_turrets(const std::vector<pcs_turret> &t){ turrets = t; }

	std::vector<pcs_dock_point> get_docking(){ return docking; }
	void set_docking(const std::vector<pcs_dock_point> &t){ docking = t; }

	std::vector<pcs_thruster> get_thrusters(){ return thrusters; }
	void set_thrusters(const std::vector<pcs_thruster> &t){ thrusters = t; }

	std::vector<pcs_shield_triangle> get_shield_mesh(){ return shield_mesh; }
	void set_shield_mesh(const std::vector<pcs_shield_triangle> &t){ shield_mesh = t; }

	std::vector<pcs_insig> get_insignia(){ return insignia; }
	void set_insignia(const std::vector<pcs_insig> &t){ insignia = t; }

	std::vector<pcs_path> get_ai_paths(){ return ai_paths; }
	void set_ai_paths(const std::vector<pcs_path> &t){ ai_paths = t; }

	std::vector<pcs_glow_array> get_glow_points(){ return light_arrays; }
	void set_glow_points(const std::vector<pcs_glow_array> &t){ light_arrays = t; }

	void set_active_model(int idx){ active_submodel = idx; }
	int get_active_model(){ return active_submodel; };

	void set_active_texture(int idx){ active_texture = idx; }
	int get_active_texture(){ return active_texture; };

	int find_LOD_root(int idx){
		if (idx<0)
			return 0;
		if (subobjects[idx].parent_sobj <0)
			return idx;
		for (unsigned int i = 0; i<header.detail_levels.size(); i++){
			if (idx == header.detail_levels[i])
				return idx;
		}
		return find_LOD_root(subobjects[idx].parent_sobj);
	}

	bool is_debris(int idx){
		if (idx<0)
			return false;
		for (unsigned int i = 0; i<header.debris_pieces.size(); i++){
			if (idx == header.debris_pieces[i])
				return true;
		}
		return false;
	}

	vector3d get_model_offset(int i){
		if (i >= (int)subobjects.size())
			return vector3d(0, 0, 0);
		if (i <0)
			return autocentering;
		return get_model_offset(subobjects[i].parent_sobj) + subobjects[i].offset;
	}

	//gets the poly count of all children of subobject idx
	size_t get_child_subobj_poly_count(int idx){
		size_t total = 0;
		//we realy need to make a proper tree
		for (unsigned int i = 0; i<subobjects.size(); i++){
			if (subobjects[i].parent_sobj == idx){
				total += subobjects[i].polygons.size() + get_child_subobj_poly_count(i);
			}
		}
		return total;
	}


	bool get_bsp_cache_status(){ return can_bsp_cache; }

	//gets the average size of the model
	float get_avg_dimintion(){
		float d[6];
		d[0] = fabs(header.min_bounding.x);
		d[1] = fabs(header.min_bounding.y);
		d[2] = fabs(header.min_bounding.z);
		d[3] = fabs(header.max_bounding.x);
		d[4] = fabs(header.max_bounding.y);
		d[5] = fabs(header.max_bounding.z);

		float avg = 0;
		for (int i = 0; i<6; i++){
			avg += d[i];
		}
		avg /= 6.0f;
		return avg;
	}

	void init_vertex_buffers(bool enabled);//sets up all vertex buffers for this model
	void make_vertex_buffers();//sets up all vertex buffers for this model
	void make_vertex_buffer(int sobj);//sets up vertex buffers for a subobject

};


#endif //_pcs_file_h_

