#if defined(_WIN32)
#include <windows.h>
#endif


#include "async_progress.h"
#include "pcs_file.h"
#include "pcs_pof_bspfuncs.h"
#include <fstream>
#include <set>
#include <cfloat>
#include "color.h"
#include "omnipoints.h"
#include <iostream>

#include <fstream>
#include <cfloat>
#include "pcs_file.h"
#include "POFHandler.h"
#include "BSPHandler.h"
#include "pcs_pof_bspfuncs.h"
#include "DAEHandler.h"
#include <boost/scoped_ptr.hpp>
#include <iostream>

//#include "pcs2.h"


#define M_PI 3.14159265359
unsigned int PCS_Model::BSP_MAX_DEPTH = 0;
unsigned int PCS_Model::BSP_CUR_DEPTH = 0;
unsigned int PCS_Model::BSP_NODE_POLYS = 1;
bool PCS_Model::BSP_COMPILE_ERROR = false;
unsigned long PCS_Model::BSP_TREE_TIME = 0; //was wxLongLong

//splits the poly at position I into two seperate polyogns, along the ith and jth vert
bool PCS_Model::split_poly(std::vector<pcs_polygon>&polys, int I, int i, int j){

	if (i>j){
		//for simplicities sake i will be less than j
		int temp = i;
		i = j;
		j = temp;
	}

	if (polys[I].verts.size() < 4 ||
		//can't split a triangle
		i == j ||
		//can't split along one vert
		j - i<2 || (i == 0 && (unsigned)j == polys[I].verts.size() - 1)){
		//one of them would have fewer than 3 verts: degenerate
		return false;
	}


	polys.push_back(polys[I]);
	//make a copy of the poly
	pcs_polygon&split1 = polys[I];
	pcs_polygon&split2 = polys[polys.size() - 1];

	//then remove the uneeded verts from the first and second
	int h;
	//remove the first i verts and everything after j
	for (h = 0; h<i; h++)
		split1.verts.erase(split1.verts.begin());
	split1.verts.resize(j - i + 1);

	//remove everything between i and j (non inclusive) from the second
	for (h = i + 1; h<j; h++)
		split2.verts.erase(split2.verts.begin() + (i + 1));

	//fix normal
	vector3d norm(0, 0, 0);

	for (h = 0; h<(int)split1.verts.size(); h++)
		norm = norm + CrossProduct(split1.verts[(h + split1.verts.size() - 1) % split1.verts.size()].point - split1.verts[h].point, split1.verts[(h + 1) % split1.verts.size()].point - split1.verts[h].point);
	split1.norm = MakeUnitVector(norm);// / float(split1.verts.size());

	norm = vector3d(0, 0, 0);

	for (h = 0; h<(int)split2.verts.size(); h++)
		norm = norm + CrossProduct(split2.verts[(h + split2.verts.size() - 1) % split2.verts.size()].point - split2.verts[h].point, split2.verts[(h + 1) % split2.verts.size()].point - split2.verts[h].point);
	split2.norm = MakeUnitVector(norm);// / float(split2.verts.size());


	return true;
}

//return true if closest point is outside either line segments
bool closest_line_pt(vector3d l1p1, vector3d l1p2, vector3d l2p1, vector3d l2p2, vector3d*closest1 = NULL, vector3d*closest2 = NULL){

	vector3d ln1 = MakeUnitVector(l1p2 - l1p1);
	vector3d ln2 = MakeUnitVector(l2p2 - l2p1);

	vector3d In = CrossProduct(ln1, ln2);

	vector3d pnA = MakeUnitVector(CrossProduct(ln1, In));
	vector3d pnB = MakeUnitVector(CrossProduct(ln2, In));

	float d;

	d = dot(pnB, ln1);
	if (d == 0.0f)return true;
	float tA = -dot(pnB, l1p1 - l2p1) / d;

	d = dot(pnA, ln2);
	if (d == 0.0f)return true;
	float tB = -dot(pnA, l2p1 - l1p1) / d;

	if (closest1){
		(*closest1) = l1p1 + ln1*tA;
	}
	if (closest2){
		(*closest2) = l2p1 + ln2*tB;
	}
	if (tA<0.0f || tB<0.0f || Distance(l1p1, l1p2) < tA || Distance(l2p1, l2p2) < tB)
		return true;
	else
		return false;
}

float full_angle(vector3d A, vector3d B, vector3d C, vector3d pln){
	A = MakeUnitVector(A - C);
	B = MakeUnitVector(B - C);

	float ang = acos(dot(A, B));

	if (dot(CrossProduct(B, A), pln) > 0)
		ang = M_PI*2.0f - ang;

	return ang;
}

void interconect_poly_on_verts(std::vector<pcs_polygon>&polys, int i, std::vector<unsigned int>&concaves, vector3d&pnorm){

	pcs_polygon&poly = polys[i];
	const unsigned int psize = polys[i].verts.size();

	if (concaves.size() == 1){
		//in this case there is no other concave point
		//so let's split it in half and try that
		//I'm thinking in most cases spliting in half will be the best result

		//but it would be better if I tried every posable split untill I found one 
		//that had at least one good polygon with the largest number of verts in it 
		//less than half the total

		//but this should be good enough

		//probly

		if (PCS_Model::split_poly(polys, i, concaves[0], (concaves[0] + poly.verts.size() / 2) % poly.verts.size())){
			PCS_Model::filter_polygon(polys, i);
			//the first split poly will be in the same place as the one we passed
		}
	}
	else if (concaves.size() > 1){
		//we have more than one concave point, try splitting between them to see if
		//we can find a split that will not cause problems
		unsigned int s;
		unsigned int t;
		for (s = 0; s < concaves.size(); s++){
			vector3d&sth = poly.verts[concaves[s]].point;
			vector3d&sthm1 = poly.verts[(concaves[s] - 1 + poly.verts.size()) % poly.verts.size()].point;
			vector3d&sthp1 = poly.verts[(concaves[s] + 1) % poly.verts.size()].point;
			for (t = s + 1; t < concaves.size(); t++){
				if ((concaves[s] + 1) % poly.verts.size() == concaves[t] || (concaves[t] + 1) % poly.verts.size() == concaves[s])
					continue;

				vector3d&tth = poly.verts[concaves[t]].point;

				//make sure we don't go into the bad place
				float  ang = full_angle(sthm1, sthp1, sth, pnorm);
				float tang = full_angle(sthm1, tth, sth, pnorm);
				if (tang <= ang)
					continue;
				//if the angle between s and t is between the angle for s and s+1, 
				//then we can dismis this pair right now

				//for every concave pair
				//test to see if the split made by it will cause problems
				//by seeing if there is an intersection with any other line
				unsigned int j;
				for (j = 0; j < poly.verts.size(); j++){
					vector3d test1, test2;
					if (concaves[s] == j || concaves[s] == (j + 1) % psize || concaves[t] == j || concaves[t] == (j + 1) % psize)
						continue;
					if (!closest_line_pt(poly.verts[j].point, poly.verts[(j + 1) % poly.verts.size()].point, sth, tth, &test1, &test2))
						break;
				}
				if (j == poly.verts.size())break;
			}
			if (t < concaves.size())break;
		}

		//if we found a good pair make the split, 
		if (s < concaves.size()){
			if (PCS_Model::split_poly(polys, i, concaves[s], concaves[t])){
				PCS_Model::filter_polygon(polys, i);
				//the first split poly will be in the same place as the one we passed
			}
		}
		else{
			//otherwise do an exaustive search for any split that will make a good poly
			//that starts with a concave point
			for (s = 0; s < concaves.size(); s++){
				vector3d&sth = poly.verts[concaves[s]].point;
				vector3d&sthm1 = poly.verts[(concaves[s] - 1 + poly.verts.size()) % poly.verts.size()].point;
				vector3d&sthp1 = poly.verts[(concaves[s] + 1) % poly.verts.size()].point;
				unsigned int T;
				for (T = 0; T < poly.verts.size() - 1; T++){
					t = (T + concaves[s] + 2) % poly.verts.size();

					if ((concaves[s] + 1) % poly.verts.size() == t || (t + 1) % poly.verts.size() == concaves[s])
						continue;

					vector3d&tth = poly.verts[t].point;

					//make sure we don't go into the bad place
					float  ang = full_angle(sthm1, sthp1, sth, pnorm);
					float tang = full_angle(sthm1, tth, sth, pnorm);
					if (tang <= ang)
						continue;
					//if the angle between s and t is between the angle for s and s+1, 
					//then we can dismis this pair right now

					//for every pair
					//test to see if the split made by it will cause problems
					//by seeing if there is an intersection with any other line
					unsigned int j;
					for (j = 0; j < poly.verts.size(); j++){
						if (concaves[s] == j || concaves[s] == (j + 1) % psize || t == j || t == (j + 1) % psize)
							continue;
						vector3d&jth = poly.verts[j].point;
						//vector3d&jthm1 = poly.verts[(j-1+poly.verts.size())%poly.verts.size()].point;
						vector3d&jthp1 = poly.verts[(j + 1) % poly.verts.size()].point;
						if (!closest_line_pt(jth, jthp1, sth, tth))
							break;
					}
					//if it made it to the end then s and t are good
					if (j == poly.verts.size())break;
				}
				//if t didn't make it to the end of the list, then we found a good pair in the t loop
				if (T < poly.verts.size() - 1)break;
			}
			if (s < concaves.size()){
				//we found a good split so make it
				if (PCS_Model::split_poly(polys, i, concaves[s], t)){
					PCS_Model::filter_polygon(polys, i);
					//the first split poly will be in the same place as the one we passed
				}
			}
			else{
				//now we are fairly screwed, as there is no way to fix this
				std::cout << "*ERROR*:uncorectable geometry encountered! \n\nTruly this is the darkest of hours.";
			}
		}
	}
}

//fixes concave polygons
void PCS_Model::filter_polygon(std::vector<pcs_polygon>&polys, int i){

	if (polys[i].verts.size() < 4)return;
	//no need

	pcs_polygon&poly = polys[i];

	unsigned int psize = polys[i].verts.size();

	vector3d avg(0, 0, 0);
	for (unsigned int j = 0; j<psize; j++){
		avg += poly.verts[j].point;
	}
	avg = avg / psize;

	std::vector<vector3d> norms(psize);
	vector3d pnorm;
	for (unsigned int j = 0; j<psize; j++){
		norms[j] = CrossProduct(
			poly.verts[(j + 1) % poly.verts.size()].point - poly.verts[j].point,
			poly.verts[(j - 1 + poly.verts.size()) % poly.verts.size()].point - poly.verts[j].point);
		if (!no_nan(norms[j])){
			norms[j] = vector3d(0, 0, 0);
		}
		pnorm += norms[j];
		norms[j] = MakeUnitVector(norms[j]);
	}
	pnorm = MakeUnitVector(pnorm);



	//go through and try to find 2 or more concave points

	std::vector<unsigned int> concaves;

	//make up a list of concave points
	for (unsigned int j = 0; j<psize; j++){
		if (dot(norms[j], pnorm) <= 0.0f){
			//this vert is concave
			concaves.push_back(j);
		}
	}

	interconect_poly_on_verts(polys, i, concaves, pnorm);

	//poly should be convex at this point, more or less
	//now lets go for coplanar


	psize = polys[i].verts.size();

	norms.resize(psize);
	//these need to be rebuilt, because if the above function did anything they are totaly invalid
	for (unsigned int j = 0; j<psize; j++){
		norms[j] = CrossProduct(
			polys[i].verts[(j + 1) % polys[i].verts.size()].point - polys[i].verts[j].point,
			polys[i].verts[(j - 1 + polys[i].verts.size()) % polys[i].verts.size()].point - polys[i].verts[j].point);
		if (!no_nan(norms[j])){
			norms[j] = vector3d(0, 0, 0);
		}
		pnorm += norms[j];
		norms[j] = MakeUnitVector(norms[j]);
	}

	pnorm = MakeUnitVector(pnorm);
	std::vector<unsigned int>&nonplanar = concaves;//lets reuse this
	nonplanar.resize(0);

	//make up a list of nonplanar points
	for (unsigned int j = 0; j<psize; j++){
		if (dot(norms[j], pnorm) <= 0.999f){
			//this vert is concave
			nonplanar.push_back(j);
		}
	}

	interconect_poly_on_verts(polys, i, nonplanar, pnorm);


	if (polys[i].verts.size() > 20){
		//fix polys with too many verts
		if (split_poly(polys, i, 0, poly.verts.size() / 2)){
			filter_polygon(polys, i);
			//the first split poly will be in the same place as the one we passed
		}
		return;
	}
}

void PCS_Model::filter_geometry(std::vector<pcs_polygon>&polys){
	for (unsigned int i = 0; i< polys.size(); i++){
		filter_polygon(polys, i);
	}
}

// PP shortcuts
#define PCS_ADD_TO_VEC(vec, var) unsigned int idx = vec.size(); \
								 vec.resize(idx+1); \
								 if (var) \
									vec[idx] = *var;

//****************************************************************************
// PCS Model File (PMF) Loader
//****************************************************************************

/*
File Version History
100: Initial
101: Added "facet_angle" to polygon
103: Save correct BSP cache
*/
int PCS_Model::LoadFromPMF(std::string filename, AsyncProgress* progress)
//PMF = PCS Model File
{
	this->Reset();
	progress->setTarget(17);
	std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
	if (!in)
		return 1;

	progress->incrementWithMessage("Checking Filesig");
	// -------- header -------- 
	char sig[] = { 'P', 'M', 'F', '1' }, fs[4];
	int ver;

	in.read(fs, 4);
	if (strncmp(fs, sig, 4))
		return 2; // filesig failed
	BFRead(ver, int)

		if (ver < PMF_MIN_VERSION || ver > PMF_MAX_VERSION)
			return 3; // version failed
	// --------  data  -------- 

	progress->incrementWithMessage("Reading Header");

	BFRead(header.max_radius, float)
		BFRead(header.min_bounding, vector3d)
		BFRead(header.max_bounding, vector3d)
		header.max_radius_override = header.max_radius;
	header.min_bounding_override = header.min_bounding;
	header.max_bounding_override = header.max_bounding;
	BFReadVector(header.detail_levels)
		BFReadVector(header.debris_pieces)
		BFRead(header.mass, float)
		BFRead(header.mass_center, vector3d)
		in.read((char*)header.MOI, sizeof(float) * 9);
	BFReadVector(header.cross_sections)
		BFRead(autocentering, vector3d)

		progress->incrementWithMessage("Reading Textures");

	unsigned int i;
	BFRead(i, int)
		textures.resize(i);
	for (i = 0; i < textures.size(); i++) {
		BFReadString(textures[i])
			cout << textures[i] << endl;
	}

		progress->incrementWithMessage("Reading Subobjects");
	BFReadAdvVector(subobjects)

		progress->incrementWithMessage("Reading Info Strings");
	BFRead(i, int)
		model_info.resize(i);
	for (i = 0; i < model_info.size(); i++)
		BFReadString(model_info[i])

		progress->incrementWithMessage("Reading Eyes");
	BFReadVector(eyes)

		progress->incrementWithMessage("Reading Specials");
	BFReadAdvVector(special)

		progress->incrementWithMessage("Reading Weapons");
	BFReadAdvVector(weapons)

		progress->incrementWithMessage("Reading Turrets");
	BFReadAdvVector(turrets)

		progress->incrementWithMessage("Reading Docking");
	BFReadAdvVector(docking)

		progress->incrementWithMessage("Reading Thrusters");
	BFReadAdvVector(thrusters)

		progress->incrementWithMessage("Reading Shields");
	BFReadVector(shield_mesh)

		progress->incrementWithMessage("Reading Insignia");
	BFReadAdvVector(insignia)

		progress->incrementWithMessage("Reading Paths");
	BFReadAdvVector(ai_paths)

		progress->incrementWithMessage("Reading Glows");
	BFReadAdvVector(light_arrays)

		progress->incrementWithMessage("Reading BSP Cache");
	if (ver >= 102)
	{
		BFReadAdvVector(bsp_cache)
			BFRead(can_bsp_cache, bool)
			if (ver == 102) {
				// XXX: bsp_data will always be junk so never use cached data.
				can_bsp_cache = false;
			}

		// new in ver 102.. but not associated with above
		BFRead(has_fullsmoothing_data, bool)

	}
	Transform(matrix(), vector3d());
	header.max_radius_overridden = header.max_radius_override != header.max_radius;
	header.min_bounding_overridden = header.min_bounding_override != header.min_bounding;
	header.max_bounding_overridden = header.max_bounding_override != header.max_bounding;
	for (size_t i = 0; i < subobjects.size(); ++i) {
		pcs_sobj& sobj = subobjects[i];
		sobj.radius_overridden = fabs(sobj.radius - sobj.radius_override) > 0.0001f;
		sobj.bounding_box_min_point_overridden = sobj.bounding_box_min_point != sobj.bounding_box_min_point_override;
		sobj.bounding_box_max_point_overridden = sobj.bounding_box_max_point != sobj.bounding_box_max_point_override;
	}

	progress->incrementProgress();

	return 0;
}

//****************************************************************************
// PCS Model File (PMF) Saver
//****************************************************************************

int PCS_Model::SaveToPMF(std::string filename, AsyncProgress* progress)
// PMF = PCS Model File
{
	progress->setTarget(17);

	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if (!out)
		return 1;

	progress->incrementWithMessage("Writing Filesig");
	// -------- header -------- 
	char sig[] = { 'P', 'M', 'F', '1' };
	int ver = PMF_VERSION;

	out.write(sig, 4);
	BFWrite(ver, int)

		progress->incrementWithMessage("Writing Header");
	// --------  data  -------- 

	BFWrite(header.max_radius_overridden ? header.max_radius_override : header.max_radius, float);
	BFWrite(header.min_bounding_overridden ? header.min_bounding_override : header.min_bounding, vector3d);
	BFWrite(header.max_bounding_overridden ? header.max_bounding_override : header.max_bounding, vector3d);
	BFWriteVector(header.detail_levels)
		BFWriteVector(header.debris_pieces)
		BFWrite(header.mass, float)
		BFWrite(header.mass_center, vector3d)
		out.write((char*)header.MOI, sizeof(float) * 9);
	BFWriteVector(header.cross_sections)
		BFWrite(autocentering, vector3d)


		progress->incrementWithMessage("Writing Textures");
	unsigned int i = textures.size();
	BFWrite(i, int)
		for (i = 0; i < textures.size(); i++)
			BFWriteString(textures[i])


			progress->incrementWithMessage("Writing Subobjects");
	BFWriteAdvVector(subobjects)

		progress->incrementWithMessage("Writing Info Strings");
	i = model_info.size();
	BFWrite(i, int)
		for (i = 0; i < model_info.size(); i++)
			BFWriteString(model_info[i])

			progress->incrementWithMessage("Writing Eyes");
	BFWriteVector(eyes)

		progress->incrementWithMessage("Writing Special");
	BFWriteAdvVector(special)

		progress->incrementWithMessage("Writing Weapons");
	BFWriteAdvVector(weapons)

		progress->incrementWithMessage("Writing Turrets");
	BFWriteAdvVector(turrets)

		progress->incrementWithMessage("Writing Docking");
	BFWriteAdvVector(docking)

		progress->incrementWithMessage("Writing Thrusters");
	BFWriteAdvVector(thrusters)

		progress->incrementWithMessage("Writing Shields");
	BFWriteVector(shield_mesh)

		progress->incrementWithMessage("Writing Insignia");
	BFWriteAdvVector(insignia)

		progress->incrementWithMessage("Writing Paths");
	BFWriteAdvVector(ai_paths)

		progress->incrementWithMessage("Writing Glows");
	BFWriteAdvVector(light_arrays)

		progress->incrementWithMessage("Writing BSP Cache");
	BFWriteAdvVector(bsp_cache)
		BFWrite(can_bsp_cache, bool)


		BFWrite(has_fullsmoothing_data, bool)

		progress->incrementProgress();
	return 0;
}

//****************************************************************************
// Resetter
//****************************************************************************

void PCS_Model::Reset()
{
	//leaveing this here in case other code uses it, it needs to be something
	active_submodel = 0;
	active_texture = -1;
	highlight_active_model = false;
	//Wireframe = false;
	//Textureless = false;
	can_bsp_cache = false;
	has_fullsmoothing_data = false;
	vbos_enabled = false;


	header.max_radius = 0;

	autocentering = header.mass_center = header.min_bounding = header.max_bounding = vector3d(0, 0, 0);

	header.detail_levels.resize(0);
	header.debris_pieces.resize(0);
	header.mass = 0.0;

	memset(header.MOI, 0, sizeof(float) * 9);

	//for (size_t i = 0; i < subobjects.size(); i++)
    // subobjects[i].destroy_vertex_buffer();

	header.cross_sections.resize(0);
	textures.resize(0);
	subobjects.resize(0);
	model_info.resize(0);
	eyes.resize(0);
	special.resize(0);
	weapons.resize(0);
	turrets.resize(0);
	docking.resize(0);
	thrusters.resize(0);
	shield_mesh.resize(0);
	insignia.resize(0);
	ai_paths.resize(0);
	light_arrays.resize(0);

	bsp_cache.resize(0);

}

//****************************************************************************
// Accessors
//****************************************************************************

float PCS_Model::GetMaxRadius()
{
	return header.max_radius;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

float PCS_Model::GetMass()
{
	return header.mass;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::GetMOI(std::vector<float>& tensor)
// float[3][3]
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			tensor[(3 * i) + j] = header.MOI[i][j];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


vector3d PCS_Model::GetMinBounding()
{
	return header.min_bounding;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

vector3d PCS_Model::GetMaxBounding()
{
	return header.max_bounding;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

vector3d PCS_Model::GetCenterOfMass()
{
	return header.mass_center;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

vector3d PCS_Model::GetAutoCenter()
{
	return autocentering;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


std::string& PCS_Model::ModelInfo(unsigned int idx)
{
	return model_info[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-



int PCS_Model::GetLODCount()
{
	return header.detail_levels.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetDebrisCount()
{
	return header.debris_pieces.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetCrossSectCount()
{
	return header.cross_sections.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetTexturesCount()
{
	return textures.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetSOBJCount()
{
	return subobjects.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetEyeCount()
{
	return eyes.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetSpecialCount()
{
	return special.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetWeaponCount()
{
	return weapons.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetTurretCount()
{
	return turrets.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetDockingCount()
{
	return docking.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetThrusterCount()
{
	return thrusters.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetShldTriCount()
{
	return shield_mesh.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetInsigniaCount()
{
	return insignia.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetPathCount()
{
	return ai_paths.size();
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::GetLightCount()
{
	return light_arrays.size();
}


//****************************************************************************
// Referencers (both Acc/Mod)
//****************************************************************************

int&					PCS_Model::LOD(unsigned int idx)
{
	return header.detail_levels[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int&					PCS_Model::Debris(unsigned int idx)
{
	return header.debris_pieces[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_crs_sect&			PCS_Model::CrossSect(unsigned int idx)
{
	return header.cross_sections[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

std::string&			PCS_Model::Texture(unsigned int idx)
{
	return textures[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_sobj&				PCS_Model::SOBJ(unsigned int idx)
{
	return subobjects[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_eye_pos&			PCS_Model::Eye(unsigned int idx)
{
	return eyes[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_special&			PCS_Model::Special(unsigned int idx)
{
	return special[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_slot&				PCS_Model::Weapon(unsigned int idx)
{
	return weapons[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_turret&				PCS_Model::Turret(unsigned int idx)
{
	return turrets[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_dock_point&			PCS_Model::Dock(unsigned int idx)
{
	return docking[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_thruster&			PCS_Model::Thruster(unsigned int idx)
{
	return thrusters[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_shield_triangle&	PCS_Model::ShldTri(unsigned int idx)
{
	return shield_mesh[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_insig&				PCS_Model::Insignia(unsigned int idx)
{
	return insignia[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_path&				PCS_Model::Path(unsigned int idx)
{
	return ai_paths[idx];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

pcs_glow_array&			PCS_Model::Light(unsigned int idx)
{
	return light_arrays[idx];
}


//****************************************************************************
// Modifiers
//****************************************************************************

void PCS_Model::SetMaxRadius(float rad)
{
	header.max_radius = rad;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::SetMass(float mass)
{
	this->header.mass = mass;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::SetMOI(float tensor[3][3])
// float[3][3]
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			header.MOI[i][j] = tensor[i][j];
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::SetMinBounding(const vector3d &bnd)
{
	header.min_bounding = bnd;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::SetMaxBounding(const vector3d &bnd)
{
	header.max_bounding = bnd;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::SetCenterOfMass(const vector3d &cnt)
{
	header.mass_center = cnt;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::SetAutoCenter(const vector3d &cnt)
{
	autocentering = cnt;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::AddModelInfo(std::string info)
{
	/*unsigned int idx = model_info.size();
	model_info.resize(idx+1);
	model_info[idx] = info;*/
	model_info.push_back(info);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::SetModelInfo(unsigned int idx, std::string info)
{
	model_info[idx] = info;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddLOD(int sobj)
{
	/*	unsigned int idx = header.detail_levels.size();
	header.detail_levels.resize(idx+1);
	header.detail_levels[idx] = sobj;*/
	header.detail_levels.push_back(sobj);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelLOD(unsigned int idx)
{
	RemoveIndex(header.detail_levels, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::AddDebris(int sobj)
{
	/*	unsigned int idx = header.debris_pieces.size();
	header.debris_pieces.resize(idx+1);
	header.debris_pieces[idx] = sobj;*/
	header.debris_pieces.push_back(sobj);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelDebris(unsigned int idx)
{
	RemoveIndex(header.debris_pieces, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddCrossSect(pcs_crs_sect *cs)
{
	PCS_ADD_TO_VEC(header.cross_sections, cs)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelCrossSect(unsigned int idx)
{
	RemoveIndex(header.cross_sections, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::maybeAddTexture(std::string txt)
{
	int index = FindInList<std::string>(textures, txt);
	if (index != -1)
		return index;
	index = textures.size();
	textures.push_back(txt);
	return index;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::AddTextures(std::string txt)
{
	//unsigned int idx = textures.size();
	//textures.resize(idx+1);
	//textures[idx] = txt;
	textures.push_back(txt);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelTextures(unsigned int idx)
{
	RemoveIndex(textures, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddSOBJ(pcs_sobj *obj)
{
	pmf_bsp_cache cache;
	cache.decache();
	bsp_cache.push_back(cache);
	if (obj)
		subobjects.push_back(*obj);
	else
	{
		pcs_sobj empty;
		subobjects.push_back(empty);
	}
	/*if (vbos_enabled) {
		subobjects.back().vertex_buffer.clear();
		subobjects.back().line_vertex_buffer.clear();
		make_vertex_buffer(subobjects.size() - 1);
	}*/
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelSOBJ(int idx)
{
	//we have to push down every reference to every subobject after this one
	//and mark each reference to this one as -1 (none, or root) or just erase it, or reparent it

	//start with subobject parents
	unsigned int i;
	for (i = 0; i<subobjects.size(); i++){
		if (subobjects[i].parent_sobj == idx)
			subobjects[i].parent_sobj = subobjects[idx].parent_sobj;//reparent to my parent
		if (subobjects[i].parent_sobj >(int)idx)
			subobjects[i].parent_sobj--;//it's actualy been moved one back, or it will be at the end of this function anyway
	}

	//LODs, if we deleted
	for (i = 0; i<header.detail_levels.size(); i++){
		if (header.detail_levels[i] == idx)
			header.detail_levels.erase(header.detail_levels.begin() + i);//does not exsist anymore
		if (i >= header.detail_levels.size())
			break;//we just deleted the last level
		if (header.detail_levels[i] >(int)idx)
			header.detail_levels[i]--;//it's actualy been moved one back, or it will be at the end of this function anyway
	}
	//debris, if we deleted
	for (i = 0; i<header.debris_pieces.size(); i++){
		if (header.debris_pieces[i] == idx)
			header.debris_pieces.erase(header.debris_pieces.begin() + i);//does not exsist anymore
		if (i >= header.debris_pieces.size())
			break;//we just deleted the last level
		if (header.debris_pieces[i] >(int)idx)
			header.debris_pieces[i]--;//it's actualy been moved one back, or it will be at the end of this function anyway
	}
	for (i = 0; i<turrets.size(); i++){
		if (turrets[i].sobj_par_phys == idx)
			turrets[i].sobj_par_phys = turrets[i].sobj_parent;
		if (turrets[i].sobj_par_phys >(int)idx)
			turrets[i].sobj_par_phys--;

		if (turrets[i].sobj_parent == idx)
			turrets.erase(turrets.begin() + i);//does not exsist anymore
		if (i >= turrets.size())
			break;//we just deleted the end
		if (turrets[i].sobj_parent > (int)idx)
			turrets[i].sobj_parent--;//it's actualy been moved one back, or it will be at the end of this function anyway
	}
	for (i = 0; i<eyes.size(); i++){
		if (eyes[i].sobj_number == idx)
			eyes[i].sobj_number = subobjects[idx].parent_sobj;
		if (i >= eyes.size())
			break;//we just deleted the end
		if (eyes[i].sobj_number >(int)idx)
			eyes[i].sobj_number--;
	}
	for (i = 0; i<light_arrays.size(); i++){
		if (light_arrays[i].obj_parent == idx)
			light_arrays[i].obj_parent = subobjects[idx].parent_sobj;
		if (i >= light_arrays.size())
			break;//we just deleted the end
		if (light_arrays[i].obj_parent >(int)idx)
			light_arrays[i].obj_parent--;
	}

	//finaly...
	if (active_submodel == idx)
		active_submodel = subobjects[idx].parent_sobj;
	if (active_submodel > (int)idx)
		active_submodel--;
	if (active_submodel < 0 && header.detail_levels.size() > 0)
		active_submodel = header.detail_levels[0];

	if (bsp_cache.size() > (unsigned)idx) {
		bsp_cache.erase(bsp_cache.begin() + idx);
	}
	RemoveIndex(subobjects, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::SetObjectChanged(unsigned int idx)
{
	if (idx >= subobjects.size())
		return;

	if (vbos_enabled) {
	//	make_vertex_buffer(idx);
	}
	if (can_bsp_cache)
		bsp_cache[idx].changed = true;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddEye(pcs_eye_pos *eye)
{
	PCS_ADD_TO_VEC(eyes, eye)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelEye(unsigned int idx)
{
	RemoveIndex(eyes, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddSpecial(pcs_special *spcl)
{
	PCS_ADD_TO_VEC(special, spcl)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelSpecial(unsigned int idx)
{
	RemoveIndex(special, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddWeapon(pcs_slot *weap)
{
	PCS_ADD_TO_VEC(weapons, weap)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelWeapon(unsigned int idx)
{
	RemoveIndex(weapons, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddTurret(pcs_turret *trrt)
{
	turrets.push_back(*trrt);
	//	PCS_ADD_TO_VEC(turrets, trrt)
	//was getting an invalid allocation size with this for some reason
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelTurret(unsigned int idx)
{
	RemoveIndex(turrets, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::AddDocking(pcs_dock_point *dock)
{
	PCS_ADD_TO_VEC(docking, dock)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelDocking(unsigned int idx)
{
	RemoveIndex(docking, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddThruster(pcs_thruster *thrust)
{
	PCS_ADD_TO_VEC(thrusters, thrust)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelThruster(unsigned int idx)
{
	RemoveIndex(thrusters, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddShldTri(pcs_shield_triangle *stri)
{
	PCS_ADD_TO_VEC(shield_mesh, stri)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelShldTri(unsigned int idx)
{
	RemoveIndex(shield_mesh, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddInsignia(pcs_insig *insig)
{
	PCS_ADD_TO_VEC(insignia, insig)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelInsignia(unsigned int idx)
{
	RemoveIndex(insignia, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddPath(pcs_path *path)
{
	PCS_ADD_TO_VEC(ai_paths, path)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelPath(unsigned int idx)
{
	RemoveIndex(ai_paths, idx);
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


void PCS_Model::AddLight(pcs_glow_array *lights)
{
	PCS_ADD_TO_VEC(light_arrays, lights)
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::DelLight(unsigned int idx)
{
	RemoveIndex(light_arrays, idx);
}


//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void PCS_Model::Calculate_Smoothing_Data(int &sobjs_comp)
{
	std::vector<std::vector<int> > covertals; // other polygons that share a vertex with us
	unsigned int i, j, k, l, m, cvc;
	bool tBool;//, afKill;
	vector3d tvect;//, *vects;
	pcs_sobj *sobj;
	//pcs_polygon *polA, *polB;
	if (has_fullsmoothing_data) return;

	sobjs_comp = 0;
	for (i = 0; i < this->subobjects.size(); i++) // subobjects
	{
		sobj = &this->subobjects[i];
		for (j = 0; j < sobj->polygons.size(); j++) // outer polygon
		{

			// ---- first.. let's see if we're faceted.. if that's the case we don't have to compile the covertals lists -------
			tBool = true;
			for (k = 0; k < sobj->polygons[j].verts.size() && tBool == true; k++)
			{
				tBool = tBool && (sobj->polygons[j].norm == sobj->polygons[j].verts[k].norm);
			}
			if (tBool)
			{
				// we're faceted
				for (k = 0; k < sobj->polygons[j].verts.size(); k++)
					sobj->polygons[j].verts[k].facet_angle = 0;
				continue;
			}

			//  ---------- find all polygons that share our verts -------
			covertals.resize(sobj->polygons[j].verts.size());
			for (k = 0; k < sobj->polygons[j].verts.size(); k++) // outer vertex
			{
				cvc = 0;
				covertals[k].resize(10); // reasonable start

				for (l = 0; l < sobj->polygons.size(); l++) // inner polygon
				{
					if (l == j)
						continue;
					for (m = 0; m < sobj->polygons[l].verts.size(); m++) // inner vertex
					{
						if (sobj->polygons[j].verts[k].point == sobj->polygons[l].verts[m].point)
						{
							covertals[k][cvc] = l;
							cvc++;

							if (cvc >= covertals[k].size())
								covertals[k].resize(cvc * 2);
						}
					} // inner vertex
				} // inner polygon

				covertals[k].resize(cvc);
			} // outer vertex

			// ok.. we have all our vertex lists covertals... now let's try to figure out what we are

			// -- trying full smooth --
			tBool = true;
			for (k = 0; k < sobj->polygons[j].verts.size() && tBool == true; k++)
			{
				tvect = sobj->polygons[j].norm;

				for (l = 0; l < covertals[k].size(); l++)
				{
					tvect += sobj->polygons[covertals[k][l]].norm;
				}
				tvect = tvect / (1 + covertals[k].size());
				tBool = tBool && (tvect == sobj->polygons[j].verts[k].norm);
			}
			if (tBool)
			{
				// we're full smooth
				for (k = 0; k < sobj->polygons[j].verts.size(); k++)
					sobj->polygons[j].verts[k].facet_angle = -1;
				continue;
			}

			// give up.. set autofacet to 32 degrees

			for (k = 0; k < sobj->polygons[j].verts.size(); k++)
				sobj->polygons[j].verts[k].facet_angle = 32;

			// ----- oh... crap... we're autofacet.... this is going to be UGLY --------
			// infact.. we're NOT going to try this.. it would be ungodly bad and probably never find the right angle to match
			//vector3d average_vectors_if_less_than_angle(int numvectors, float angle, vector3d src, vector3d *vects);
			/*afKill = false;
			for (k = 0; k < 360 && !afKill; k) // progressively increase the angle until we find it
			{
			tBool = true;
			for (l = 0; l < sobj->polygons[j].verts.size(); l++)
			{
			vects = new vector3d[1+covertals[l].size()];
			vects[0] = sobj->polygons[j].norm;

			for (m = 0; m < covertals[l].size(); m++)
			{
			vects[m+1] = sobj->polygons[covertals[k][l]].norm
			}

			delete[] vects;
			}
			} //angle loop
			*/


		} // outer polygon

		sobjs_comp++;
	} // subobjects

	has_fullsmoothing_data = true;
}


//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

vector3d PCS_Model::OffsetFromParent(int ObjNum)
{
	vector3d RetVal = MakeVector(0, 0, 0);
	int parnum;

	parnum = subobjects[ObjNum].parent_sobj;
	if (parnum == -1)
		return RetVal;
	else
	{
		RetVal = subobjects[ObjNum].offset;
		return RetVal + OffsetFromParent(parnum);
	}
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

int PCS_Model::FindTexture(std::string name)
{
	for (unsigned int i = 0; i < textures.size(); i++)
		if (textures[i] == name)
			return i;
	return -1;
}

bobboau::matrix PCS_Model::moi_recalculate(int Xres, int Zres){
	float xw = (header.max_bounding.x - header.min_bounding.x) / (Xres);
	float zw = (header.max_bounding.z - header.min_bounding.z) / (Zres);

	bobboau::matrix ret = bobboau::matrix(vector3d(), vector3d(), vector3d());

	for (int x = 0; x<Xres; x++){
		for (int z = 0; z<Zres; z++){
			ret = ret + (moi_recalculate(header.min_bounding.x + x*xw + xw / 2.0f, header.min_bounding.z + z*zw + zw / 2.0f, xw, zw));
		}
	}

	/*
	ret = moi_recalculate((header.max_bounding.x + header.min_bounding.x)/2
	,(header.max_bounding.z + header.min_bounding.z)/2
	,(header.max_bounding.x - header.min_bounding.x)
	,(header.max_bounding.z - header.min_bounding.z)
	,header.min_bounding.y,header.max_bounding.y-header.min_bounding.y);
	*/
	ret = ret * header.mass;
	//	return ret;
	return ret.invert();
}

bobboau::matrix PCS_Model::moi_recalculate(float X, float Z, float xw, float zw){
	std::vector<vector3d> cpoints;
	if (!moi_colide(cpoints, X, Z))
		return bobboau::matrix(vector3d(), vector3d(), vector3d());//return null
	//moi colide ensures there is an even number of colisions and orders them into pairs
	bobboau::matrix ret = moi_recalculate(X, Z, xw, zw, cpoints[0].y, cpoints[1].y - cpoints[0].y);
	for (unsigned int i = 1; i<cpoints.size() / 2; i++){
		ret = ret + moi_recalculate(X, Z, xw, zw, cpoints[i].y, cpoints[i + 1].y - cpoints[i].y);
	}
	return ret;
}

bobboau::matrix PCS_Model::moi_recalculate(float X, float Z, float xw, float zw, float Y, float dy){
	bobboau::matrix ret;
	float m = 1.0;//dy*xw*zw/2.0;
	float ftm = 4.0 / 3.0*m;
	float X2 = X*X;
	float Y2 = Y*Y;
	float Z2 = Z*Z;
	float xw2 = xw*xw;
	float zw2 = zw*zw;
	float dy2 = dy*dy;
	ret.a2d[0][0] = ftm*(3.0*dy*Y + dy2 + zw2 + 3.0*Y2 + 12.0*Z2);
	ret.a2d[1][1] = ftm*(xw2 + zw2 + 12.0*X2 + 12.0*Z2);
	ret.a2d[2][2] = ftm*(3.0*dy*Y + dy2 + xw2 + 12.0*X2 + 3.0*Y2);

	ret.a2d[0][1] = ret.a2d[1][0] = -4.0*m*X*(dy + 2.0*Y);
	ret.a2d[0][2] = ret.a2d[2][0] = -16.0*m*X*Z;
	ret.a2d[2][1] = ret.a2d[1][2] = -4.0*m*Z*(dy + 2.0*Y);

	return ret;
}

vector3d poly_min(pcs_polygon&poly){
	vector3d ret;
	ret = poly.verts[0].point;
	for (unsigned int i = 1; i<poly.verts.size(); i++){
		vector3d&point = poly.verts[i].point;
		if (point.x < ret.x)ret.x = point.x;
		if (point.y < ret.y)ret.y = point.y;
		if (point.z < ret.z)ret.z = point.z;
	}
	return ret;
}

vector3d poly_max(pcs_polygon&poly){
	vector3d ret;
	ret = poly.verts[0].point;
	for (unsigned int i = 1; i<poly.verts.size(); i++){
		vector3d&point = poly.verts[i].point;
		if (point.x > ret.x)ret.x = point.x;
		if (point.y > ret.y)ret.y = point.y;
		if (point.z > ret.z)ret.z = point.z;
	}
	return ret;
}

int moicpcf(const void*a, const void*b){
	return ((vector3d*)(a))->y < ((vector3d*)(b))->y;
}

int point_face(vector3d *checkp, std::vector<pcs_vertex> verts, vector3d * norm1){
	std::vector<vector3d> v(verts.size());
	for (unsigned int i = 0; i< verts.size(); i++){
		v[i] = verts[i].point;
	}
	return point_face(checkp, v, norm1);
}

bool PCS_Model::moi_colide(std::vector<vector3d>&cpoints, float x, float z){
	pcs_sobj&model = subobjects[header.detail_levels[0]];
	for (unsigned int i = 0; i<model.polygons.size(); i++){
		pcs_polygon&poly = model.polygons[i];
		vector3d min = poly_min(poly);
		vector3d max = poly_max(poly);
		if (x>max.x || x<min.x)continue;
		if (z>max.z || z<min.z)continue;

		bool sucsess;
		vector3d cpoint = plane_line_intersect(poly.verts[0].point, poly.norm, vector3d(x, 0.0f, z), vector3d(0.0f, 1.0f, 0.0f), &sucsess);
		if (!sucsess)continue;

		if (!point_face(&cpoint, poly.verts, &poly.norm))
			continue;

		cpoints.push_back(cpoint);
	}

	if (cpoints.size() % 2)
		cpoints.resize(0);

	if (cpoints.size() > 0){
		qsort(&cpoints[0], cpoints.size(), sizeof(vector3d), moicpcf);
		return true;
	}
	else{
		return false;
	}
}

void PCS_Model::Transform(const matrix& transform, const vector3d& translation) {
	std::set<int> dock_paths;
	header.min_bounding = vector3d(FLT_MAX, FLT_MAX, FLT_MAX);
	header.max_bounding = vector3d(FLT_MIN, FLT_MIN, FLT_MIN);
	header.max_radius = 0.0f;
	for (std::vector<pcs_sobj>::iterator it = subobjects.begin(); it < subobjects.end(); ++it) {
		if (it->parent_sobj == -1) {
			it->Transform(*this, (int)(it - subobjects.begin()), transform, translation, true, true);
		}
	}
	if (header.max_radius == 0.0f) {
		header.min_bounding = vector3d(0.0f, 0.0f, 0.0f);
		header.max_bounding = vector3d(0.0f, 0.0f, 0.0f);
	}
	for (std::vector<pcs_special>::iterator it = special.begin(); it < special.end(); ++it) {
		it->Transform(*this, transform, translation);
	}
	for (std::vector<pcs_slot>::iterator it = weapons.begin(); it < weapons.end(); ++it) {
		it->Transform(transform, translation);
	}
	for (std::vector<pcs_dock_point>::iterator it = docking.begin(); it < docking.end(); ++it) {
		it->Transform(*this, transform, translation);
		for (std::vector<int>::iterator jt = it->paths.begin(); jt < it->paths.end(); ++jt) {
			dock_paths.insert(*jt);
		}
	}
	for (std::vector<pcs_thruster>::iterator it = thrusters.begin(); it < thrusters.end(); ++it) {
		it->Transform(transform, translation);
	}
	for (std::vector<pcs_shield_triangle>::iterator it = shield_mesh.begin(); it < shield_mesh.end(); ++it) {
		it->Transform(*this, transform, translation);
	}
	for (std::vector<pcs_insig>::iterator it = insignia.begin(); it < insignia.end(); ++it) {
		it->Transform(transform, translation);
	}
	for (std::vector<pcs_path>::iterator it = ai_paths.begin(); it < ai_paths.end(); ++it) {
		if (it->parent.empty() && dock_paths.find((int)(it - ai_paths.begin())) != dock_paths.end()) {
			it->Transform(transform, translation);
		}
	}
	header.mass_center = transform * header.mass_center + translation;
	header.mass *= fabs(transform.determinant());
	header.max_radius_override *= fabs(transform.determinant());
	header.min_bounding_override = transform * header.min_bounding_override + translation;
	header.max_bounding_override = transform * header.max_bounding_override + translation;
}


//****************************************************************************
// Parallax Object File (POF) Saver
//****************************************************************************

int PCS_Model::SaveToPOF(std::string filename, AsyncProgress* progress)
{
	PCS_Model::BSP_MAX_DEPTH = 0;
	PCS_Model::BSP_NODE_POLYS = 1;
	PCS_Model::BSP_TREE_TIME = 0;
	PCS_Model::BSP_COMPILE_ERROR = false;
	POF poffile;
	unsigned int i, j, k, l;
	progress->setTarget(6 + light_arrays.size() + ai_paths.size() + insignia.size() + shield_mesh.size() +
		thrusters.size() + docking.size() + turrets.size() + weapons.size() + special.size() +
		eyes.size() + model_info.size() + subobjects.size() + textures.size());
	char cstringtemp[256];


	// Update Progress
	progress->incrementWithMessage("Writing Header Pt1");

	// --------- convert cross sections --------- 
	std::vector<cross_section> sections;
	sections.resize(header.cross_sections.size());

	for (i = 0; i < header.cross_sections.size(); i++)
	{
		sections[i].depth = header.cross_sections[i].depth;
		sections[i].radius = header.cross_sections[i].radius;
	}
	poffile.HDR2_Set_CrossSections(header.cross_sections.size(), sections);

	// Update Progress
	progress->incrementWithMessage("Writing Header Pt2");

	// --------- ACEN --------- 
	poffile.ACEN_Set_acen(POFTranslate(autocentering));


	// Update Progress
	progress->incrementWithMessage("Writing Acen");
	// --------- TXTR --------- 


	// Update Progress
	progress->incrementWithMessage("Writing Textures");
	for (i = 0; i < textures.size(); i++)
		poffile.TXTR_AddTexture(textures[i]);


	// --------- Sub object Consversion ---------

	//wxLongLong time = wxGetLocalTimeMillis();
	bool bsp_compiled = false;
	header.max_radius = 0.0f;
	for (i = 0; i < subobjects.size(); i++)
	{
		// Update Progress
		//if (subobjects[i].name == "debris08")
		//	sprintf(cstringtemp, "Submodel %d: %s SENTINAL!", i, subobjects[i].name.c_str());
		//else
		sprintf(cstringtemp, "Submodel %d: %s", i, subobjects[i].name.c_str());
		progress->incrementWithMessage(cstringtemp);

		//memset((char *)&obj, 0, sizeof(OBJ2)); this is NO LONGER ALLOWED - obj2 now contains an class w/ vtable
		boost::scoped_ptr<OBJ2> obj(new OBJ2);
		obj->submodel_number = i;
		if (!PMFObj_to_POFObj2(i, *obj, bsp_compiled, header.max_radius))
		{
			return 2; // error occured in bsp splitting!
		}
		poffile.OBJ2_Add(*obj); // takes over object management - including pointers
	}
	//time = wxGetLocalTimeMillis() - time;

	// we succeeded in compiling - let's cache the result
	can_bsp_cache = true;
	bsp_cache.resize(subobjects.size());
	for (i = 0; i < subobjects.size(); i++)
		poffile.OBJ2_Get_BSPData(i, bsp_cache[i].bsp_data);


	// --------- ---------------------- ---------


	int idx = GetModelInfoCount();
	char cstrtmp[256];
	//wxString strtmp = PCS2_VERSION;
	//sprintf(cstrtmp, "PMFSaveToPOF: Compiled on %s with %s\nmax BSP depth was %d\nmost polys in a single node was %d\nTotal Compile time was %ldms, tree generation time was %ldms", std::string(strtmp.mb_str()).c_str(), std::string(PCS2_COMP_VERSION.mb_str()).c_str(), PCS_Model::BSP_MAX_DEPTH, PCS_Model::BSP_NODE_POLYS, time.ToLong(), PCS_Model::BSP_TREE_TIME.ToLong());

	bool found = false;
	for (i = 0; i < model_info.size() && !found; i++)
	{
		if (strstr(model_info[i].c_str(), "PMFSaveToPOF") != NULL)
		{
			found = true;
			if (bsp_compiled) // update the string
				model_info[i] = cstrtmp;
		}
	}

	if (!found)
		AddModelInfo(cstrtmp);

	j = 0;
	for (i = 0; i < model_info.size(); i++)
		j += model_info[i].length() + 1;

	boost::scoped_ptr<char> pinf(new char[j]);
	memset(pinf.get(), 0, j);
	j = 0;

	for (i = 0; i < model_info.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing String %d", i);
		progress->incrementWithMessage(cstringtemp);

		strncpy(pinf.get() + j, model_info[i].c_str(), model_info[i].length());
		j += model_info[i].length() + 1;
	}
	poffile.PINF_Set(pinf.get(), j);

	if (found)
		model_info.resize(idx); // back down to size

	// ---------  EYE --------- 

	for (i = 0; i < eyes.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Eye %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.EYE_Add_Eye(eyes[i].sobj_number,
			POFTranslate(eyes[i].sobj_offset),
			POFTranslate(eyes[i].normal));
	}


	// --------- SPCL --------- 	

	for (i = 0; i < special.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Special %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.SPCL_AddSpecial(special[i].name, special[i].properties,
			POFTranslate(special[i].point), special[i].radius);
	}

	k = l = 0;
	// --------- weapons (GPNT/MPNT) --------- 
	for (i = 0; i < weapons.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Weapon %d", i);
		progress->incrementWithMessage(cstringtemp);
		if (weapons[i].type == GUN)
		{
			poffile.GPNT_AddSlot();

			for (j = 0; j < weapons[i].muzzles.size(); j++)
			{
				poffile.GPNT_AddPoint(k, POFTranslate(weapons[i].muzzles[j].point),
					POFTranslate(weapons[i].muzzles[j].norm));
			}
			k++;
		}
		else
		{
			poffile.MPNT_AddSlot();

			for (j = 0; j < weapons[i].muzzles.size(); j++)
			{
				poffile.MPNT_AddPoint(l, POFTranslate(weapons[i].muzzles[j].point),
					POFTranslate(weapons[i].muzzles[j].norm));
			}
			l++;
		}
	}

	// --------- turrets TGUN/TMIS --------- 
	k = l = 0;

	for (i = 0; i < turrets.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Turret %d", i);
		progress->incrementWithMessage(cstringtemp);
		if (turrets[i].type == GUN)
		{
			poffile.TGUN_Add_Bank(turrets[i].sobj_parent,
				turrets[i].sobj_par_phys,
				POFTranslate(turrets[i].turret_normal));
			for (j = 0; j < turrets[i].fire_points.size(); j++)
			{
				poffile.TGUN_Add_FirePoint(k, POFTranslate(turrets[i].fire_points[j]));
			}
			k++;
		}
		else
		{
			poffile.TMIS_Add_Bank(turrets[i].sobj_parent,
				turrets[i].sobj_par_phys,
				POFTranslate(turrets[i].turret_normal));
			for (j = 0; j < turrets[i].fire_points.size(); j++)
			{
				poffile.TMIS_Add_FirePoint(l, POFTranslate(turrets[i].fire_points[j]));
			}
			l++;
		}
	}

	// --------- docking --------- 
	for (i = 0; i < docking.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Docking %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.DOCK_Add_Dock(docking[i].properties);

		for (j = 0; j < docking[i].dockpoints.size(); j++)
		{
			poffile.DOCK_Add_Point(i, POFTranslate(docking[i].dockpoints[j].point),
				POFTranslate(docking[i].dockpoints[j].norm));
		}

		for (j = 0; j < docking[i].paths.size(); j++)
		{
			poffile.DOCK_Add_SplinePath(i, docking[i].paths[j]);
		}
	}

	// --------- thrusters --------- 
	for (i = 0; i < thrusters.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Thruster %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.FUEL_Add_Thruster(thrusters[i].properties);

		for (j = 0; j < thrusters[i].points.size(); j++)
		{
			poffile.FUEL_Add_GlowPoint(i, thrusters[i].points[j].radius,
				POFTranslate(thrusters[i].points[j].pos),
				POFTranslate(thrusters[i].points[j].norm));
		}
	}

	// --------- shield_mesh --------- 
	int fcs[3], nbs[3];
	std::vector<vector3d> points(shield_mesh.size() * 3);
	vector3d shldtvect;

	// make points list
	l = 0;
	for (i = 0; i < shield_mesh.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			bool found = false;
			for (unsigned int k = 0; k < points.size(); k++)
			{
				if (points[k] == POFTranslate(shield_mesh[i].corners[j]))
				{
					found = true;
					break;
				}
			}
			if (!found)
			{
				if (l >= points.size())
					points.resize(points.size() * 2);
				points[l] = POFTranslate(shield_mesh[i].corners[j]);
				l++;
			}
		}
	}
	points.resize(l);

	// translate shield mesh
	for (i = 0; i < shield_mesh.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Shield Tri %d", i);
		progress->incrementWithMessage(cstringtemp);
		// adds points to list if need, determine fcs[]
		for (j = 0; j < 3; j++)
		{
			shldtvect = POFTranslate(shield_mesh[i].corners[j]);
			fcs[j] = FindInList(points, shldtvect);
		}

		// determine neighbors (nbs[])
		j = 0;
		for (k = 0; k < shield_mesh.size() && j < 3; k++)
		{
			if (Neighbor(shield_mesh[i], shield_mesh[k]) && i != k)
			{
				nbs[j] = k;
				j++;
			}
		}
		// add
		poffile.SHLD_Add_Face(POFTranslate(shield_mesh[i].face_normal), fcs, nbs);
	}

	// add points
	// Update Progress
	progress->incrementWithMessage("Writing Shield Points");
	for (i = 0; i < points.size(); i++)
		poffile.SHLD_Add_Vertex(points[i]);


	progress->incrementWithMessage("Writing Shield Collision Tree");
	// --------------- generate shield collision tree ----------------
	if (poffile.SHLD_Count_Faces() > 0)
	{
		std::vector<pcs_polygon> shldmesh(poffile.SHLD_Count_Faces());

		// make a pcs_polygon mesh from the shields
		for (i = 0; i < shldmesh.size(); i++)
		{
			shldmesh[i].verts.resize(3);

			poffile.SHLD_Get_Face(i, shldmesh[i].norm, fcs, nbs);

			for (j = 0; j < 3; j++)
			{
				poffile.SHLD_Get_Vertex(fcs[j], shldmesh[i].verts[j].point);
				shldmesh[i].verts[j].norm = shldmesh[i].norm;
			}

			shldmesh[i].centeroid = PolygonCenter(shldmesh[i]);
		}

		// make the tree
		vector3d smin, smax;
		std::unique_ptr<bsp_tree_node> shld_root = MakeTree(shldmesh, smax, smin);

		// pack the tree
		int sldc_size = CalcSLDCTreeSize(shld_root.get());
		std::vector<char> sldc;
		sldc.resize(sldc_size);

		PackTreeInSLDC(shld_root.get(), 0, &sldc.front(), sldc_size);

		poffile.SLDC_SetTree(std::move(sldc)); // POFHandler will steal our copy of the buffer
	}

	// --------- insignia --------- 

	vector3d uv, vv;
	float *u = (float *)&uv, *v = (float *)&vv;
	for (i = 0; i < insignia.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Insignia %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.INSG_Add_insignia(insignia[i].lod, POFTranslate(insignia[i].offset));

		for (j = 0; j < insignia[i].faces.size(); j++)
		{
			for (k = 0; k < 3; k++)
			{
				while ((l = poffile.INST_Find_Vert(i, POFTranslate(insignia[i].faces[j].verts[k]))) == (unsigned)-1)
				{
					poffile.INSG_Add_Insig_Vertex(i, POFTranslate(insignia[i].faces[j].verts[k]));
				}
				fcs[k] = l;
				u[k] = insignia[i].faces[j].u[k];
				v[k] = insignia[i].faces[j].v[k];
			}
			poffile.INSG_Add_Insig_Face(i, fcs, uv, vv);
		}
	}

	// --------- ai_paths --------- 
	for (i = 0; i < ai_paths.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Path %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.PATH_Add_Path(ai_paths[i].name, ai_paths[i].parent);

		for (j = 0; j < ai_paths[i].verts.size(); j++)
		{
			poffile.PATH_Add_Vert(i, POFTranslate(ai_paths[i].verts[j].pos),
				ai_paths[i].verts[j].radius);
		}
	}

	// --------- light arrays --------- 
	pcs_glow_array *gla;
	for (i = 0; i < light_arrays.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Writing Glow %d", i);
		progress->incrementWithMessage(cstringtemp);
		gla = &light_arrays[i];
		poffile.GLOW_Add_LightGroup(gla->disp_time, gla->on_time, gla->off_time,
			gla->obj_parent, gla->LOD, gla->type, gla->properties);
		for (j = 0; j < gla->lights.size(); j++)
		{
			poffile.GLOW_Add_GlowPoint(i, gla->lights[j].radius,
				POFTranslate(gla->lights[j].pos),
				POFTranslate(gla->lights[j].norm));
		}
	}

	// --------- HDR2 --------- 

	// let's make new BBoxes based on the subobjects
	vector3d minbox, maxbox, tmpmin, tmpmax;
	poffile.OBJ2_Get_BoundingMax(0, maxbox);
	poffile.OBJ2_Get_BoundingMin(0, minbox);

	for (i = 1; i < poffile.OBJ2_Count(); i++)
	{
		vector3d sobj_offset(POFTranslate(OffsetFromParent(i)));
		poffile.OBJ2_Get_BoundingMax(i, tmpmax);
		poffile.OBJ2_Get_BoundingMin(i, tmpmin);
		ExpandBoundingBoxes(maxbox, minbox, tmpmax + sobj_offset);
		ExpandBoundingBoxes(maxbox, minbox, tmpmin + sobj_offset);

	}

	for (i = 0; i < poffile.SHLD_Count_Vertices(); i++)
	{
		poffile.SHLD_Get_Vertex(i, tmpmax);
		ExpandBoundingBoxes(maxbox, minbox, tmpmax);
	}

	// update our bounding boxes
	//axis doesn't matter on the bounding boxes - it's all negative/positive values! i'm a faqing moron
	poffile.HDR2_Set_MinBound(header.min_bounding_overridden ? header.min_bounding_override : minbox);
	poffile.HDR2_Set_MaxBound(header.max_bounding_overridden ? header.max_bounding_override : maxbox);
	this->header.max_bounding = minbox;
	this->header.min_bounding = maxbox;

	poffile.HDR2_Set_MaxRadius(header.max_radius_overridden ? header.max_radius_override : header.max_radius);
	poffile.HDR2_Set_Details(header.detail_levels.size(), header.detail_levels);
	poffile.HDR2_Set_Debris(header.debris_pieces.size(), header.debris_pieces);
	poffile.HDR2_Set_Mass(header.mass);
	poffile.HDR2_Set_MassCenter(POFTranslate(header.mass_center));
	poffile.HDR2_Set_MomentInertia(header.MOI);
	poffile.HDR2_Set_SOBJCount(GetSOBJCount());

	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);

	if (!poffile.SavePOF(out))
		return 1;

	return 0;
}


//****************************************************************************
// Parallax Object File (POF) Loader
//****************************************************************************

int PCS_Model::LoadFromPOF(std::string filename, AsyncProgress* progress)
{
	this->Reset();
	char cstringtemp[256];
	progress->setMessage("Opening and Reading POF");
	progress->Notify();

	std::ifstream infile(filename.c_str(), std::ios::in | std::ios::binary);
	if (!infile)
		return 1;

	POF poffile(infile);
	progress->setTarget(4 + poffile.SPCL_Count() + poffile.GPNT_SlotCount() + poffile.MPNT_SlotCount() + poffile.TGUN_Count_Banks() + poffile.TMIS_Count_Banks() +
		poffile.DOCK_Count_Docks() + poffile.FUEL_Count_Thrusters() + poffile.INSG_Count_Insignia() + poffile.PATH_Count_Paths() +
		poffile.GLOW_Count_LightGroups() + poffile.OBJ2_Count());

	// Update Progress
	progress->incrementWithMessage("Getting Header");

	header.max_radius = poffile.HDR2_Get_MaxRadius();
	header.max_radius_override = header.max_radius;
	header.min_bounding = poffile.HDR2_Get_MinBound();
	header.max_bounding = poffile.HDR2_Get_MaxBound();
	POFTranslateBoundingBoxes(header.min_bounding, header.max_bounding);
	header.min_bounding_override = header.min_bounding;
	header.max_bounding_override = header.max_bounding;

	unsigned int i, j, k;
	int scratch; // useless variable - for legacy remnant argument
	poffile.HDR2_Get_Details(scratch, header.detail_levels);
	poffile.HDR2_Get_Debris(scratch, header.debris_pieces);

	header.mass = poffile.HDR2_Get_Mass();
	header.mass_center = POFTranslate(poffile.HDR2_Get_MassCenter());
	poffile.HDR2_Get_MomentInertia(header.MOI);

	// --------- convert cross sections --------- 
	std::vector<cross_section> sections;
	poffile.HDR2_Get_CrossSections(scratch, sections);

	if (scratch != -1)
	{
		header.cross_sections.resize(scratch);

		for (i = 0; i < header.cross_sections.size(); i++)
		{
			header.cross_sections[i].depth = sections[i].depth;
			header.cross_sections[i].radius = sections[i].radius;
		}
	}
	// --------- ACEN --------- 
	// Update Progress
	progress->incrementWithMessage("Getting ACEN, TXTR, PINF, EYE");
	autocentering = POFTranslate(poffile.ACEN_Get_acen());

	// --------- TXTR --------- 
	textures.resize(poffile.TXTR_Count_Textures());
	std::string tmp_test;
	for (i = 0; i < textures.size(); i++)
	{
		tmp_test = poffile.TXTR_GetTextures(i);
		cout << tmp_test << endl;
		textures[i] = tmp_test;
	}

	// --------- ---------------------- ---------

	model_info = poffile.PINF_Get();

	can_bsp_cache = false;
	//for (i = 0; i < model_info.size(); i++)
	//{
	//	if ( //strstr(model_info[i].c_str(), "BSPGEN") || // Volition's Compiler - caching revoked 2008-02-11 by Kazan because cannot gaurantee that tagged models actually game from V's compiler
	//		//strstr(model_info[i].c_str(), "POF-CS Compiler v1.3.4") ||  // removed PCS1 from recognized cacheable list 2008-01-10 - Kazan
	//		strstr(model_info[i].c_str(), ""))
	//	{
	//		can_bsp_cache = true;
	//		break;
	//	}
	//}


	// ---------  EYE --------- 
	eyes.resize(poffile.EYE_Count_Eyes());

	for (i = 0; i < eyes.size(); i++)
	{
		poffile.EYE_Get_Eye(i, eyes[i].sobj_number, eyes[i].sobj_offset, eyes[i].normal);
		eyes[i].sobj_offset = POFTranslate(eyes[i].sobj_offset);
		eyes[i].normal = POFTranslate(eyes[i].normal);
	}

	// --------- SPCL --------- 
	special.resize(poffile.SPCL_Count());

	for (i = 0; i < special.size(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Special %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.SPCL_Get_Special(i, special[i].name, special[i].properties,
			special[i].point, special[i].radius);
		special[i].point = POFTranslate(special[i].point);
	}


	// --------- weapons (GPNT/MPNT) --------- 
	weapons.resize(poffile.GPNT_SlotCount() + poffile.MPNT_SlotCount());

	for (i = 0; i < poffile.GPNT_SlotCount(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Gun Point %d", i);
		progress->incrementWithMessage(cstringtemp);
		weapons[i].type = GUN;
		weapons[i].muzzles.resize(poffile.GPNT_PointCount(i));

		for (j = 0; j < poffile.GPNT_PointCount(i); j++)
		{
			poffile.GPNT_GetPoint(i, j, weapons[i].muzzles[j].point,
				weapons[i].muzzles[j].norm);
			weapons[i].muzzles[j].point = POFTranslate(weapons[i].muzzles[j].point);
			weapons[i].muzzles[j].norm = POFTranslate(weapons[i].muzzles[j].norm);
		}
	}

	k = poffile.GPNT_SlotCount();
	for (i = 0; i < poffile.MPNT_SlotCount(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Missile Point %d", i);
		progress->incrementWithMessage(cstringtemp);
		weapons[i + k].type = MISSILE;
		weapons[i + k].muzzles.resize(poffile.MPNT_PointCount(i));

		for (j = 0; j < poffile.MPNT_PointCount(i); j++)
		{
			poffile.MPNT_GetPoint(i, j, weapons[i + k].muzzles[j].point,
				weapons[i + k].muzzles[j].norm);
			weapons[i + k].muzzles[j].point = POFTranslate(weapons[i + k].muzzles[j].point);
			weapons[i + k].muzzles[j].norm = POFTranslate(weapons[i + k].muzzles[j].norm);
		}
	}

	// --------- turrets TGUN/TMIS --------- 
	//Matt  -- turrets are allocated here in the resize fn
	turrets.resize(poffile.TGUN_Count_Banks() + poffile.TMIS_Count_Banks());

	for (i = 0; i < poffile.TGUN_Count_Banks(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Gun Turret %d", i);
		progress->incrementWithMessage(cstringtemp);
		turrets[i].type = GUN;
		//turrets[i].name = string("GunTurret" + i);
		//turrets[i].turretIndex = i;
		//stuff parent idx and parent phys idx and normal into turret
		poffile.TGUN_Get_Bank(i, turrets[i].sobj_parent, turrets[i].sobj_par_phys, turrets[i].turret_normal);
		//translate the normal
		turrets[i].turret_normal = POFTranslate(turrets[i].turret_normal);
		//allocate fire points
		turrets[i].fire_points.resize(poffile.TGUN_Count_Points(i));
		
		//populate firepoints
		for (j = 0; j < poffile.TGUN_Count_Points(i); j++)
		{
			poffile.TGUN_Get_FirePoint(i, j, turrets[i].fire_points[j]);
			turrets[i].fire_points[j] = POFTranslate(turrets[i].fire_points[j]);
		}
	}

	k = poffile.TGUN_Count_Banks(); //array offset
	for (i = 0; i < poffile.TMIS_Count_Banks(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Missile Turret %d (%d)", i, i + k);
		progress->incrementWithMessage(cstringtemp);
		turrets[i + k].type = 1; //maybe this should be missile?
		poffile.TMIS_Get_Bank(i, turrets[i + k].sobj_parent, turrets[i + k].sobj_par_phys, turrets[i + k].turret_normal);

		turrets[i + k].turret_normal = POFTranslate(turrets[i + k].turret_normal);

		turrets[i + k].fire_points.resize(poffile.TMIS_Count_Points(i));

		for (j = 0; j < poffile.TMIS_Count_Points(i); j++)
		{
			poffile.TMIS_Get_FirePoint(i, j, turrets[i + k].fire_points[j]);
			turrets[i + k].fire_points[j] = POFTranslate(turrets[i + k].fire_points[j]);
		}
	}

	// --------- docking --------- 
	docking.resize(poffile.DOCK_Count_Docks());

	for (i = 0; i < poffile.DOCK_Count_Docks(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Docking Point %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.DOCK_Get_DockProps(i, docking[i].properties);

		docking[i].dockpoints.resize(poffile.DOCK_Count_Points(i));

		for (j = 0; j < poffile.DOCK_Count_Points(i); j++)
		{
			poffile.DOCK_Get_Point(i, j, docking[i].dockpoints[j].point,
				docking[i].dockpoints[j].norm);
			docking[i].dockpoints[j].point = POFTranslate(docking[i].dockpoints[j].point);
			docking[i].dockpoints[j].norm = POFTranslate(docking[i].dockpoints[j].norm);

		}

		docking[i].paths.resize(poffile.DOCK_Count_SplinePaths(i));

		for (j = 0; j < poffile.DOCK_Count_SplinePaths(i); j++)
		{
			poffile.DOCK_Get_SplinePath(i, j, docking[i].paths[j]);
		}
	}

	// --------- thrusters --------- 
	thrusters.resize(poffile.FUEL_Count_Thrusters());

	for (i = 0; i < poffile.FUEL_Count_Thrusters(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Thruster %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.FUEL_Get_ThrusterProps(i, thrusters[i].properties);

		thrusters[i].points.resize(poffile.FUEL_Count_Glows(i));
		for (j = 0; j < poffile.FUEL_Count_Glows(i); j++)
		{
			poffile.FUEL_Get_GlowPoint(i, j, thrusters[i].points[j].radius,
				thrusters[i].points[j].pos,
				thrusters[i].points[j].norm);
			thrusters[i].points[j].pos = POFTranslate(thrusters[i].points[j].pos);
			thrusters[i].points[j].norm = POFTranslate(thrusters[i].points[j].norm);
		}
	}

	// --------- shield_mesh --------- 
	// Update Progress
	progress->incrementWithMessage("Getting Shields");
	shield_mesh.resize(poffile.SHLD_Count_Faces());
	int fcs[3], nbs[3];
	for (i = 0; i < shield_mesh.size(); i++)
	{
		poffile.SHLD_Get_Face(i, shield_mesh[i].face_normal, fcs, nbs);

		shield_mesh[i].face_normal = POFTranslate(shield_mesh[i].face_normal); \

			for (j = 0; j < 3; j++)
			{
				poffile.SHLD_Get_Vertex(fcs[j], shield_mesh[i].corners[j]);
				shield_mesh[i].corners[j] = POFTranslate(shield_mesh[i].corners[j]);
			}
	}

	// --------- insignia --------- 
	insignia.resize(poffile.INSG_Count_Insignia());

	vector3d uv, vv;
	float *u = (float *)&uv, *v = (float *)&vv;

	for (i = 0; i < poffile.INSG_Count_Insignia(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Insignia %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.INSG_Get_Insignia(i, insignia[i].lod, insignia[i].offset);

		insignia[i].offset = POFTranslate(insignia[i].offset);

		insignia[i].faces.resize(poffile.INSG_Count_Faces(i));

		for (j = 0; j < poffile.INSG_Count_Faces(i); j++)
		{
			poffile.INSG_Get_Insig_Face(i, j, fcs, uv, vv);

			for (k = 0; k < 3; k++)
			{
				poffile.INSG_Get_Insig_Vertex(i, fcs[k], insignia[i].faces[j].verts[k]);
				insignia[i].faces[j].verts[k] = POFTranslate(insignia[i].faces[j].verts[k]);

				insignia[i].faces[j].u[k] = u[k];
				insignia[i].faces[j].v[k] = v[k];
			}
		}
	}

	// --------- ai_paths --------- 
	ai_paths.resize(poffile.PATH_Count_Paths());
	for (i = 0; i < poffile.PATH_Count_Paths(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Path %d", i);
		progress->incrementWithMessage(cstringtemp);
		poffile.PATH_Get_Path(i, ai_paths[i].name, ai_paths[i].parent);

		ai_paths[i].verts.resize(poffile.PATH_Count_Verts(i));

		for (j = 0; j < poffile.PATH_Count_Verts(i); j++)
		{
			poffile.PATH_Get_Vert(i, j, ai_paths[i].verts[j].pos,
				ai_paths[i].verts[j].radius);
			ai_paths[i].verts[j].pos = POFTranslate(ai_paths[i].verts[j].pos);
		}
	}

	// --------- light arrays --------- 
	light_arrays.resize(poffile.GLOW_Count_LightGroups());
	pcs_glow_array *gla;

	for (i = 0; i < poffile.GLOW_Count_LightGroups(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Glow Array %d", i);
		progress->incrementWithMessage(cstringtemp);
		gla = &light_arrays[i];
		poffile.GLOW_Get_Group(i, gla->disp_time, gla->on_time, gla->off_time,
			gla->obj_parent, gla->LOD, gla->type, gla->properties);
		gla->lights.resize(poffile.GLOW_Count_Glows(i));

		for (j = 0; j < poffile.GLOW_Count_Glows(i); j++)
		{
			poffile.GLOW_Get_GlowPoint(i, j, gla->lights[j].radius,
				gla->lights[j].pos,
				gla->lights[j].norm);
			gla->lights[j].pos = POFTranslate(gla->lights[j].pos);
			gla->lights[j].norm = POFTranslate(gla->lights[j].norm);
		}
	}

	// --------- Sub object Consversion ---------
	subobjects.resize(poffile.OBJ2_Count());

	if (can_bsp_cache)
		bsp_cache.resize(poffile.OBJ2_Count());
	for (i = 0; i < poffile.OBJ2_Count(); i++)
	{
		// Update Progress
		sprintf(cstringtemp, "Getting Object %d", i);
		//progress->incrementWithMessage(cstringtemp);

		pcs_sobj* obj = &subobjects[i];
		poffile.OBJ2_Get_Parent(i, obj->parent_sobj);
		poffile.OBJ2_Get_Radius(i, obj->radius);
		obj->radius_override = obj->radius;

		poffile.OBJ2_Get_Offset(i, obj->offset);
		obj->offset = POFTranslate(obj->offset);

		poffile.OBJ2_Get_GeoCenter(i, obj->geometric_center);
		obj->geometric_center = POFTranslate(obj->geometric_center);

		poffile.OBJ2_Get_BoundingMin(i, obj->bounding_box_min_point);

		poffile.OBJ2_Get_BoundingMax(i, obj->bounding_box_max_point);

		POFTranslateBoundingBoxes(obj->bounding_box_min_point, obj->bounding_box_max_point);
		obj->bounding_box_min_point_override = obj->bounding_box_min_point;
		obj->bounding_box_max_point_override = obj->bounding_box_max_point;

		poffile.OBJ2_Get_Name(i, obj->name);
		poffile.OBJ2_Get_Props(i, obj->properties);
		int type;
		poffile.OBJ2_Get_MoveType(i, type); // -1 = none, 1 = rotate
		switch (type)
		{
		case 1:
			obj->movement_type = ROTATE;
			break;
		default:
			obj->movement_type = MNONE;
		}
		poffile.OBJ2_Get_MoveAxis(i, type); // -1 = none, 1 = X, 2 = Z, 3 = Y
		switch (type)
		{
		case 0:
			obj->movement_axis = MV_X;
			break;
		case 1:
			obj->movement_axis = MV_Z;
			break;
		case 2:
			obj->movement_axis = MV_Y;
			break;
		default:
			obj->movement_axis = ANONE;
		}

		// fast method
		int bspsz;
		char *bspdata = NULL;
		poffile.OBJ2_Get_BSPDataPtr(i, bspsz, bspdata);

		//OBJ2_Get_BSPData

		obj->polygons.resize(100); // good starting size;

		unsigned int used_polygons = 0;
		BSP_DefPoints points;
		BSPTransPMF(0, (unsigned char *)bspdata, points, obj->polygons, used_polygons);
		obj->polygons.resize(used_polygons); // resize to exact size

		if (can_bsp_cache)
			poffile.OBJ2_Get_BSPData(i, bsp_cache[i].bsp_data);

		//if (bspdata)
		//	delete[] bspdata;
		//bspdata = NULL;
	}
	Transform(matrix(), vector3d());
	header.max_radius_overridden = fabs(header.max_radius - header.max_radius_override) > 0.0001f;
	header.max_bounding_overridden = header.max_bounding != header.max_bounding_override;
	header.min_bounding_overridden = header.min_bounding != header.min_bounding_override;
	for (size_t i = 0; i < subobjects.size(); ++i) {
		pcs_sobj& sobj = subobjects[i];
		sobj.radius_overridden = fabs(sobj.radius - sobj.radius_override) > 0.0001f;
		sobj.bounding_box_min_point_overridden = sobj.bounding_box_min_point != sobj.bounding_box_min_point_override;
		sobj.bounding_box_max_point_overridden = sobj.bounding_box_max_point != sobj.bounding_box_max_point_override;
	}
	return 0;
}

//****************************************************************************
// Parallax Object File (POF) Conversion Auxiliery Functions
//****************************************************************************

bool Neighbor(pcs_shield_triangle &face1, pcs_shield_triangle &face2)
{
	int CommonVerts = 0;

	if (face1.corners[0] == face2.corners[0] ||
		face1.corners[0] == face2.corners[1] ||
		face1.corners[0] == face2.corners[2])
		CommonVerts++;

	if (face1.corners[1] == face2.corners[0] ||
		face1.corners[1] == face2.corners[1] ||
		face1.corners[1] == face2.corners[2])
		CommonVerts++;

	if (face1.corners[2] == face2.corners[0] ||
		face1.corners[2] == face2.corners[1] ||
		face1.corners[2] == face2.corners[2])
		CommonVerts++;

	return (CommonVerts > 1 && CommonVerts < 3); // if CV=3 we're looking at ourself :D
}


//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

bool PCS_Model::PMFObj_to_POFObj2(int src_num, OBJ2 &dst, bool &bsp_compiled, float& model_radius)
{

	pcs_sobj &src = subobjects[src_num];

	dst.submodel_parent = src.parent_sobj;
	dst.offset = POFTranslate(src.offset);
	dst.geometric_center = POFTranslate(src.geometric_center);
	dst.submodel_name = APStoString(src.name);
	dst.properties = APStoString(src.properties);

	switch (src.movement_type)
	{
	case ROTATE:
		dst.movement_type = 1;
		break;
	default:
		dst.movement_type = -1;
	}
	switch (src.movement_axis)
	{
	case MV_X:
		dst.movement_axis = 0;
		break;
	case MV_Z:
		dst.movement_axis = 1;
		break;
	case MV_Y:
		dst.movement_axis = 2;
		break;
	default:
		dst.movement_axis = -1;
	}

	dst.reserved = 0;



	if (!can_bsp_cache || bsp_cache[src_num].changed)
	{

		// convert them to POF axis
		std::vector<pcs_polygon> clean_list = src.polygons;
		for (size_t i = 0; i < clean_list.size(); i++)
		{
			clean_list[i].norm = POFTranslate(clean_list[i].norm);
			for (size_t j = 0; j < clean_list[i].verts.size(); j++)
			{
				clean_list[i].verts[j].point = POFTranslate(clean_list[i].verts[j].point);
				clean_list[i].verts[j].norm = POFTranslate(clean_list[i].verts[j].norm);
			}
			clean_list[i].centeroid = PolygonCenter(clean_list[i]);
		}

		// BSP Compilation!
		// assemble points list
		std::vector<bsp_vert> points_list;
		std::vector<vector3d> pnts;
		std::unordered_map<vector3d, int> point_to_index;
		std::unordered_map<vector3d, int> normal_to_index;
		for (size_t i = 0; i < pnts.size(); i++) {
			point_to_index.insert(std::make_pair(pnts[i], i));
		}
		bsp_vert temp;
		points_list.reserve(clean_list.size());
		for (size_t i = 0; i < clean_list.size(); i++)
		{
			for (size_t j = 0; j < clean_list[i].verts.size(); j++)
			{
				auto point = point_to_index.find(clean_list[i].verts[j].point);
				if (point == point_to_index.end()) {
					point_to_index.insert(std::make_pair(clean_list[i].verts[j].point, points_list.size()));
					points_list.emplace_back();
					points_list.back().point = clean_list[i].verts[j].point;
					pnts.push_back(points_list.back().point);
				}
				auto normal = normal_to_index.find(clean_list[i].verts[j].norm);
				if (normal == normal_to_index.end()) {
					points_list[normal_to_index.size() / 128].norms.push_back(clean_list[i].verts[j].norm);
					normal_to_index.insert(std::make_pair(clean_list[i].verts[j].norm, normal_to_index.size()));
				}
			}
		}



		// create our defpoints
		BSP_DefPoints points;
		MakeDefPoints(points, points_list);
		vector3d AvgNormal;

		// create tree
		std::unique_ptr<bsp_tree_node> root = MakeTree(clean_list, dst.bounding_box_max_point, dst.bounding_box_min_point);

		// allocate buffer and write the defpoints
		dst.bsp_data.resize(points.head.size + CalculateTreeSize(root.get(), clean_list));

		if (points.Write(&dst.bsp_data.front()) != points.head.size)
			return false; // calculation error

		//std::ofstream bsp_debug("c:\\bsp.txt");
		//DebugPrintTree(root, bsp_debug);

		// pack the tree
		int error_flags = 0;
		PackTreeInBSP(root.get(), points.head.size, &dst.bsp_data.front(), clean_list, normal_to_index, point_to_index, points, dst.geometric_center, dst.bsp_data.size(), error_flags);

		// we got errors!
		if (error_flags != BSP_NOERRORS)
			return false;

		// update the bsp_compiled to be true
		bsp_compiled = true;

		// update the BSP cache with the new result
		if (can_bsp_cache)
		{
			// clear the saved - stale cache
			bsp_cache[src_num].decache();
			bsp_cache[src_num].bsp_data = dst.bsp_data;
			bsp_cache[src_num].changed = false;
		}


	}
	else // Used cached copy!
	{
		dst.bsp_data = bsp_cache[src_num].bsp_data;
	}
	dst.radius = 0.0f;
	dst.bounding_box_max_point = vector3d(FLT_MIN, FLT_MIN, FLT_MIN);
	dst.bounding_box_min_point = vector3d(FLT_MAX, FLT_MAX, FLT_MAX);

	vector3d global_offset(OffsetFromParent(src_num));
	for (unsigned int i = 0; i<src.polygons.size(); i++){
		for (unsigned int j = 0; j<src.polygons[i].verts.size(); j++){
			ExpandBoundingBoxes(dst.bounding_box_max_point, dst.bounding_box_min_point, src.polygons[i].verts[j].point);
			float norm = Magnitude(src.polygons[i].verts[j].point);
			if (norm > dst.radius) {
				dst.radius = norm;
			}
			float global_norm = Magnitude(src.polygons[i].verts[j].point + global_offset);
			if (global_norm > model_radius) {
				model_radius = global_norm;
			}

		}
	}
	if (dst.radius == 0.0f) {
		dst.bounding_box_max_point = vector3d();
		dst.bounding_box_min_point = vector3d();
	}
	if (src.radius_overridden) {
		dst.radius = src.radius_override;
	}
	if (src.bounding_box_min_point_overridden) {
		dst.bounding_box_min_point = src.bounding_box_min_point_override;
	}
	if (src.bounding_box_max_point_overridden) {
		dst.bounding_box_max_point = src.bounding_box_max_point_override;
	}
	POFTranslateBoundingBoxes(dst.bounding_box_min_point, dst.bounding_box_max_point);
	return true;
}


int PCS_Model::LoadFromDAE(std::string filename, AsyncProgress* progress, bool mirror_x, bool mirror_y, bool mirror_z) {
	DAEHandler dae_handler(filename, this, progress, mirror_x, mirror_y, mirror_z);
	return dae_handler.populate();
}

int PCS_Model::SaveToDAE(std::string filename, AsyncProgress* progress, int helpers, int props_as_helpers) {
	DAESaver dae_handler(filename, this, helpers, props_as_helpers, progress);
	return dae_handler.save();
}