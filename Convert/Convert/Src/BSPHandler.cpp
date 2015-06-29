
#if !defined(_WIN32)
#include "compat/filelength.h"
#include "compat/strncmp.h"
#endif
#include "BSPHandler.h"
#include <cstdio>
#include <memory.h>
#include <iostream>

using namespace std;
//********************************************************************************************************************

std::string BSP::DataIn(char *buffer, int size)
{
	std::string stats;
	char *localptr = buffer;
	char cstemp[64];
	bool go = true;
	BSP_BoundBox	bnd;
	BSP_DefPoints	pnt;
	BSP_FlatPoly	fpol;
	BSP_SortNorm	snrm;
	BSP_TmapPoly	tpol;

	BSP_BlockHeader head;

	stats = "";
	//0 - EOF, 1 - DEFPOINTS, 2 - FLATPOLY, 3 - IMAPPOLY, 4- SHORNORM, 5- BOUNDBOX
	while (go)
	{


		head.Read(localptr);
		localptr += 8;


		if (localptr - buffer >= size)
		{
			go = false;
			break;
		}
		switch (head.id)
		{
		case 0:
			//stats += "EOF\n";
			//stats += "\n";
			//go = false;
			break;

		case 1:

			memset(cstemp, 0, 64);
			sprintf(cstemp, "[%d bytes (%d bytes unused)]", head.size, head.size - pnt.Read(localptr, head));
			localptr += (head.size - 8);//8 for the already added header

			//stats += "DEFPOINTS ";
			//stats += cstemp;
			//stats += "\n";

			/*
			memset(cstemp, 0, 64);
			sprintf(cstemp, "%d", pnt.n_norms);
			stats += "     +n_norms: ";
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "%d", pnt.n_verts);
			stats += "     +n_verts: ";
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "%d", pnt.offset);
			stats += "     +offset: +";
			stats += cstemp;
			stats += "\n";

			Collector = 0;
			for (i = 0; i < pnt.n_verts; i++)
			{
			if (i < pnt.n_norms)
			Collector += int(unsigned char(pnt.norm_counts[i]));
			}

			memset(cstemp, 0, 64);
			sprintf(cstemp, "%d bytes", sizeof(vector) * (Collector + pnt.n_verts));
			stats += "     +vertex data size: ";
			stats += cstemp;
			stats += "\n";*/
			//stats += "\n";

			Add_DefPoints(pnt);
			break;

		case 2:
			//stats += "FLATPOLY  ";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "[%d bytes (%d bytes unused)]", head.size, head.size - fpol.Read(localptr, head));
			localptr += (head.size - 8);//8 for the already added header

			/*memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ Red: %d, Green: %d, Blue: %d, Pad: %d", fpol.red, fpol.green, fpol.blue, fpol.pad);
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ Normal: (%f, %f, %f), Radius %f", fpol.normal.x, fpol.normal.y, fpol.normal.z, fpol.radius);
			stats += cstemp;
			stats += "\n";

			for (i = 0; i < fpol.nverts; i++)
			{

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ Point: (%f, %f, %f),",
			points[0].vertex_data[fpol.verts[i].vertnum].vertex.x,
			points[0].vertex_data[fpol.verts[i].vertnum].vertex.y,
			points[0].vertex_data[fpol.verts[i].vertnum].vertex.z);
			stats += cstemp;

			memset(cstemp, 0, 64);
			sprintf(cstemp, "  Normal: (%f, %f, %f),",
			points[0].vertex_data[fpol.verts[i].vertnum].norms[fpol.verts[i].normnum].x,
			points[0].vertex_data[fpol.verts[i].vertnum].norms[fpol.verts[i].normnum].y,
			points[0].vertex_data[fpol.verts[i].vertnum].norms[fpol.verts[i].normnum].z);
			stats += cstemp;

			stats += "\n";


			}*/

			//stats += "\n";

			Add_FlatPoly(fpol);
			break;

		case 3:
			//	stats += "IMAPPOLY ";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "[%d bytes (%d bytes unused)]", head.size, head.size - tpol.Read(localptr, head));
			localptr += (head.size - 8);//8 for the already added header
			//stats += cstemp;
			//stats += "\n";



			/*memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ Normal: (%f, %f, %f), Radius %f", tpol.normal.x, tpol.normal.y, tpol.normal.z, tpol.radius);
			stats += cstemp;
			stats += "\n";

			for (i = 0; i < tpol.nverts; i++)
			{

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ Point: (%f, %f, %f),",
			points[0].vertex_data[tpol.verts[i].vertnum].vertex.x,
			points[0].vertex_data[tpol.verts[i].vertnum].vertex.y,
			points[0].vertex_data[tpol.verts[i].vertnum].vertex.z);
			stats += cstemp;

			memset(cstemp, 0, 64);
			sprintf(cstemp, "  Normal: (%f, %f, %f),",
			points[0].vertex_data[tpol.verts[i].vertnum].norms[tpol.verts[i].normnum].x,
			points[0].vertex_data[tpol.verts[i].vertnum].norms[tpol.verts[i].normnum].y,
			points[0].vertex_data[tpol.verts[i].vertnum].norms[tpol.verts[i].normnum].z);
			stats += cstemp;

			stats += "\n";


			}*/

			//stats += "\n";

			Add_TmapPoly(tpol);
			break;

		case 4:
			stats += "SORTNORM ";


			memset(cstemp, 0, 64);
			sprintf(cstemp, "[%d bytes (%d bytes unused)] ", head.size, head.size - snrm.Read(localptr, head));
			localptr += (head.size - 8);//8 for the already added header
			//stats += cstemp;
			//stats += "\n";

			/*memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ BackOffset: %d, FrontOffset: %d, OnlineOffset: %d, ",
			snrm.back_offset, snrm.front_offset, snrm.online_offset);
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ PostListOffset: %d, PreListOffset: %d, ",
			snrm.back_offset, snrm.front_offset, snrm.online_offset);
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ Plane_Normal: (%f, %f, %f)",
			snrm.plane_normal.x, snrm.plane_normal.y, snrm.plane_normal.z);
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ Plane_Point: (%f, %f, %f)",
			snrm.plane_point.x, snrm.plane_point.y, snrm.plane_point.z);
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ BoundBoxMin: (%f, %f, %f)",
			snrm.min_bounding_box_point.x, snrm.min_bounding_box_point.y, snrm.min_bounding_box_point.z);
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ BoundBoxMax: (%f, %f, %f)",
			snrm.max_bounding_box_point.x, snrm.max_bounding_box_point.y, snrm.max_bounding_box_point.z);
			stats += cstemp;
			stats += "\n";
			stats += "\n";*/

			Add_SortNorm(snrm);
			break;

		case 5:
			//stats += "BOUNDBOX ";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "[%d bytes (%d bytes unused)]", head.size, head.size - bnd.Read(localptr, head));
			localptr += (head.size - 8);//8 for the already added header
			/*stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ BoundBoxMin: (%f, %f, %f)", bnd.min_point.x, bnd.min_point.y, bnd.min_point.z);
			stats += cstemp;
			stats += "\n";

			memset(cstemp, 0, 64);
			sprintf(cstemp, "   @+ BoundBoxMax: (%f, %f, %f)", bnd.max_point.x, bnd.max_point.y, bnd.max_point.z);
			stats += cstemp;
			stats += "\n";
			stats += "\n";*/


			Add_BoundBox(bnd);
			break;

		default:
			stats += "WARNING: Encountered UNKNOWN: Check data format and/or pointer alignment!\n";
			go = false;
			break;

		}
	}

	return stats;
}

//********************************************************************************************************************


std::ostream& BSP::BSPDump(std::ostream &os) // dumps human readable BSP information into ostream;
{
	int i;
	for (i = 0; i < numpoints; i++)
	{
		os << "DEFPOINTS " << i << ":" << std::endl;
		os << "@N_VERTS: " << points[i].n_verts << std::endl;
		os << "@N_NORMS: " << points[i].n_norms << std::endl;
		os << "@OFFSET: " << points[i].offset << std::endl;

		for (int j = 0; j < points[i].n_verts; j++)
		{
			os << "  *Vertex[" << j << "] = " << points[i].vertex_data[j].vertex << std::endl;
			for (int k = 0; k < points[i].norm_counts[j]; k++)
				os << "       +Norm[" << k << "] = " << points[i].vertex_data[j].norms[k] << std::endl;
		}

	}

	os << " --------- POLYGONS --------- " << std::endl;

	for (i = 0; i < numtpolys; i++)
	{
		os << "@POLYGON: " << i << std::endl;
		os << "   *Norm: " << tpolys[i].normal << std::endl;
		os << "   *Center: " << tpolys[i].center << std::endl;
		os << "   *Radius: " << tpolys[i].radius << std::endl;
		os << "   *Texture: " << tpolys[i].tmap_num << std::endl;
		os << "   *NumVerts: " << tpolys[i].nverts << std::endl;
		for (int j = 0; j < tpolys[i].nverts; j++)
		{
			os << "       #vertnum: " << tpolys[i].verts[j].vertnum << std::endl;
			os << "       #normnum: " << tpolys[i].verts[j].normnum << std::endl;
			os << "       #u: " << tpolys[i].verts[j].u << std::endl;
			os << "       #v: " << tpolys[i].verts[j].v << std::endl;
		}
	}
	return os;
}

//********************************************************************************************************************

void BSP::Add_BoundBox(BSP_BoundBox bound)
{
	bounders.push_back(bound);

}

bool BSP::Del_BoundBox(int index)
{
	if (index < 0 || index >= numbounders)
		return false;
	bounders.erase(bounders.begin() + index);
	return true;
}

//********************************************************************************************************************


void BSP::Add_DefPoints(BSP_DefPoints pnts)
{
	points.push_back(pnts);
}

bool BSP::Del_DefPoints(int index)
{
	if (index < 0 || index >= numpoints)
		return false;
	points.erase(points.begin() + index);
	return true;
}

//********************************************************************************************************************

void BSP::Add_FlatPoly(BSP_FlatPoly fpol)
{
	fpolys.push_back(fpol);
}

bool BSP::Del_FlatPoly(int index)
{
	if (index < 0 || index >= numfpolys)
		return false;
	fpolys.erase(fpolys.begin() + index);
	return true;
}

//********************************************************************************************************************

void BSP::Add_SortNorm(BSP_SortNorm sn)
{
	snorms.push_back(sn);
}

bool BSP::Del_SortNorm(int index)
{
	if (index < 0 || index >= numsnorms)
		return false;
	snorms.erase(snorms.begin() + index);
	return true;
}

//********************************************************************************************************************

void BSP::Add_TmapPoly(BSP_TmapPoly tpol)
{
	tpolys.push_back(tpol);
}


bool BSP::Del_TmapPoly(int index)
{
	if (index < 0 || index >= numtpolys)
		return false;
	tpolys.erase(tpolys.begin() + index);
	return true;
}

std::ostream& operator<<(std::ostream &os, BSP_TmapPoly tpoly)
{
	os << "Textured Polygon" << std::endl;
	os << "@center: " << tpoly.center << std::endl;
	os << "@normal: " << tpoly.normal << std::endl;
	os << "@radius: " << tpoly.radius << std::endl;
	os << "@nverts: " << tpoly.nverts << std::endl;
	os << "@tmap_num: " << tpoly.tmap_num << std::endl;
	return os;
}

std::ostream& operator<<(std::ostream &os, BSP_FlatPoly fpoly)
{

	os << "NonTextured Polygon" << std::endl;
	os << "@center: " << fpoly.center << std::endl;
	os << "@normal: " << fpoly.normal << std::endl;
	os << "@radius: " << fpoly.radius << std::endl;
	os << "@nverts: " << fpoly.nverts << std::endl;
	os << "@RGB: " << int(fpoly.red) << ","
		<< int(fpoly.green) << "," << int(fpoly.blue) << std::endl;
	os << "@Pad: " << int(fpoly.pad) << std::endl;
	return os;
}
