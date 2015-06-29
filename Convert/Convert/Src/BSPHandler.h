#include "BSPDataStructs.h"
#include <ios>


#if !defined(_BSP_HANDLER_H_)
#define _BSP_HANDLER_H_


class BSP
{
public:
	std::vector<BSP_BoundBox> bounders;
	std::vector<BSP_DefPoints> points;
	std::vector<BSP_FlatPoly> fpolys;
	std::vector<BSP_SortNorm> snorms;
	std::vector<BSP_TmapPoly> tpolys;

	int numbounders, numpoints, numfpolys, numsnorms, numtpolys;
	void Clear()
	{
		bounders.clear();
		points.clear();
		fpolys.clear();
		snorms.clear();
		tpolys.clear();
	}

	BSP()
	{
		Clear();
	}
	BSP(char *buffer, int size)
	{
		Clear();
		DataIn(buffer, size);
	}
	//~BSP();

	//--------------------------------
	std::string DataIn(char *buffer, int size);
	std::ostream& BSPDump(std::ostream &os); // dumps human readable BSP information into ostream;

	int Count_Bounding()
	{
		return bounders.size();
	}

	int Count_Points()
	{
		return points.size();
	}

	int Count_FlatPolys()
	{
		return fpolys.size();
	}

	int Count_SortNorms()
	{
		return snorms.size();
	}

	int Count_TmapPolys()
	{
		return tpolys.size();
	}

	//--------------------------------

	void Add_BoundBox(BSP_BoundBox bound);
	bool Del_BoundBox(int index);

	void Add_DefPoints(BSP_DefPoints pnts);
	bool Del_DefPoints(int index);

	void Add_FlatPoly(BSP_FlatPoly fpol);
	bool Del_FlatPoly(int index);

	void Add_SortNorm(BSP_SortNorm sn);
	bool Del_SortNorm(int index);

	void Add_TmapPoly(BSP_TmapPoly tpol);
	bool Del_TmapPoly(int index);

};


std::ostream& operator<<(std::ostream &os, BSP_TmapPoly tpoly);
std::ostream& operator<<(std::ostream &os, BSP_FlatPoly fpoly);
bool operator==(BSP_TmapPoly &a, BSP_TmapPoly &b);
bool operator==(BSP_FlatPoly &a, BSP_FlatPoly &b);

#endif //_BSP_HANDLER_H_
