#include "BSPDataStructs.h"
#include <memory.h>
#include <cmath>

//********************************************************************************************************************

int BSP_BlockHeader::Read(char *buffer) // Binary mode!
{
	memcpy(&id, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&size, buffer, sizeof(int));
	buffer += sizeof(int);

	return MySize();
}

int BSP_BlockHeader::Write(char *buffer)
{
	memcpy(buffer, &id, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &size, sizeof(int));
	buffer += sizeof(int);

	return 8;
}


//********************************************************************************************************************

int BSP_DefPoints::Read(char *buffer, BSP_BlockHeader hdr)
{
	char *temp = buffer;

	head = hdr;
	memcpy(&n_verts, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&n_norms, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&offset, buffer, sizeof(int));
	buffer += sizeof(int);

	norm_counts.resize(n_verts);
	memcpy(&norm_counts.front(), buffer, n_verts);
	buffer += n_verts;

	// ok this is based off my GUESS .. i hope i guessed correctly
	vertex_data.resize(n_verts);

	//realign the buffer
	buffer = temp;
	buffer += (offset - 8);
	// buffer now in alignment to read the vertex data;

	for (int i = 0; i < n_verts; i++)
	{
		memcpy(&vertex_data[i].vertex, buffer, sizeof(vector3d));
		buffer += sizeof(vector3d);

		if (int(norm_counts[i]) > 0)
		{
			vertex_data[i].norms.resize(norm_counts[i]);
			memcpy(&vertex_data[i].norms.front(), buffer, sizeof(vector3d) * int(norm_counts[i]));
			size_t offset = normals.size();
			normals.resize(normals.size() + norm_counts[i]);
			memcpy(&normals.front() + offset, buffer, sizeof(vector3d)* int(norm_counts[i]));
			buffer += (sizeof(vector3d) * norm_counts[i]);
		}
		else
			vertex_data[i].norms.clear();
	}
	//vertex_data

	return MySize();
}

int BSP_DefPoints::Write(char *buffer)
{
	int Collector = 0;
	char *tbuff = buffer;

	head.size = MySize();
	tbuff += head.Write(tbuff); //size of the header;

	memcpy(tbuff, &n_verts, sizeof(int));
	tbuff += sizeof(int);

	memcpy(tbuff, &n_norms, sizeof(int));
	tbuff += sizeof(int);

	memcpy(tbuff, &offset, sizeof(int));
	tbuff += sizeof(int);

	memcpy(tbuff, &norm_counts.front(), n_verts);
	tbuff += n_verts;

	// ok this is based off my GUESS .. i hope i guessed correctly
	for (int i = 0; i < n_verts; i++)
	{
		memcpy(tbuff, &vertex_data[i].vertex, sizeof(vector3d));
		tbuff += sizeof(vector3d);

		if (norm_counts[i])
		{
			memcpy(tbuff, &vertex_data[i].norms.front(), sizeof(vector3d) * ((int)(unsigned char)norm_counts[i]));
			tbuff += (sizeof(vector3d) * norm_counts[i]);
			Collector += int(norm_counts[i]);
		}
	}
	//vertex_data

	return (unsigned int)(tbuff - buffer);
}


int BSP_DefPoints::MySize()
{
	/*
	struct BSP_DefPoints // 1 - Defines Verticies
	{
	BSP_BlockHeader head;			//8			|8
	int n_verts;					//4			|12
	int n_norms;					//4			|16
	int offset;						//4			|20			|from start of chunk to vertex data
	unsigned char *norm_counts;		//n_norms	|20+n_norms	|n_verts in size;
	vertdata *vertex_data;			//Equation	|20+n_norms+Equation

	};*/
	int Collector = 0;

	for (int i = 0; i < n_verts; i++)
	{
		Collector += (int)((unsigned char)norm_counts[i]);
	}

	return (20 + n_verts + (sizeof(vector3d) * (Collector + n_verts)));
}

//********************************************************************************************************************

int Flat_vertex::Read(char *buffer)
{
	memcpy(&vertnum, buffer, sizeof(short));
	buffer += sizeof(short);

	memcpy(&normnum, buffer, sizeof(short));
	buffer += sizeof(short);

	return MySize();
}

int Flat_vertex::Write(char *buffer)
{
	memcpy(buffer, &vertnum, sizeof(short));
	buffer += sizeof(short);

	memcpy(buffer, &normnum, sizeof(short));
	buffer += sizeof(short);

	return 4;
}

//********************************************************************************************************************

int BSP_FlatPoly::Read(char *buffer, BSP_BlockHeader hdr)
{
	head = hdr;

	memcpy(&normal, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(&center, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(&radius, buffer, sizeof(float));
	buffer += sizeof(float);

	memcpy(&nverts, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&red, buffer, sizeof(byte));
	buffer += sizeof(byte);

	memcpy(&green, buffer, sizeof(byte));
	buffer += sizeof(byte);

	memcpy(&blue, buffer, sizeof(byte));
	buffer += sizeof(byte);

	memcpy(&pad, buffer, sizeof(byte));
	buffer += sizeof(byte);

	verts.resize(nverts);

	for (int i = 0; i < nverts; i++)
	{
		verts[i].Read(buffer);
		buffer += 4; // sizeof short*2 - which is the size of a Flat_vertex
	}

	return MySize();

}


vector3d BSP_FlatPoly::MyCenter(std::vector<vector3d> Verts)
{

	float TotalArea = 0, triarea;
	vector3d Centroid = MakeVector(0, 0, 0), midpoint;


	for (size_t i = 0; i < Verts.size(); i++)
	{
		midpoint = Verts[i] + Verts[i + 1] + Verts[i + 2];
		midpoint = midpoint / 3;

		// Area of Triangle defined by P1, P2, P3 = vector3d(crossProduct(p2-p1, p3-p1)).magnitude()/2

		triarea = Magnitude(CrossProduct(Verts[i + 1] - Verts[i],
			Verts[i + 2] - Verts[i])); // this needs to be area * 2
		midpoint = triarea*midpoint;
		TotalArea += triarea;
		Centroid += midpoint;

	}

	Centroid = float(1.0 / TotalArea) * Centroid;
	return Centroid;

}

float BSP_FlatPoly::MyRadius(vector3d center, std::vector<vector3d> Verts)
{
	float RetVal = 0;


	vector3d max;
	max.x = Abs(Verts[0].x);
	max.y = Abs(Verts[0].y);
	max.z = Abs(Verts[0].z);

	for (int i = 0; i < nverts; i++)
	{
		if (Abs(Verts[i].x) > max.x)
			max.x = Abs(Verts[i].x);

		if (Abs(Verts[i].y) > max.y)
			max.y = Abs(Verts[i].y);

		if (Abs(Verts[i].z) > max.z)
			max.z = Abs(Verts[i].z);
	}

	RetVal = ((max.x - Abs(center.x)) * (max.x - Abs(center.x))) +
		((max.y - Abs(center.y)) * (max.y - Abs(center.y))) +
		((max.z - Abs(center.z)) * (max.z - Abs(center.z)));
	RetVal = float(sqrt(RetVal));

	return RetVal;
}



int BSP_FlatPoly::Write(char *buffer)
{
	buffer += head.Write(buffer);

	memcpy(buffer, &normal, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(buffer, &center, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(buffer, &radius, sizeof(float));
	buffer += sizeof(float);

	memcpy(buffer, &nverts, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &red, sizeof(byte));
	buffer += sizeof(byte);

	memcpy(buffer, &green, sizeof(byte));
	buffer += sizeof(byte);

	memcpy(buffer, &blue, sizeof(byte));
	buffer += sizeof(byte);

	memcpy(buffer, &pad, sizeof(byte));
	buffer += sizeof(byte);

	for (int i = 0; i < nverts; i++)
	{
		buffer += verts[i].Write(buffer);
	}

	return (head.MySize() + (4 /*sizeof Flat_vertex*/ * nverts)
		+ (sizeof(byte) * 4) + sizeof(int) + sizeof(float) + (sizeof(vector3d) * 2));


}


//********************************************************************************************************************

int Tmap_vertex::Read(char *buffer)
{
	memcpy(&vertnum, buffer, sizeof(short));
	buffer += sizeof(short);

	memcpy(&normnum, buffer, sizeof(short));
	buffer += sizeof(short);

	memcpy(&u, buffer, sizeof(float));
	buffer += sizeof(float);

	memcpy(&v, buffer, sizeof(float));
	buffer += sizeof(float);

	return MySize();
}


int Tmap_vertex::Write(char *buffer)
{
	memcpy(buffer, &vertnum, sizeof(short));
	buffer += sizeof(short);

	memcpy(buffer, &normnum, sizeof(short));
	buffer += sizeof(short);

	memcpy(buffer, &u, sizeof(float));
	buffer += sizeof(float);

	memcpy(buffer, &v, sizeof(float));
	buffer += sizeof(float);

	return ((sizeof(short) * 2) + (sizeof(float) * 2));
}

//********************************************************************************************************************

int BSP_TmapPoly::Read(char *buffer, BSP_BlockHeader hdr)
{
	head = hdr;

	memcpy(&normal, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(&center, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(&radius, buffer, sizeof(float));
	buffer += sizeof(float);

	memcpy(&nverts, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&tmap_num, buffer, sizeof(int));
	buffer += sizeof(int);

	verts.resize(nverts);

	for (int i = 0; i < nverts; i++)
	{
		verts[i].Read(buffer);
		buffer += 12; //12 = sizeof(Tmap_vertex)
	}

	return MySize();

}

vector3d BSP_TmapPoly::MyCenter(std::vector<vector3d> Verts)
{

	float TotalArea = 0, triarea;
	vector3d Centroid = MakeVector(0, 0, 0), midpoint;


	for (int i = 0; i < nverts - 2; i++)
	{
		midpoint = Verts[i] + Verts[i + 1] + Verts[i + 2];
		midpoint = midpoint / 3;

		// Area of Triangle defined by P1, P2, P3 = vector3d(crossProduct(p2-p1, p3-p1)).magnitude()/2

		triarea = Magnitude(CrossProduct(Verts[i + 1] - Verts[i],
			Verts[i + 2] - Verts[i])); // this needs to be area * 2
		midpoint = triarea*midpoint;
		TotalArea += triarea;
		Centroid += midpoint;

	}

	Centroid = float(1.0 / TotalArea) * Centroid;
	return Centroid;
}

float BSP_TmapPoly::MyRadius(vector3d center, std::vector<vector3d> Verts)
{
	float RetVal = 0;


	vector3d max;
	max.x = Abs(Verts[0].x);
	max.y = Abs(Verts[0].y);
	max.z = Abs(Verts[0].z);

	for (int i = 0; i < nverts; i++)
	{
		if (Abs(Verts[i].x) > max.x)
			max.x = Abs(Verts[i].x);

		if (Abs(Verts[i].y) > max.y)
			max.y = Abs(Verts[i].y);

		if (Abs(Verts[i].z) > max.z)
			max.z = Abs(Verts[i].z);
	}

	RetVal = ((max.x - Abs(center.x)) * (max.x - Abs(center.x))) +
		((max.y - Abs(center.y)) * (max.y - Abs(center.y))) +
		((max.z - Abs(center.z)) * (max.z - Abs(center.z)));
	RetVal = float(sqrt(RetVal));

	return RetVal;
}


int BSP_TmapPoly::Write(char *buffer)
{
	buffer += head.Write(buffer);

	memcpy(buffer, &normal, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(buffer, &center, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(buffer, &radius, sizeof(float));
	buffer += sizeof(float);

	memcpy(buffer, &nverts, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &tmap_num, sizeof(int));
	buffer += sizeof(int);

	for (int i = 0; i < nverts; i++)
	{

		buffer += verts[i].Write(buffer); //12 = sizeof(Tmap_vertex)
	}

	return ((12 * nverts) + (sizeof(int) * 2) + sizeof(float) + (sizeof(vector3d) * 2));
}

//********************************************************************************************************************

int BSP_SortNorm::Read(char *buffer, BSP_BlockHeader hdr)
{
	head = hdr;

	memcpy(&plane_normal, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(&plane_point, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(&reserved, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&front_offset, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&back_offset, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&prelist_offset, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&postlist_offset, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&online_offset, buffer, sizeof(int));
	buffer += sizeof(int);

	memcpy(&min_bounding_box_point, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(&max_bounding_box_point, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	return MySize();
}

int BSP_SortNorm::Write(char *buffer)
{
	buffer += head.Write(buffer); // i was getting a all kinds of problems because i forgot to write this

	memcpy(buffer, &plane_normal, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(buffer, &plane_point, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(buffer, &reserved, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &front_offset, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &back_offset, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &prelist_offset, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &postlist_offset, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &online_offset, sizeof(int));
	buffer += sizeof(int);

	memcpy(buffer, &min_bounding_box_point, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(buffer, &max_bounding_box_point, sizeof(vector3d));
	buffer += sizeof(vector3d);

	return 72;
}


//********************************************************************************************************************

int BSP_BoundBox::Read(char *buffer, BSP_BlockHeader hdr)
{
	head = hdr;
	memcpy(&min_point, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(&max_point, buffer, sizeof(vector3d));
	buffer += sizeof(vector3d);

	return MySize();
}

int BSP_BoundBox::Write(char *buffer)
{
	buffer += head.Write(buffer);

	memcpy(buffer, &min_point, sizeof(vector3d));
	buffer += sizeof(vector3d);

	memcpy(buffer, &max_point, sizeof(vector3d));
	buffer += sizeof(vector3d);

	return 30;
}

//********************************************************************************************************************


bool operator==(BSP_TmapPoly &a, BSP_TmapPoly &b)
{
	bool Ret = ((a.center == b.center) && (a.radius == b.radius) && (a.normal == b.normal) && (a.nverts == b.nverts));

	if (Ret)
	{
		for (int i = 0; i < a.nverts; i++)
		{
			if (!(a.verts[i].normnum == b.verts[i].normnum &&
				a.verts[i].vertnum == b.verts[i].vertnum &&
				a.verts[i].u == b.verts[i].u &&
				a.verts[i].v == b.verts[i].v))
				return false;

		}
	}
	else return false;
	return true;
}



bool operator==(BSP_FlatPoly &a, BSP_FlatPoly &b)
{
	bool Ret = ((a.center == b.center) && (a.radius == b.radius) && (a.normal == b.normal) && (a.nverts == b.nverts));

	if (Ret)
	{
		for (int i = 0; i < a.nverts; i++)
		{
			if (!(a.verts[i].normnum == b.verts[i].normnum &&
				a.verts[i].vertnum == b.verts[i].vertnum))
				return false;

		}
	}
	else
		return false;
	return true;
}
