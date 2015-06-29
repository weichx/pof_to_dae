#include "vector3d.h"
#include <cmath>
#include <memory.h>
#include <cstdio>
#include <iostream>


//returns the average of all vectors that are within the given angle of the src vector
vector3d average_vectors_if_less_than_angle(int numvectors, float angle, vector3d src, vector3d *vects)
{
	vector3d RetVal = MakeVector(0, 0, 0);
	int nv = 0;

	for (int i = 0; i < numvectors; i++)
	{
		if (Angle(src, vects[i])<angle){
			RetVal += vects[i];
			nv++;
		}
	}

	RetVal = RetVal / float(nv);

	return RetVal;
}


//****************************************************************************************************************


float Magnitude(const vector3d &vec)
{
	return (float)sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

//****************************************************************************************************************

vector3d CrossProduct(const vector3d &one, const vector3d &two)
{
	return MakeVector(float(double(one.y)*double(two.z) - double(two.y)*double(one.z)),
		float(double(two.x)*double(one.z) - double(one.x)*double(two.z)),
		float(double(one.x)*double(two.y) - double(two.x)*double(one.y)));
}

float dot(const vector3d& A, const vector3d& B){
	return A.x*B.x + A.y*B.y + A.z*B.z;
}

vector3d FigureNormal(vector3d &one, vector3d &two, vector3d &three)
{
	vector3d ret;
	double v1x, v1y, v1z, v2x, v2y, v2z;

	/*
	C.x = A.y*B.z - A.z*B.y;
	C.y = A.z*B.x - A.x*B.z;
	C.z = A.x*B.y - A.y*B.x;

	I found Cross some information on the net on how to do this.. he crossed some things
	*/

	v1x = double(one.x) - double(two.x);
	v1y = double(one.y) - double(two.y);
	v1z = double(one.z) - double(two.z);

	v2x = double(three.x) - double(two.x);
	v2y = double(three.y) - double(two.y);
	v2z = double(three.z) - double(two.z);

	ret.x = float(((v1y * v2z) - (v1z * v2y)));
	ret.y = float(((v1z * v2x) - (v1x * v2z)));
	ret.z = float(((v1x * v2y) - (v1y * v2x)));

	return ret;

}
//****************************************************************************************************************


float VectorMagnitude(vector3d &v)
{
	return Magnitude(v);
}

//****************************************************************************************************************


float Angle(vector3d &v1, vector3d &v2)
{
	float Dp1 = VectorMagnitude(v1)*VectorMagnitude(v2); // missing the angle calc
	float Dp2 = (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);

	float Ang = Dp2 / Dp1;
	if (Ang > 1.00000f)
		Ang = 1.00000f;
	if (Ang < -1.00000f)
		Ang = -1.00000f;
	Ang = (float)acos(Ang); //inverse cosine ang (cos-1 aka arccos)

	return Ang;
}


//****************************************************************************************************************

vector3d& operator+=(vector3d &one, const vector3d &two)
{
	one = one + two;
	return one;
}

vector3d& operator-=(vector3d &one, const vector3d &two)
{
	one = one - two;
	return one;
}


//****************************************************************************************************************

vector3d AverageVectors(int numvectors, vector3d *vects)
{
	vector3d RetVal = MakeVector(0, 0, 0);

	for (int i = 0; i < numvectors; i++)
	{
		RetVal += vects[i];
	}

	RetVal = RetVal / float(numvectors);

	return RetVal;
}

//****************************************************************************************************************

vector3d MakeVector(float ax, float ay, float az)
{
	vector3d temp;
	temp.x = ax;
	temp.y = ay;
	temp.z = az;

	return temp;
}

//****************************************************************************************************************

vector3d MakeUnitVector(const vector3d &vect)
{
	double VMag = sqrt(double(vect.x * vect.x)
		+ double(vect.y * vect.y)
		+ double(vect.z * vect.z));

	return vector3d(float(double(vect.x) / VMag),
		float(double(vect.y) / VMag),
		float(double(vect.z) / VMag));


}
//****************************************************************************************************************

vector3d SafeMakeUnitVector(const vector3d &vect)
{
	double VMag = sqrt(double(vect.x * vect.x)
		+ double(vect.y * vect.y)
		+ double(vect.z * vect.z));

	if (VMag < 1e-5) {
		return vector3d();
	}
	return vector3d(float(double(vect.x) / VMag),
		float(double(vect.y) / VMag),
		float(double(vect.z) / VMag));


}
//****************************************************************************************************************


bool operator>(const vector3d &one, const vector3d &two)
{
	return ((one.x >= two.x) && (one.y >= two.y) && (one.z >= two.z) && !(one == two));
}

bool operator<(const vector3d &one, const vector3d &two)
{
	return ((one.x <= two.x) && (one.y <= two.y) && (one.z <= two.z) && !(one == two));
}

bool operator>=(const vector3d &one, const vector3d &two)
{
	return ((one.x >= two.x) && (one.y >= two.y) && (one.z >= two.z));
}


bool operator<=(const vector3d &one, const vector3d &two)
{
	return ((one.x <= two.x) && (one.y <= two.y) && (one.z <= two.z));
}

bool operator==(vector3d &one, vector3d &two)
{
	return one.x == two.x && one.y == two.y && one.z == two.z;
}

bool operator==(const vector3d &one, const vector3d &two)
{
	return one.x == two.x && one.y == two.y && one.z == two.z;
}

bool operator!=(const vector3d &one, const vector3d &two)
{
	return !(one == two);
}

bool operator!=(vector3d &one, vector3d &two)
{
	return !(one == two);
}

//****************************************************************************************************************


vector3d operator+(const vector3d &one, const vector3d &two)
{
	vector3d RetVal;
	RetVal.x = one.x + two.x;
	RetVal.y = one.y + two.y;
	RetVal.z = one.z + two.z;
	return RetVal;
}


vector3d operator-(const vector3d &one, const vector3d &two)
{
	vector3d Vtemp;
	Vtemp.x = one.x - two.x;
	Vtemp.y = one.y - two.y;
	Vtemp.z = one.z - two.z;

	return Vtemp;
}

//****************************************************************************************************************


float Distance(const vector3d &one, const vector3d &two)
{
	float A = (one.x - two.x);
	float B = (one.y - two.y);
	float C = (one.z - two.z);
	float D = (A*A) + (B*B) + (C*C);

	return float(sqrt(D));
}


//****************************************************************************************************************


std::ostream & operator << (std::ostream & os, vector3d &vec)
{
	os << vec.x << " " << vec.y << " " << vec.z;
	return os;
}


std::istream & operator >> (std::istream & in, vector3d &vec)
{
	in >> vec.x;
	in >> vec.y;
	in >> vec.z;

	return in;
}


//****************************************************************************************************************

vector3d operator *(float scalar, const vector3d &v)
{
	vector3d Ret;
	Ret.x = scalar * v.x;
	Ret.y = scalar * v.y;
	Ret.z = scalar * v.z;
	return Ret;
}

//diferent parameter order
vector3d operator *(const vector3d &v, float scalar)
{
	vector3d Ret;
	Ret.x = scalar * v.x;
	Ret.y = scalar * v.y;
	Ret.z = scalar * v.z;
	return Ret;
}


vector3d operator /(const vector3d &v, float scalar)
{
	vector3d Ret;
	Ret.x = v.x / scalar;
	Ret.y = v.y / scalar;
	Ret.z = v.z / scalar;
	return Ret;
}


float Abs(float n)
{
	if (n < 0)
		n = -n;
	return n;
}

void ExpandBoundingBoxes(vector3d &max, vector3d &min, const vector3d &cur)
{
	if (cur.x > max.x)
		max.x = cur.x;
	if (cur.y > max.y)
		max.y = cur.y;
	if (cur.z > max.z)
		max.z = cur.z;

	if (cur.x < min.x)
		min.x = cur.x;
	if (cur.y < min.y)
		min.y = cur.y;
	if (cur.z < min.z)
		min.z = cur.z;

}

float FindObjectRadius(const vector3d &max, const vector3d &min, const vector3d &center)
{
	vector3d temp_vector;
	// Set Radius
	if (fabs(max.x - center.x) >
		fabs(center.x - min.x))
		temp_vector.x = fabs(max.x - center.x);
	else
		temp_vector.x = fabs(center.x - min.x);

	if (fabs(max.y - center.y) >
		fabs(center.y - min.y))
		temp_vector.y = fabs(max.y - center.y);
	else
		temp_vector.y = fabs(center.y - min.y);

	if (fabs(max.z - center.z) >
		fabs(center.z - min.z))
		temp_vector.z = fabs(max.z - center.z);
	else
		temp_vector.z = fabs(center.z - min.z);
	return Magnitude(temp_vector);
}



//returns distance of p from the line defined as starting at lp and 
//going in the direction of ln, is negitive if the point is behind 
//the start of the line
float point_line_distance(const vector3d& p, const vector3d& lp, vector3d& ln)
{

	ln = MakeUnitVector(ln);

	vector3d nearest;
	float comp;

	comp = dot(p - lp, ln);
	nearest = lp + ln*comp;

	if (comp < 0.0f)
		return -1.0f*Distance(nearest, p);
	else
		return Distance(nearest, p);

}


vector3d plane_line_intersect(vector3d plane_pt, vector3d plane_norm, vector3d line_pt, vector3d line_norm, bool*success){

	vector3d w = line_pt - plane_pt;

	float d = -dot(plane_norm, line_norm);

	if (d == 0.0f){
		if (success)*success = false;
		return vector3d(0, 0, 0);
	}
	if (success)*success = true;
	vector3d temp(line_norm);
	return line_pt + temp * (dot(plane_norm, w) / d);
}

//returns the closest point on the given line
vector3d closest_point_on_line(vector3d p, vector3d lp, vector3d ln){
	ln = MakeUnitVector(ln);

	return lp + ln*dot(p - lp, ln);
}


float matrix_determinant_from_vectors(vector3d v1, vector3d v2, vector3d v3)
{
	float ans;
	ans = v1.x * v2.y * v3.z;
	ans += v2.x * v3.y * v1.z;
	ans += v3.x * v1.y * v2.z;
	ans -= v1.z * v2.y * v3.x;
	ans -= v2.z * v3.y * v1.x;
	ans -= v3.z * v1.y * v2.x;

	return ans;
}

vector3d closest_point_between_lines(vector3d p1, vector3d v1, vector3d p2, vector3d v2)
{
	vector3d cross, delta;
	cross = CrossProduct(v1, v2);
	delta = p2 - p1;

	float denominator = Magnitude(cross);
	denominator *= denominator;
	float num_s, s;

	if (denominator > 1e-10) {
		num_s = matrix_determinant_from_vectors(delta, v2, cross);
		s = num_s / denominator;
		return p1 + v1*s;
	}
	else{
		return vector3d(0, 0, 0);
	}

}


#define delta 0.0001f
#define	UNINITIALIZED_VALUE	-1234567.8f
int ij_table[3][2] = {
	{ 2, 1 },          //pos x biggest
	{ 0, 2 },          //pos y biggest
	{ 1, 0 },          //pos z biggest
};

int point_face(vector3d *checkp, std::vector<vector3d> verts, vector3d * norm1)
{
	float *norm, *P;
	vector3d t;
	int i0, i1, i2;

	norm = (float *)norm1;

	//project polygon onto plane by finding largest component of normal
	t.x = fabsf(norm[0]);
	t.y = fabsf(norm[1]);
	t.z = fabsf(norm[2]);

	if (t.x > t.y) if (t.x > t.z) i0 = 0; else i0 = 2;
	else if (t.y > t.z) i0 = 1; else i0 = 2;

	if (norm[i0] > 0.0f) {
		i1 = ij_table[i0][0];
		i2 = ij_table[i0][1];
	}
	else {
		i1 = ij_table[i0][1];
		i2 = ij_table[i0][0];
	}

	// Have i0, i1, i2
	P = (float *)checkp;

	float u0, u1, u2, v0, v1, v2, alpha = UNINITIALIZED_VALUE;
	double beta;

	int inter = 0, i = 2;

	u0 = P[i1] - verts[0][i1];
	v0 = P[i2] - verts[0][i2];

	do {

		u1 = verts[i - 1][i1] - verts[0][i1];
		u2 = verts[i][i1] - verts[0][i1];

		v1 = verts[i - 1][i2] - verts[0][i2];
		v2 = verts[i][i2] - verts[0][i2];


		if ((u1 >-delta) && (u1<delta))	{
			beta = u0 / u2;
			if ((beta >= 0.0f) && (beta <= 1.0f))	{
				alpha = float((v0 - beta*v2) / v1);
				inter = ((alpha >= 0.0f) && (alpha + beta <= 1.0f));
			}
		}
		else {

			beta = (v0*u1 - u0*v1) / (v2*u1 - u2*v1);
			if ((beta >= 0.0f) && (beta <= 1.0f))	{
				alpha = float((u0 - beta*u2) / u1);
				inter = ((alpha >= 0.0f) && (alpha + beta <= 1.0f));
			}


		}

	} while ((!inter) && ((unsigned int)(++i) < verts.size()));

	return inter;
}
