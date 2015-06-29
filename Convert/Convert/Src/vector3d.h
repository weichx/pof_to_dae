#include <ios>
#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_
#include <string>
#include <vector>
#include <cstring>

// COB and POF have this in common
struct vector3d {
	float x, y, z;

	vector3d() : x(0.0), y(0.0), z(0.0) {}
	vector3d(float ax, float ay, float az) : x(ax), y(ay), z(az) {}
	/*vector3d(double ax, double ay, double az)
	{
	x = float(ax);
	y = float(ay);
	az = float(az);
	}*/
	//vector3d(vector3d &v) : x(v.x), y(v.y), z(v.z) {}
	//operator=(vector3d &v) { x=v.x; y=v.x; z=v.z; }
	float&operator[](int i){
		if (i % 3 == 0)return x;
		if (i % 3 == 1)return y;
		if (i % 3 == 2)return z;
		static float error;
		return error;//should not happen
	}
	float operator[](int i) const {
		if (i % 3 == 0)return x;
		if (i % 3 == 1)return y;
		if (i % 3 == 2)return z;
		static float error;
		return error;//should not happen
	}
};

//return true if non NAN
inline bool no_nan(float f){
	return f>0.0f || f <= 0.0f;
}
//return true if none of the values are NAN
inline bool no_nan(vector3d&vec){
	return no_nan(vec.x) && no_nan(vec.y) && no_nan(vec.z);
}

float VectorMagnitude(vector3d &v);

std::ostream & operator << (std::ostream & os, vector3d &vec);
std::istream & operator >> (std::istream & in, vector3d &vec);
bool operator==(vector3d &one, vector3d &two);
bool operator==(const vector3d &one, const vector3d &two);
bool operator!=(const vector3d &one, const vector3d &two);
bool operator!=(vector3d &one, vector3d &two);

vector3d operator-(const vector3d &one, const vector3d &two);
vector3d operator+(const vector3d &one, const vector3d &two);
vector3d& operator+=(vector3d &one, const vector3d &two);
vector3d& operator-=(vector3d &one, const vector3d &two);

vector3d MakeUnitVector(const vector3d &vect);
vector3d SafeMakeUnitVector(const vector3d &vect);

vector3d operator *(float scalar, const vector3d &v);
vector3d operator *(const vector3d &v, float scalar);
vector3d operator /(const vector3d &v, float scalar);
float Distance(const vector3d &one, const vector3d &two);

bool operator>(const vector3d &one, const vector3d &two);
bool operator<(const vector3d &one, const vector3d &two);

bool operator>=(const vector3d &one, const vector3d &two);
bool operator<=(const vector3d &one, const vector3d &two);

vector3d MakeVector(float ax, float ay, float az);
float Abs(float n);

vector3d AverageVectors(int numvectors, vector3d *vects);
float Angle(vector3d &v1, vector3d &v2);

vector3d CrossProduct(const vector3d &one, const vector3d &two);
float dot(const vector3d& A, const vector3d& B);
float Magnitude(const vector3d &vec);
inline bool null_vec(vector3d &a){ return a.x == 0.0f && a.y == 0.0f && a.z == 0.0f; }

vector3d average_vectors_if_less_than_angle(int numvectors, float angle, vector3d src, vector3d *vects);

vector3d FigureNormal(vector3d &one, vector3d &two, vector3d &three);
float FindObjectRadius(const vector3d &max, const vector3d &min, const vector3d &center);
void ExpandBoundingBoxes(vector3d &max, vector3d &min, const vector3d &cur);

//returns distance of p from the line defined as starting at lp and
//going in the direction of ln, is negitive if the point is behind
//the start of the line
float point_line_distance(const vector3d& p, const vector3d& lp, vector3d& ln);

//gives the point were the
//success is true unless the line is paralel to the plane,
//in this case the return value is {0,0,0}
vector3d plane_line_intersect(vector3d plane_pt, vector3d plane_norm, vector3d line_pt, vector3d line_norm, bool*success = NULL);

//returns the closest point on line 2 to line 1
vector3d closest_point_between_lines(vector3d p1, vector3d v1, vector3d p2, vector3d v2);

//returns the closes point on he line tp p
vector3d closest_point_on_line(vector3d p, vector3d lp, vector3d ln);

//given a point on the poly plane, find if it is in the poly
//stolen from FS2 source, which was in turn stolen from
//From Graphics Gems I, "An efficient Ray-Polygon intersection", p390
//returns non zero on sucsess
int point_face(vector3d *checkp, std::vector<vector3d> verts, vector3d * norm1);

namespace std {
	template<>
	struct hash<vector3d> 	{
		typedef vector3d argument_type;
		typedef std::size_t value_type;

		value_type operator()(argument_type const& v) const {
			value_type const h1(std::hash<float>()(v.x));
			value_type const h2(std::hash<float>()(v.y));
			value_type const h3(std::hash<float>()(v.z));
			return h1 ^ (h2 << 7) ^ (h3 << 13);
		}
	};
}


//there is probly a better place for this, but as of now it's only used in the header
namespace bobboau {

	struct matrix{
		matrix(){
			memset(a2d, 0, sizeof(float) * 9);
			a2d[0][0] = 1;
			a2d[1][1] = 1;
			a2d[2][2] = 1;
		}
		matrix(float f[3][3]){
			memcpy(a2d, f, sizeof(float) * 9);
		}
		matrix(vector3d r, vector3d u, vector3d f){
			memcpy(a2d, &r, sizeof(float) * 3);
			memcpy(a2d[1], &u, sizeof(float) * 3);
			memcpy(a2d[2], &f, sizeof(float) * 3);
		}
		matrix(const matrix&m){
			memcpy(a2d, &m, sizeof(float) * 9);
		}
		matrix operator=(const matrix&m){
			memcpy(a2d, &m, sizeof(float) * 9);
			return (*this);
		}

		float a2d[3][3];

		matrix operator * (const float&f){
			return matrix(vector3d(a2d[0][0] * f, a2d[0][1] * f, a2d[0][2] * f),
				vector3d(a2d[1][0] * f, a2d[1][1] * f, a2d[1][2] * f),
				vector3d(a2d[2][0] * f, a2d[2][1] * f, a2d[2][2] * f)
				);
		}

		matrix operator + (const matrix&m){
			return matrix(vector3d(a2d[0][0] + m.a2d[0][0], a2d[0][1] + m.a2d[0][1], a2d[0][2] + m.a2d[0][2]),
				vector3d(a2d[1][0] + m.a2d[1][0], a2d[1][1] + m.a2d[1][1], a2d[1][2] + m.a2d[1][2]),
				vector3d(a2d[2][0] + m.a2d[2][0], a2d[2][1] + m.a2d[2][1], a2d[2][2] + m.a2d[2][2])
				);
		}
		matrix invert(){
			matrix ret;
			double d = -a2d[0][2] * a2d[1][1] * a2d[2][0] + a2d[0][1] * a2d[1][2] * a2d[2][0] + a2d[0][2] * a2d[1][0] * a2d[2][1] - a2d[0][0] * a2d[1][2] * a2d[2][1] - a2d[0][1] * a2d[1][0] * a2d[2][2] + a2d[0][0] * a2d[1][1] * a2d[2][2];

			ret.a2d[0][0] = float(double(-a2d[1][2] * a2d[2][1] + a2d[1][1] * a2d[2][2]) / d);
			ret.a2d[0][1] = float(double(a2d[0][2] * a2d[2][1] - a2d[0][1] * a2d[2][2]) / d);
			ret.a2d[0][2] = float(double(-a2d[0][2] * a2d[1][1] + a2d[0][1] * a2d[1][2]) / d);

			ret.a2d[1][0] = float(double(a2d[1][2] * a2d[2][0] - a2d[1][0] * a2d[2][2]) / d);
			ret.a2d[1][1] = float(double(-a2d[0][2] * a2d[2][0] + a2d[0][0] * a2d[2][2]) / d);
			ret.a2d[1][2] = float(double(a2d[0][2] * a2d[1][0] - a2d[0][0] * a2d[1][2]) / d);

			ret.a2d[2][0] = float(double(-a2d[1][1] * a2d[2][0] + a2d[1][0] * a2d[2][1]) / d);
			ret.a2d[2][1] = float(double(a2d[0][1] * a2d[2][0] - a2d[0][0] * a2d[2][1]) / d);
			ret.a2d[2][2] = float(double(-a2d[0][1] * a2d[1][0] + a2d[0][0] * a2d[1][1]) / d);

			return ret;
		}
	};

} //namespace bobboau

#endif
