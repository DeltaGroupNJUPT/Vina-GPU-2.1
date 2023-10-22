/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_COMMON_H
#define VINA_COMMON_H

#include <cassert>
#include <string>
#include <limits>
#include <utility> // pair
#include <algorithm> // too common
#include <vector> // used in typedef, and commonly used overall
#include <cmath> // commonly used
#include <iostream> // various debugging everywhere
#include <fstream> // print_coords
#include <iomanip> // to_string
#include <sstream> // to_string
#include <string> // probably included by the above anyway, common anyway

#include <boost/serialization/vector.hpp> // can't come before the above two - wart fixed in upcoming Boost versions
#include <boost/serialization/base_object.hpp> // movable_atom needs it - (derived from atom)
#include <boost/filesystem/path.hpp> // typedef'ed

#include "macros.h"

typedef double fl;

template<typename T>					//C++模板,形参(不需要明确定义变量类型)
T sqr(T x) {
	return x*x;							//X^2
}

const fl not_a_num = std::sqrt(fl(-1)); // FIXME? check ????似乎是强制跳出?保护作用？

typedef std::size_t sz;
typedef std::pair<fl, fl> pr;

/*
vec结构体
data[0]~[2]
*/
struct vec {
	fl data[3];
	vec() {
#ifndef NDEBUG
		data[0] = data[1] = data[2] = not_a_num;
#endif
	}
	vec(fl x, fl y, fl z) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}

	//通过operator[]重构获取data[]的值
	const fl& operator[](sz i) const { assert(i < 3); return data[i]; }
	      fl& operator[](sz i)       { assert(i < 3); return data[i]; }

	//data数组各元素平方和
    fl norm_sqr() const {
		return sqr(data[0]) + sqr(data[1]) + sqr(data[2]);
	}

	//data数组各元素平方和，开根号(空间向量模)
	fl norm() const {
		return std::sqrt(norm_sqr());
	}

	/*
	function:结构体直接调用函数 实现两个数组的按位乘法后求和
	demo:vec a(1,2,3)
		 vec b(1,2,3)
		 c = a.operator*(b)		// c = 1*1+2*2+3*3=12
	*/
	fl operator*(const vec& v) const {
		return data[0] * v[0] + data[1] * v[1] + data[2] * v[2];
	}

	/*
	function:重构+=运算操作,vec1+=vec2
	demo:vec a, b, c;
		 a = { 1,2,3 };
		 b = { 2,2,2 };
		 c = a+=b;			//c.data={3,4,5}
	*/
	const vec& operator+=(const vec& v) {
		data[0] += v[0];
		data[1] += v[1];
		data[2] += v[2];
		return *this;
	}

	/*
	function:重构-=运算操作,vec1-=vec2
	demo:vec a, b, c;
		 a = { 1,2,3 };
		 b = { 3,3,3 };
		 c = a-=b;			//c.data={2,1,0}
	*/
	const vec& operator-=(const vec& v) {
		data[0] -= v[0];
		data[1] -= v[1];
		data[2] -= v[2];
		return *this;
	}

	/*
	function:重构+=运算操作, vec+=s
	demo:vec a,c;
		 double b;
		 a = { 1,2,3 };
		 b = 2.00
		 c = a+=b;			//c.data={3,4,5}
	*/
	const vec& operator+=(fl s) {
		data[0] += s;
		data[1] += s;
		data[2] += s;
		return *this;
	}

	/*
	function:重构-=运算操作, vec-=s
	demo:vec a,c;
		double b;
		a = { 1,2,3 };
		b = 1.00
		c = a+=b;			//c.data={0,1,2}
	*/
	const vec& operator-=(fl s) {
		data[0] -= s;
		data[1] -= s;
		data[2] -= s;
		return *this;
	}

	/*
	function:重构+运算操作,vec1+vec2
	demo:vec a, b, c;
	a = { 1,2,3 };
	b = { 2,2,2 };
	c = a+b;			//c.data={3,4,5}
	*/
	vec operator+(const vec& v) const {
		return vec(data[0] + v[0],
			       data[1] + v[1],
				   data[2] + v[2]);
	}

	/*
	function:重构-运算操作,vec1-vec2
	demo:vec a, b, c;
	a = { 2,2,3 };
	b = { 2,2,2 };
	c = a-b;			//c.data={0,0,1}
	*/
	vec operator-(const vec& v) const {
		return vec(data[0] - v[0],
			       data[1] - v[1],
				   data[2] - v[2]);
	}

	/*
	function:重构*=运算操作, vec*=s
	demo:vec a,c;
		double b;
		a = { 1,2,3 };
		b = 2.00
		c = a*=b;			//c.data={2,4,6}
	*/
	const vec& operator*=(fl s) {
		data[0] *= s;
		data[1] *= s;
		data[2] *= s;
		return *this;
	}
	/*
	将data数组均赋值为s
	*/
	void assign(fl s) {
		data[0] = data[1] = data[2] = s;
	}
	//固定范围size为3
	sz size() const { return 3; }

	//数据序列化
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & data;
	}
};

/*
function:s * v[0], s * v[1], s * v[2]
	     值依次存入vec结构体的data[0]~data[2]
*/
inline vec operator*(fl s, const vec& v) {
	return vec(s * v[0], s * v[1], s * v[2]);
}

/*
function:a[1]*b[2] - a[2]*b[1]
		 a[2]*b[0] - a[0]*b[2]
		 a[0]*b[1] - a[1]*b[0]
值依次存入vec结构体的data[0]~data[2]
*/
inline vec cross_product(const vec& a, const vec& b) {
	return vec(a[1]*b[2] - a[2]*b[1],
		       a[2]*b[0] - a[0]*b[2],
		       a[0]*b[1] - a[1]*b[0]);
}

/*
function:a[0] * b[0]
		 a[1] * b[1]
		 a[2] * b[2]
值依次存入vec结构体的data[0]~data[2]
*/
inline vec elementwise_product(const vec& a, const vec& b) {
	return vec(a[0] * b[0],
		       a[1] * b[1],
			   a[2] * b[2]);
}


/*
mat结构体
data[0]~[8]
*/
struct mat {
	fl data[9];
	mat() {
#ifndef NDEBUG
		data[0] = data[1] = data[2] =
		data[3] = data[4] = data[5] =
		data[6] = data[7] = data[8] = not_a_num;
#endif
	}
	// column-major
	const fl& operator()(sz i, sz j) const { assert(i < 3); assert(j < 3); return data[i + 3*j]; }
	      fl& operator()(sz i, sz j)       { assert(i < 3); assert(j < 3); return data[i + 3*j]; }

	mat(fl xx, fl xy, fl xz,
		fl yx, fl yy, fl yz,
		fl zx, fl zy, fl zz) {

		data[0] = xx; data[3] = xy; data[6] = xz;
		data[1] = yx; data[4] = yy; data[7] = yz;
		data[2] = zx; data[5] = zy; data[8] = zz;
	}

	/*
	function:结构体直接调用函数 实现两个数组的按位乘法后求和
	demo:mat a(1,2...,9)
	vec b(1,2,3)
	vec c = a.operator*(b)		// c.data[0] = 1*1+4*2+7*3
								   c.data[1] = 2*1+5*2+8*3
								   c.data[2] = 3*1+6*2+9*3
	*/
	vec operator*(const vec& v) const {
		return vec(data[0]*v[0] + data[3]*v[1] + data[6]*v[2], 
			       data[1]*v[0] + data[4]*v[1] + data[7]*v[2],
				   data[2]*v[0] + data[5]*v[1] + data[8]*v[2]);
	}

	/*
	function:重构*=运算操作, mat*=s
	demo:mat a,c;
	double b;
	a = { 1,2...,9 };
	b = 2.00
	c = a*=b;			//c.data={ 2,4,...,18}
	*/
	const mat& operator*=(fl s) {
		VINA_FOR(i, 9)
			data[i] *= s;
		return *this;
	}
};

//vec结构体型 vector
typedef std::vector<vec> vecv;

typedef std::pair<vec, vec> vecp;
//double型 vector
typedef std::vector<fl> flv;
//pair<int,int>型//类似map
typedef std::vector<pr> prv;
//unsigned int 型 vector
typedef std::vector<sz> szv;
typedef boost::filesystem::path path;

//???内部错误？周期误差？
struct internal_error {
	std::string file;
	unsigned line;
	internal_error(const std::string& file_, unsigned line_) : file(file_), line(line_) {}
};

#ifdef NDEBUG
	#define VINA_CHECK(P) do { if(!(P)) throw internal_error(__FILE__, __LINE__); } while(false)
#else
	#define VINA_CHECK(P) assert(P)
#endif

const fl pi = fl(3.1415926535897931);

/*
input: double x; uint max_sz
function: 若0<x<max_sz则 强制整形并输出
*/
inline sz fl_to_sz(fl x, sz max_sz) { // return a value in [0, max_sz]
    if(x <= 0) return 0;
    if(x >= max_sz) return max_sz;
	sz tmp = static_cast<sz>(x);		
	return (std::min)(tmp, max_sz); // sz -> fl cast loses precision. 'min' is to guard against returning values out of range
}

const fl fl_tolerance = fl(0.001);

/*
input: double a; double b
function: 计算|a-b|是否小于0.001返回布尔逻辑
*/
inline bool eq(fl a, fl b) {
	return std::abs(a - b) < fl_tolerance; 
}

/*
input: vec a; vec b
function: 计算a,b对应结构体中的data数组中
          对应位的差的绝对值否小于0.001返回布尔逻辑
*/
inline bool eq(const vec& a, const vec& b) {
	return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);
}

template<typename T>
/*
input: 可变类型的vector a和vector b
function:若a!=b的深度不等或其中元素存在不相等的情况
		 则返回false
*/
bool eq(const std::vector<T>& a, const std::vector<T>& b) {
	if(a.size() != b.size()) return false;
	VINA_FOR_IN(i, a)
		if(!eq(a[i], b[i]))
			return false;
	return true;
}


const fl max_fl = (std::numeric_limits<fl>::max)();							//当前处理器下double型变量最大值
const sz max_sz = (std::numeric_limits<sz>::max)();							//当前处理器下uint型变量最大值
const unsigned max_unsigned = (std::numeric_limits<unsigned>::max)();
const fl epsilon_fl = std::numeric_limits<fl>::epsilon();

const vec zero_vec(0, 0, 0);
const vec max_vec(max_fl, max_fl, max_fl);

/*
input: double x;
function: 判断x<0.1*max_fl(当前处理旗下double型变量最大值)
*/
inline bool not_max(fl x) {
	return (x < 0.1 * max_fl);
}


/*
input: 结构体vec a; vec b
function: 运算vec中逐个data数组元素差的平方
demo:vec a(4,4,4)
	 vec b(1,2,3)
	 c = vec_distance_sqr(a,b)
   //c = (4-1)^2+(4-2)^2+(4-3)^2=14 
*/
inline fl vec_distance_sqr(const vec& a, const vec& b) {
	//"\"为换行符
	return sqr(a[0] - b[0]) + \
		   sqr(a[1] - b[1]) + \
		   sqr(a[2] - b[2]);
}

/*
input: 结构体vec a
function: 运算vec中逐个data数组元素的平方和
demo:vec a(1,2,3)
	 b = sqr(a)
//   b = 1^2+2^2+3^2=14
*/
inline fl sqr(const vec& v) {
	return sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
}

/*
input: vector x; vector y
funtion:返回原有x向量深度值,并将y向量内容全部插入x向量末尾
demo: vector<int>x(2)=[1,2]   
	  vector<int>y(3)=[1,2,3]
	  unsigned a = vector_append(x,y)
	  a = 2
	  //执行后x = [1,2,1,2,3]
*/
template<typename T>							//C++模板,形参(不需要明确定义变量类型)
sz vector_append(std::vector<T>& x, const std::vector<T>& y) { // return old size
	sz old_size = x.size();
	x.insert(x.end(), y.begin(), y.end());
	return old_size;
}

/*
input: vector v
funtion:返回向量v中最小值坐标
demo: vector<int>v(4)=[2,2,1,3]
	  unsigned a = find_min(v)
	  //执行后a = 2
*/
template<typename T>
sz find_min(const std::vector<T>& v) { // returns v.size() i.e. 0 for empty vectors; the index of the smallest elt otherwise
	sz tmp = v.size();
	VINA_FOR_IN(i, v)
		if(i == 0 || v[i] < v[tmp])
			tmp = i;
	return tmp;
}

/*
input: double x
funtion:将x进行特征重放缩(标准化),范围为[-π,+π]
*/
inline void normalize_angle(fl& x) { // subtract or add enough 2*pi's to make x be in [-pi, pi]
	if     (x >  3*pi) { // very large
		fl n = ( x - pi) / (2*pi); // how many 2*pi's do you want to subtract?
		x -= 2*pi*std::ceil(n); // ceil can be very slow, but this should not be called often
		normalize_angle(x);
	}
	else if(x < -3*pi) { // very small
		fl n = (-x - pi) / (2*pi); // how many 2*pi's do you want to add?
		x += 2*pi*std::ceil(n); // ceil can be very slow, but this should not be called often
		normalize_angle(x);
	}
	else if(x >    pi) { // in (   pi, 3*pi]
		x -= 2*pi;
	}
	else if(x <   -pi) { // in [-3*pi,  -pi)
		x += 2*pi;
	}
	assert(x >= -pi && x <= pi);
	// in [-pi, pi]
}

/* 
input: double x 
funtion:将x进行特征重放缩(标准化),范围为[-π,+π]并返回放缩后的值
*/
inline fl normalized_angle(fl x) {
	normalize_angle(x);
	return x;
}

/*
input: 任意类型 x; long long width; char fill
funtion:输出width长度的字符串,
		若width<0则打印x的内容,
		否则打印fill+x的内容(fill的打印长度=width-x.size())
demo0: x = 222; width = -1; fill = '*';
	   to_string最终输出 "222"
demo1: x = 222; width = 2; fill = '*';
	   to_string最终输出 "222"
demo2: x = 222; width = 4; fill = '*';
	   to_string最终输出 "*222"
*/
template<typename T>
std::string to_string(const T& x, std::streamsize width = 0, char fill = ' ') { // default 0 width means no restrictions on width
	std::ostringstream out;
	out.fill(fill);
	if(width > 0)
		out << std::setw(width);
	out << x;
	return out.str();
}

/*
input: 任意类型vector向量 v
funtion:将v中所有元素累加求和并返回结果
demo: std::vector<unsigned int>a(3)=[1,2,3];
	  sum(a)返回结果为6
*/
template<typename T>
T sum(const std::vector<T>& v) {
	T acc = 0;
	VINA_FOR_IN(i, v)
		acc += v[i];
	return acc;
}

// multiply pK by this to get free energy in kcal/mol:
	// K = exp(E/RT)  -- lower K and lower E == better binder
	// pK = -log10(K)   => K = 10^(-pK)
	// E = RT ln(K) = RT ln (10^(-pK)) = - RT * ln(10) * pK
const fl pK_to_energy_factor = -8.31 /* RT in J/K/mol */ * 0.001 /* kilo */ * 300 /* K */ / 4.184 /* J/cal */ * std::log(10.0); //  -0.6 kcal/mol * log(10) = -1.38

/*
input: double pK
funtion:-8.31*pK
*/
inline fl pK_to_energy(fl pK) { return pK_to_energy_factor * pK; }

//打印double型x
inline void print(fl x, std::ostream& out = std::cout) {
	out << x;
}

//打印unsigned型x
inline void print(sz x, std::ostream& out = std::cout) {
	out << x;
}

/*
逐一打印vec v中内容,并以,分割
demo: vec a(1,2,3)
	  print(a)
	  // (1,2,3)
*/
inline void print(const vec& v, std::ostream& out = std::cout) {
	out << "(";
	VINA_FOR_IN(i, v) {
		if(i != 0) 
			out << ", ";
		print(v[i], out);
	}
	out << ")";
}


/* 
逐一打印vec v中内容,并以,分割 
demo: std::vector a(2,1)	  
	  print(a) 	  
	 // [1,1]
*/
template<typename T>
void print(const std::vector<T>& v, std::ostream& out = std::cout) {
	out << "[";
	VINA_FOR_IN(i, v) {
		if(i != 0) 
			out << " ";
		print(v[i], out);
	}
	out << "]";
}

/*
打印x变量
*/
template<typename T>
void printnl(const T& x, std::ostream& out = std::cout) {
	print(x, out);
	out << '\n';
}

/*
input: string str,start;
function: str字符长度大于等于start字符长度
		  且str.substr(0, start.size()) == start			//??? str不就等于str.substr(0, start.size())
		  则返回布尔值
*/
inline bool starts_with(const std::string& str, const std::string& start) {
	return str.size() >= start.size() && str.substr(0, start.size()) == start;
}

/*
input: vector v; 任意类型element
function: element是否存在于v中
*/
template<typename T>
bool has(const std::vector<T>& v, const T& element) {
	return std::find(v.begin(), v.end(), element) != v.end();
}

#endif
