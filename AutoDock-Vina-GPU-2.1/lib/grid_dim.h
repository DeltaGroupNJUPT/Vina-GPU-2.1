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

#ifndef VINA_GRID_DIM_H
#define VINA_GRID_DIM_H

#include <boost/array.hpp>

#include "common.h"


/*结构体：
* double： begin、end    Uint：n
* function1：参数 double：begin、end    Uint：n  赋值0
* function2：返回 end- begin的值  double  只读
* function3： 判断n>0    bool  只读
* function4: 输入：两个结构体 输出：bool  判断两个结构体的n是否一样 参数begin和end的差的绝对值是否<0.01
* function5：输入：两个结构体组（三个结构体组成） 输出：bool  判断结构体组中的三对数组中的a[i]和b[i]的差的绝对值是否<0.001
* function6：输入：结构体组（三个结构体组成）  std::ostream类型引用，默认值为std::cout   功能：打印
* function7：输入：结构体组（三个结构体组成）  输出 ：长度3数组 double  将三个参数begin存入数组
* function8：输入：结构体组（三个结构体组成）   输出 ：长度3数组 double   将三个参数end存入数组
* 侵入式序列化     private部分
*/
struct grid_dim {
	fl begin;
	fl end;
	sz n; // number of intervals == number of sample points - 1
	grid_dim() : begin(0), end(0), n(0) {}
	fl span() const { return end - begin; }
	bool enabled() const { return (n > 0); }
private:                                           //侵入式序列化
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & begin;
		ar & end;
		ar & n;
	}
};

inline bool eq(const grid_dim& a, const grid_dim& b) {
	return a.n == b.n && eq(a.begin, b.begin) && eq(a.end, b.end);//判断两个结构体的n是否一样 参数begin和end的差的绝对值是否<0.01
}

typedef boost::array<grid_dim, 3> grid_dims;

inline bool eq(const grid_dims& a, const grid_dims& b) {
	return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);//判断结构体组中的三对数组中的a[i]和b[i]的差的绝对值是否<0.001
}       

inline void print(const grid_dims& gd, std::ostream& out = std::cout) {
	VINA_FOR_IN(i, gd)
		std::cout << gd[i].n << " [" << gd[i].begin << " .. " << gd[i].end << "]\n";
}//    打印n[b..e]   

inline vec grid_dims_begin(const grid_dims& gd) {
	vec tmp;
	VINA_FOR_IN(i, gd)
		tmp[i] = gd[i].begin;
	return tmp;
}//将三个begin存入数组

inline vec grid_dims_end(const grid_dims& gd) {
	vec tmp;
	VINA_FOR_IN(i, gd)
		tmp[i] = gd[i].end;
	return tmp;
}//将三个end存入数组

#endif
