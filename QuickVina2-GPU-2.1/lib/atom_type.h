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

#ifndef VINA_ATOM_TYPE_H
#define VINA_ATOM_TYPE_H

#include "atom_constants.h"
#include "triangular_matrix_index.h"

/*
结构体
*/
struct atom_type {
	enum t {EL, AD, XS, SY};					//枚举型变量t 
	sz el, ad, xs, sy;							//uint变量

	/*初始化 el=11; ad=20; xs=17; sy=18*/
	atom_type() : el(EL_TYPE_SIZE), ad(AD_TYPE_SIZE), xs(XS_TYPE_SIZE), sy(SY_TYPE_SIZE) {}

	//枚举型t对应返回值
	sz get(t atom_typing_used) const {
		switch(atom_typing_used) {
			case EL: return el;
			case AD: return ad;
			case XS: return xs;
			case SY: return sy;
			default: assert(false); return max_sz;
		}
	}

	//是否为氢(ad==6||ad==12);
	bool is_hydrogen() const {
		return ad_is_hydrogen(ad);
	}

	//是否为杂原子(ad != 1 && ad != 0  && ad != 6 && ad != 12 && ad < 20;)或xs(X-Score)==16;
	bool is_heteroatom() const {
		return ad_is_heteroatom(ad) || xs == XS_TYPE_Met_D;
	}

	//是否为可接受类型(ad<20 || xs==16)
	bool acceptable_type() const {
		return ad < AD_TYPE_SIZE || xs == XS_TYPE_Met_D;
	}

	/*
	function: el赋值函数
	若ad==20且xs==16 则 el=10
	否则el=ad输入后的返回值
	*/
	void assign_el() {
		el = ad_type_to_el_type(ad);					//等于ad输入后对应返回值
		if(ad == AD_TYPE_SIZE && xs == XS_TYPE_Met_D)	
			el = EL_TYPE_Met;
	}

	//el是否为同一元素
	bool same_element(const atom_type& a) const { // does not distinguish metals or unassigned types
		return el == a.el;
	}

	/*
	double型 共价半径
	function:若ad<20则返回;
			 若xs==16则返回1.75
	*/
	fl covalent_radius() const {
		if(ad < AD_TYPE_SIZE)        return ad_type_property(ad).covalent_radius;
		else if(xs == XS_TYPE_Met_D) return metal_covalent_radius;
		VINA_CHECK(false);           return 0; // never happens - placating the compiler
	}
	/*
	double型 最佳共价键长度
	function: 计算出 共价半径+x.共价半径
	*/
	fl optimal_covalent_bond_length(const atom_type& x) const {
		return covalent_radius() + x.covalent_radius();
	}

	//参数序列化
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & el;
		ar & ad;
		ar & xs;
		ar & sy;
	}
};

/*
原子数类型
EL型返回11;
AD型返回20;
XS型返回17;
SY型返回18
*/
inline sz num_atom_types(atom_type::t atom_typing_used) {
	switch(atom_typing_used) {
		case atom_type::EL: return EL_TYPE_SIZE;
		case atom_type::AD: return AD_TYPE_SIZE;
		case atom_type::XS: return XS_TYPE_SIZE;
		case atom_type::SY: return SY_TYPE_SIZE;
		default: assert(false); return max_sz;
	}
}

/*
获取类型对索引
function:返回三角矩阵坐标
		 若i对应结构体的枚举型值<=j 则返回i + j*(j+1)/2
		 否则返回j + i*(i+1)/2
*/
inline sz get_type_pair_index(atom_type::t atom_typing_used, const atom_type& a, const atom_type& b) { // throws error if any arg is unassigned in the given typing scheme
	sz n = num_atom_types(atom_typing_used);			//n=原子数类型返回数据

	sz i = a.get(atom_typing_used); VINA_CHECK(i < n);	//i=a对应枚举型返回数据el, ad, xs, sy
	sz j = b.get(atom_typing_used); VINA_CHECK(j < n);	//j=b对应枚举型返回数据el, ad, xs, sy

	if(i <= j) return triangular_matrix_index(n, i, j);
	else       return triangular_matrix_index(n, j, i);
}

#endif
