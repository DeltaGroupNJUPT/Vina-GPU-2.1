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

#ifndef VINA_ATOM_H
#define VINA_ATOM_H

#include "atom_base.h"

/*
input: uint i; bool in_grid;
function:输入一个参数和对应布尔逻辑，
		 输出原子索引排序
demo: atom_index(1,ture);
*/
struct atom_index {
	sz i;
	bool in_grid;
	atom_index() : i(max_sz), in_grid(false) {}							//atom_index初始化
	atom_index(sz i_, bool in_grid_) : i(i_), in_grid(in_grid_) {}		//结构体重构
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & i;
		ar & in_grid;
	}
};

/*
input:atom_index_i; atom_index_j
function:判断两个atom_index是否完全相等，并输出布尔逻辑的结果
*/
inline bool operator==(const atom_index& i, const atom_index& j) {
	return i.i == j.i && i.in_grid == j.in_grid;
}

/*
键
input: struct connected_atom_index; double length; bool rotatable
function: 输入原子索引，长度，是否可旋转，将其序列化
*/
struct bond {
	atom_index connected_atom_index;
	fl length;
	bool rotatable;
	bond() : length(0), rotatable(false) {}
	bond(const atom_index& connected_atom_index_, fl length_, bool rotatable_) : connected_atom_index(connected_atom_index_), length(length_), rotatable(rotatable_) {}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & connected_atom_index;
		ar & length;
		ar & rotatable;
	}
};

/*
function:将atom_base, coords(vec)， bonds(bond)序列化
*/
struct atom : public atom_base {	//atom实质继承了atom_type
    vec coords;						//coords重写 vec结构体			
	std::vector<bond> bonds;
	atom() : coords(max_vec) {}		//通过max_vec对coords(vec)初始化
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & boost::serialization::base_object<atom_base>(*this);
		ar & coords;
		ar & bonds;
	}
};

typedef std::vector<atom> atomv;

#endif
