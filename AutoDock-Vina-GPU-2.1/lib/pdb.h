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

#ifndef VINA_PDB_H
#define VINA_PDB_H

#include "common.h"


/*
*定义结构体pdb_atom，数据成员包括：
* 无符号类型整数id
* string类型name，residue_name，element
* int类型residue_id
* 结构体vec（矢量）类型coords，该结构里有一个含三个double类型数据的数组
* double类型b_factor
*/
struct pdb_atom {
	unsigned id;
	std::string name;
	std::string residue_name;
	int residue_id;
	vec coords;
	fl b_factor;
	std::string element;
};


/*
* 结构体pdb，数据成员包括：
* 容器atoms，该容器包含若干结构体pdb_atom
* 成员函数check()声明
*/
struct pdb {
	std::vector<pdb_atom> atoms;
	void check(fl min_distance) const;
};


//函数parse_pdb声明
pdb parse_pdb(const path& name); // can throw parse_error

#endif
