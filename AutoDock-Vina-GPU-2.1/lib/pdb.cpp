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

#include "pdb.h"
#include "parse_error.h"
#include "file.h"
#include "convert_substring.h"


/*
* 结构体pdb成员函数check定义
* input:双精度浮点数min_distance，最小距离
* 该函数作用：结构体pdb的数据成员容器atoms中，对任意两个结构体的矢量距离小于输入双精度浮点数min_distance的平方的，都输出提醒
* output：无
*/
void pdb::check(fl min_distance) const {
	VINA_FOR_IN(i, atoms) { //循环atoms容器里所有pdb_atom结构体,分别对应容器下标为i
		const pdb_atom& a = atoms[i];//a引用atoms容器里第i个pdb_atom结构体
		VINA_RANGE(j, i+1, atoms.size()) {//循环atoms容器里所有从第i+1个pdb_atom结构体剩下的pdb_atom结构体（包括i+1)，分别对应容器下标为j
			const pdb_atom& b = atoms[j];//b引用atoms容器里第j个pdb_atom结构体
			fl d2 = vec_distance_sqr(a.coords, b.coords);//计算结构体a,b中对应coords矢量的平方距离，返回给double类型数据d2
			
			//如果d2小于输入的双精度浮点数min_distance的平方
			//输出提醒：
			//The distance between 结构体a的id:结构体a的name:结构体a的element and 结构体b的id:结构体b的name:结构体b的element is d2的开方（换行）
			if(d2 < sqr(min_distance)) {
				std::cout << "The distance between " 
					<< a.id << ":" << a.name << ":" << a.element
					<< " and " 
					<< b.id << ":" << b.name << ":" << b.element
					<< " is " << std::sqrt(d2) << '\n';
			}
		}
	}
}


/*
* 函数：string_to_pdb_atom
* input:string类型str
* output：结构体pdb_atom
*/

pdb_atom string_to_pdb_atom(const std::string& str) {
	//如果string类型str，字符数小于66，抛出异常类bad_conversion()，终止当前程序的执行
	if(str.size() < 66) throw bad_conversion(); // b-factor is in 61-66
	pdb_atom tmp;//定义一个结构体pdb_atom的对象tmp
	
	tmp.id           = convert_substring<unsigned>   (str,  7, 11);     //将string类型str从( 7-1)到(11-1)下标的字符（忽略前导空格）转换成unsigned       类型数据赋予对象tmp的成员id          
	tmp.name         = convert_substring<std::string>(str, 13, 16);     //将string类型str从(13-1)到(16-1)下标的字符（忽略前导空格）转换成std::string    类型数据赋予对象tmp的成员name        
	tmp.residue_id   = convert_substring<int>        (str, 23, 26);     //将string类型str从(23-1)到(26-1)下标的字符（忽略前导空格）转换成int            类型数据赋予对象tmp的成员residue_id  
	tmp.residue_name = convert_substring<std::string>(str, 18, 20);     //将string类型str从(18-1)到(20-1)下标的字符（忽略前导空格）转换成string         类型数据赋予对象tmp的成员residue_name
	tmp.coords[0]    = convert_substring<fl>         (str, 31, 38);     //将string类型str从(31-1)到(38-1)下标的字符（忽略前导空格）转换成double         类型数据赋予对象tmp的成员coords[0]   
	tmp.coords[1]    = convert_substring<fl>         (str, 39, 46);     //将string类型str从(39-1)到(46-1)下标的字符（忽略前导空格）转换成double         类型数据赋予对象tmp的成员coords[1]   
	tmp.coords[2]    = convert_substring<fl>         (str, 47, 54);     //将string类型str从(47-1)到(54-1)下标的字符（忽略前导空格）转换成double         类型数据赋予对象tmp的成员coords[2]   
	tmp.b_factor     = convert_substring<fl>         (str, 61, 66);     //将string类型str从(61-1)到(66-1)下标的字符（忽略前导空格）转换成double         类型数据赋予对象tmp的成员b_factor    
	tmp.element      = convert_substring<std::string>(str, 77, 78);     //将string类型str从(77-1)到(78-1)下标的字符（忽略前导空格）转换成string         类型数据赋予对象tmp的成员element     
	//返回结构体对象tmp
	return tmp;
}

/*
* 函数：parse_pdb
* input:path类，定义对象时，参数传递文件路径
* output：读取文件数据并保存其中的结构体pdb
*/
pdb parse_pdb(const path& name) {
	//定义一个结构体ifile的对象in
	ifile in(name);

	//定义一个结构体pdb为tmp
	pdb tmp;

	//定义string类型str
	std::string str;
	
	//定义无符号整型count，初始值为0
	unsigned count = 0;

	//从输入流in读入字符，存到string变量str中，每次读入一个新行，直到读入了文件结束标志
	while(std::getline(in, str)) {
	//count用于保存行号
		++count;
	//判断当前字符串是否以另外一个给定的子字符串“ATOM  ”或者"HETATM"开头，是则执行if语句里的操作，否则继续读取下一行
		if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				tmp.atoms.push_back(string_to_pdb_atom(str));
				//将pdb_atom结构体对象string_to_pdb_atom存入tmp的数据成员，容器atoms
			}
			catch(...) { // bad_conversion, presumably; but can lexical_cast throw its own errors?
				throw parse_error(name, count, "ATOM syntax incorrect");
				//在将str转换成string_to_pdb_atom数据成员的过程发成异常，抛出异常类parse_error()，终止当前程序的执行
			}
		}
	}
	return tmp;
	//返回读入文件数据后的pdb对象tmp
}
