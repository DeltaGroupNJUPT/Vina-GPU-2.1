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

#ifndef VINA_TERMS_H
#define VINA_TERMS_H

#include <boost/ptr_container/ptr_vector.hpp> 
#include "model.h"

/*
	定义一个字符串name，
	virtual ~term()虚函数；
	*/

struct term {
	std::string name;
	virtual ~term() {}
};


/*	distance_additive继承结构体term，可以访问赋值term中的public成员；
	定义一个double cutoff，
	声明了一个double虚函数eval,输入为两个atom_base结构体a、b和一个double r;
	virtual ~distance_additive()虚函数
*/

struct distance_additive : public term {                                        //继承term结构体的public
	fl cutoff;
	distance_additive(fl cutoff_) : cutoff(cutoff_) {}
	virtual fl eval(const atom_base& a, const atom_base& b, fl r) const = 0;
	virtual ~distance_additive() {}
};


/*	usable继承结构体distance_additive，可以访问和赋值distance_additive的public成员；
	这里该结构体调用初始化，将double cuoff赋值给distance_additive结构体中的double cuoff，
	并对定义为atom_type类中枚举类型t的atom_typing_used赋值为t枚举元素中的XS=2；
	定义了eval函数，输入为一个atom_base结构体a、atom_base结构体以及double r，结构体atom_base继承了结构体atom_type，
	所以a.get（）调用了atom_type的成员，可以得知a.get(XS)=b.get(XS)=XS_TYPE_SIZE，double r=r；
	定义了一个double虚函数eval可以重写,输入为两个unint t1、t2和一个double r；
	virtual ~usable()虚函数

*/
struct usable : public distance_additive {
	atom_type::t atom_typing_used;             //凡是被说明为atom_type类中枚举类型t的类型变量的取值只能是枚举元素中的某一个值。
	usable(fl cutoff_) : distance_additive(cutoff_), atom_typing_used(atom_type::XS) {}     //XS=2
	fl eval(const atom_base& a, const atom_base& b, fl r) const { // should not be overriden
		return eval(a.get(atom_typing_used), b.get(atom_typing_used), r);       //get（）函数是结构体atom_base成员，
																				//因为atom_typing_use=XS，所以返回XS_TYPE_SIZE
	}
	virtual fl eval(sz t1, sz t2, fl r) const { return 0; } 
	virtual ~usable() {}
};



/* additive继承结构体term
定义一个double cutoff，然后被赋值为计算机能识别的double数可表示的最大值00DB17EE；
定义了一个double虚函数eval,输入为model结构体m和两个atom_base结构体i、j;
virtual ~additive() 虚函数
*/

struct additive : public term {
	fl cutoff;
	additive() : cutoff(max_fl) {}                                    //计算机能识别的double数可表示的最大值，我的计算机能识别的最大double数是00DB17EE
	virtual fl eval(const model& m, const atom_index& i, const atom_index& j) const = 0;
	virtual ~additive() {}
};



/* intermolecular继承结构体term

声明了一个double虚函数eval,输入为model结构体m；

*/


struct intermolecular : public term {
	virtual fl eval(const model& m) const = 0;
};




/*
		定义了num_tors、num_rotors、num_heavy_atoms、num_hydrophobic_atoms、
		ligand_max_num_h_bonds、num_ligands、ligand_lengths_sum的double类型数据；
		声明了函数operator flv()和函数 get_names()
		input——m：model结构体；
		声明了num_bonded_heavy_atoms（）和atom_rotors()函数，不过没有定义；
		然后就是侵入式序列化，可以将num_tors、num_rotors、num_heavy_atoms、num_hydrophobic_atoms、
		ligand_max_num_h_bonds、num_ligands、ligand_lengths_sum这些结构体成员存储在文件中或者通过网络发出去，
		调用时会具体给出传输类型。
		https://blog.csdn.net/zj510/article/details/8105408
*/
struct conf_independent_inputs {
	fl num_tors;
	fl num_rotors;
	fl num_heavy_atoms;
	fl num_hydrophobic_atoms;
	fl ligand_max_num_h_bonds;
	fl num_ligands;
	fl ligand_lengths_sum;
	operator flv() const;                        //声明函数operator flv()
	conf_independent_inputs(const model& m);
	std::vector<std::string> get_names() const;  //声明函数 get_names()
	conf_independent_inputs();
private:
	unsigned num_bonded_heavy_atoms(const model& m, const atom_index& i) const; // FIXME? - could be static, but I don't feel like declaring function friends
	unsigned atom_rotors(const model& m, const atom_index& i) const; // the number of rotatable bonds to heavy ligand atoms

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & num_tors;
		ar & num_rotors;
		ar & num_heavy_atoms;
		ar & num_hydrophobic_atoms;
		ar & ligand_max_num_h_bonds;
		ar & num_ligands;
		ar & ligand_lengths_sum;
	}
};


/*
定义两个虚函数eval()和size()（纯虚函数） 
变量in:conf_independent_inputs结构体；x：double；it：迭代器；
flv::const_iterator 迭代器定义，迭代器是一种检查vector内元素并遍历元素的数据类型，
const_iterator只能读取vector中的元素，不能用于改变其值。
*/
struct conf_independent : public term {
	virtual fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& it) const = 0;
	virtual sz size() const = 0; // how many parameters does it take
}; 



/*
inter——enabled：布尔型向量，fun：任意类型的指针向量，这个和vector类似不过效率更高；
			 term_set继承distance_additive，
			 term_set继承usable，
			 term_set继承additive，
			 term_set继承intermolecular，
			 term_set继承conf_independent。
*/

template<typename T>
struct term_set {
	std::vector<bool> enabled;
	boost::ptr_vector<T> fun; // FIXME? const T?    //任意类型的指针向量fun,这个和vector类似不过效率更高

	/*
		input——e:unsigned，f——任意类型指针，enabled：布尔型向量，fun：任意类型的指针向量，这个和vector类似不过效率更高；
		该函数作用是在enabled向量最末尾插入0（e小于等于0）、1（e大于0），
		并在fun最末尾添加输入f。
	*/
	void add(unsigned e, T* f) { // FIXME? const T* ?
		enabled.push_back(e > 0);      //向量最末尾插入bool类型数据
		fun.push_back(f);
	}

	/*
	inter——enabled：布尔型向量；
	output——tmp:unsigned int 初始值为0；
	该函数返回enabled中大于等于1的值的个数。
	*/
	sz num_enabled() const {
		sz tmp = 0;
		VINA_FOR_IN(i, enabled)  // for (uint VINA_MACROS_TMP = enabled.size, i = 0; i<VINA_MACROS_TMP; i++)
			if(enabled[i])      //计算enabled向量中为1的值
				++tmp;
		return tmp;
	}

/*
input——enabled_only：bool；out——字符串向量；
enabled：布尔型向量，fun：任意类型的指针向量，这个和vector类似不过效率更高；
先确定enabled向量的元素个数是否和fun向量元素相同，是的话继续否则终止；
最后如果enabled_only为0或者enabled[i]大于等于1则在out字符串向量末尾添加fun[i].name（term结构体中name）
*/
	void get_names(bool enabled_only, std::vector<std::string>& out) const { // appends to "out"
		VINA_CHECK(enabled.size() == fun.size());       //assert
		VINA_FOR_IN(i, fun)								// for (uint VINA_MACROS_TMP = fun.size, i = 0; i<VINA_MACROS_TMP; i++)
			if(!enabled_only || enabled[i])             
				out.push_back(fun[i].name);              //term结构体中name
	}



	/*
	input——in:const_iterator迭代器，out：double向量
	inter——enabled：布尔型向量，fun：任意类型的指针向量，这个和vector类似不过效率更高；
	先确定enabled向量的元素个数是否和fun向量元素相同，是的话继续否则终止；
	最后当enabled中的值大于等于1时将迭代器里第i元素添加到out向量最末尾；
	in++循环enabled.size次，不在if语句中
	*/
	void filter(flv::const_iterator& in, flv& out) const {
		VINA_CHECK(enabled.size() == fun.size());   // for (uint VINA_MACROS_TMP = enabled.size, i = 0; i<VINA_MACROS_TMP; i++)
		VINA_FOR_IN(i, enabled) {
			if(enabled[i])
				out.push_back(*in);                 //将迭代器里第i元素添加到out向量最末尾
			++in;
		}
	}


	/*
	inter——tmp：double
	enabled：布尔型向量，fun：任意类型的指针向量，这个和vector类似不过效率更高；
	通过一个循环寻找fun[i].cutoff中的最大值并返回（additive结构体中cutoff）
	*/
	fl max_cutoff() const {
		fl tmp = 0; 
		VINA_FOR_IN(i, fun)                          // for (uint VINA_MACROS_TMP = fun.size, i = 0; i<VINA_MACROS_TMP; i++)
			tmp = (std::max)(tmp, fun[i].cutoff);   //additive结构体中cutoff
		return tmp;
	}

	/*
	    fun：任意类型的指针向量，这个和vector类似不过效率更高
		该函数返回fun向量元素个数
	*/
	sz size() const { return fun.size(); }

	/*
		input——unint i；
		fun：任意类型的指针向量，这个和vector类似不过效率更高
		返回fun[i]的值；
	
	*/
	const T& operator[](sz i) const { return fun[i]; }
};




/*
	定义了double类型的向量e、i；
	定义size()函数返回e向量和i向量的元素总个数；
	定义num_weights（）函数，返回e和i向量中最大的元素个数；
	声明了eval（）函数，返回类型为double；
	然后就是侵入式序列化，可以将e、i这些结构体成员存储在文件中或者通过网络发出去，
	调用时会具体给出传输类型。
	https://blog.csdn.net/zj510/article/details/8105408

*/
struct factors {
	flv e; // external
	flv i; // internal
	sz size() const { return e.size() + i.size(); }       //返回e向量和i向量的元素总个数
	//sz num_weights() const { return (std::max)(e.size(), i.size()); } // FIXME? compiler bug? getting warnings here
	sz num_weights() const { return (e.size() > i.size()) ? e.size() : i.size(); }
	fl eval(const flv& weights, bool include_internal) const;      //声明
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & e;
		ar & i;
	}
};


/*
			 term_set继承distance_additive，然后构造一个distance_additive_terms结构体
			 term_set继承usable，然后构造一个usable_terms结构体
			 term_set继承additive，然后构造一个additive_terms结构体
			 term_set继承intermolecular，然后构造一个intermolecular_terms结构体
			 term_set继承conf_independent，然后构造一个conf_independent_terms结构体

*/

struct terms {
	term_set<distance_additive> distance_additive_terms;    //term_set继承distance_additive，然后构造一个distance_additive_terms结构体
	term_set<usable>            usable_terms;				//term_set继承usable，然后构造一个usable_terms结构体
	term_set<additive>          additive_terms;				//term_set继承additive，然后构造一个additive_terms结构体
	term_set<intermolecular>    intermolecular_terms;       //term_set继承intermolecular，然后构造一个intermolecular_terms结构体
	term_set<conf_independent>  conf_independent_terms;     //term_set继承conf_independent，然后构造一个conf_independent_terms结构体    

	// the class takes ownership of the pointer with 'add'
	void add(unsigned e, distance_additive* p) { distance_additive_terms.add(e, p); }//	在fun最末尾添加p，在enable最末尾添加0，1
	void add(unsigned e, usable* p)            {            usable_terms.add(e, p); }//	在fun最末尾添加p，在enable最末尾添加0，1
	void add(unsigned e, additive* p)          {          additive_terms.add(e, p); }//	在fun最末尾添加p，在enable最末尾添加0，1
	void add(unsigned e, intermolecular* p)    {    intermolecular_terms.add(e, p); }//	在fun最末尾添加p，在enable最末尾添加0，1
	void add(unsigned e, conf_independent* p)  {  conf_independent_terms.add(e, p); }//	在fun最末尾添加p，在enable最末尾添加0，1

	std::vector<std::string> get_names(bool enabled_only) const; // does not include conf-independent 声明函数
	sz size_internal() const;      //声明，该函数返回distance_additive_terms结构体中fun向量元素个数、
									//usable_terms结构体中fun向量元素个数、
	                                //additive_terms中fun向量元素个数的总和
	sz size() const { return size_internal() + intermolecular_terms.size(); }//// does not include conf-independent 声明函数
								   //该函数返回distance_additive_terms结构体中fun向量元素个数、
								   //usable_terms结构体中fun向量元素个数、
								   //additive_terms中fun向量元素个数、
									//intermolecular_terms的向量元素个数总和
	sz size_conf_independent(bool enabled_only) const; // number of parameters does not necessarily equal the number of operators
	fl max_r_cutoff() const;             //声明
	flv evale(const model& m) const;    //声明
	flv evali(const model& m) const;    //声明
	flv evale_robust(const model& m) const;//声明
	factors eval(const model& m) const;    //返回结构体的函数shming
	fl eval_conf_independent(const conf_independent_inputs& in, fl x, flv::const_iterator& it) const;  //声明
	flv filter_external(const flv& v) const;  //声明
	flv filter_internal(const flv& v) const;  //声明
	factors filter(const factors& f) const;   //声明
	void display_info() const;      //用于打印一些重要的信息
private:
	void eval_additive_aux(const model& m, const atom_index& i, const atom_index& j, fl r, flv& out) const; // out is added to
	//声明
};

#endif
