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

#ifndef VINA_MODEL_H
#define VINA_MODEL_H

#include <boost/optional.hpp> // for context

#include "file.h"
#include "tree.h"
#include "matrix.h"
#include "precalculate.h"
#include "igrid.h"
#include "grid_dim.h"


/*
*声明一个结构体类型interacting_pair，数据成员包括unsigned int类型的type_pair_index，a，b
*/

struct interacting_pair {
	sz type_pair_index;
	sz a;
	sz b;
	interacting_pair(sz type_pair_index_, sz a_, sz b_) : type_pair_index(type_pair_index_), a(a_), b(b_) {}
	//构造函数，当创建结构体对应对象时，参数列表传入三个unsigned int类型数时，按顺序分别初始化type_pair_index，a，b
};

typedef std::vector<interacting_pair> interacting_pairs;
//给interacting_pair类型容器取别名interacting_pairs


//定义一个结构体命名为parsed_line，其可于一个单元内存储两个相异对象，分别是string类
//和unsigned int类型的optional类（它很像一个仅能存放一个元素的容器，它实现了"未初始化"的概念：
//如果元素未初始化，那么容器就是空的，否则，容器内就是有效的，已经初始化的值。）

typedef std::pair<std::string, boost::optional<sz> > parsed_line;


typedef std::vector<parsed_line> context;
//给parsed_line类型容器取别名context

/*ligand类,继承flexible_body, atom_range
*成员包括unsigned类型degrees_of_freedom，interacting_pairs容器pairs，context容器cont
*构造函数residue传入一个flexible_body对象引用f初始化flexible_body，和一个unsigned数degrees_of_freedom，给基类atom_range传入参数0，0
*声明一个函数set_range()
*/
struct ligand : public flexible_body, atom_range {
	unsigned degrees_of_freedom; // can be different from the apparent number of rotatable bonds, because of the disabled torsions
	interacting_pairs pairs;
	context cont;
	ligand(const flexible_body& f, unsigned degrees_of_freedom_) : flexible_body(f), atom_range(0, 0), degrees_of_freedom(degrees_of_freedom_) {}
	void set_range();
};



/*residue类,继承main_branch
*构造函数residue传入一个main_branch对象引用m初始化main_branch
*/
struct residue : public main_branch {
	residue(const main_branch& m) : main_branch(m) {}
};

enum distance_type {DISTANCE_FIXED, DISTANCE_ROTOR, DISTANCE_VARIABLE};
//定义枚举类型distance_type，枚举成员包括DISTANCE_FIXED, DISTANCE_ROTOR, DISTANCE_VARIABLE三个字面值常量

//用distance_type实例化strictly_triangular_matrix类模板，并取别名为distance_type_matrix
typedef strictly_triangular_matrix<distance_type> distance_type_matrix;

struct non_cache;                 // 向前声明结构体non_cache
struct naive_non_cache;           // 向前声明结构体naive_non_cache
struct cache;                     // 向前声明结构体cache
struct szv_grid;                  // 向前声明结构体szv_grid
struct terms;                     // 向前声明结构体terms
struct conf_independent_inputs;   // 向前声明结构体conf_independent_inputs
struct pdbqt_initializer;         // 向前声明结构体pdbqt_initializer - only declared in parse_pdbqt.cpp
struct model_test;                // 向前声明结构体model_test



//定义结构体model

struct model {
    //公有成员
	//声明append函数，传入model对象的引用
	void append(const model& m);

	//定义函数atom_typing_used，返回一个t类型的数据m_atom_typing_used
	atom_type::t atom_typing_used() const { return m_atom_typing_used; }

	sz num_movable_atoms() const { return m_num_movable_atoms; } 
	//定义一个返回unsigned int类型数据的函数num_movable_atoms，返回私有成员m_num_movable_atoms
	sz num_internal_pairs() const; 
	//定义一个函数num_other_pairs，返回other_pairs容器的大小（unsigned int类型数据）
	sz num_other_pairs() const { return other_pairs.size(); }
	//定义一个函数num_other_pairs，返回ligands容器的大小（unsigned int类型数据）
	sz num_ligands() const { return ligands.size(); }
	//定义一个函数num_flex，返回flex容器的大小（unsigned int类型数据）
	sz num_flex() const { return flex.size(); }

	//定义一个函数ligand_degrees_of_freedom，传入一个unsigned int类型数据ligand_number，返回ligands容器第ligand_number个元素的成员degrees_of_freedom（unsigned int类型）
	sz ligand_degrees_of_freedom(sz ligand_number) const { return ligands[ligand_number].degrees_of_freedom; }
	//声明一个函数ligand_longest_branch，传入一个unsigned int类型数据ligand_number，输出一个unsigned int类型返回值
	sz ligand_longest_branch(sz ligand_number) const;
	//声明一个函数ligand_length，传入一个unsigned int类型数据ligand_number，输出一个unsigned int类型返回值
	sz ligand_length(sz ligand_number) const;

	//声明一个常量函数get_movable_atom_types，传入一个t类型数据atom_typing_used_，返回一个unsigned int容器
	szv get_movable_atom_types(atom_type::t atom_typing_used_) const;

	//声明一个常量函数get_size，返回一个conf_size类
	conf_size get_size() const;
	//声明一个常量函数get_initial_conf，返回一个conf类
	conf get_initial_conf() const; // torsions = 0, orientations = identity, ligand positions = current

	//声明一个常量函数movable_atoms_box，输入double类型add_to_each_dimension，granularity，其中granularity默认为0.375，返回一个grid_dims类
	grid_dims movable_atoms_box(fl add_to_each_dimension, fl granularity = 0.375) const;


	/*
	* 函数：write_flex，常量函数
	* input：path类对象常量引用name，string类型常量引用remark
	* output: 无
	* function：调用write_context函数，传入model成员flex_context，name和remark
	*/
	void write_flex  (                  const path& name, const std::string& remark) const { write_context(flex_context, name, remark); }
	
	/*
	* 函数：write_ligand，常量函数
	* input：unsigned int类型ligand_number，path类对象常量引用name，string类型常量引用remark
	* output: 无
	* function：检查ligand_number是否小于model成员ligands的大小，是则调用write_context函数，传入model成员ligands的ligand_number个元素成员cont，name和remark，否则报错，终止程序的执行
	*/
	void write_ligand(sz ligand_number, const path& name, const std::string& remark) const { VINA_CHECK(ligand_number < ligands.size()); write_context(ligands[ligand_number].cont, name, remark); }
	
	
	
	/*
	* 函数：write_structure，常量函数
	* input：文件输出流ofile类out
	* output: 无
	* function：1.遍历ligands容器元素，分别调用write_context函数，传入ligands容器第i个元素的成员cont容器和out
	* 2.判断flex容器大小是否大于0，是则调用write_context函数，传入flex_context容器和out
	*/
	void write_structure(ofile& out) const {
		VINA_FOR_IN(i, ligands)
			write_context(ligands[i].cont, out);
		if(num_flex() > 0) // otherwise remark is written in vain
			write_context(flex_context, out);
	}

	/*
	* 函数：重载write_structure，常量函数
	* input：文件输出流ofile类out，string类remark
	* output: 无
	* function：1.向文件输出流输出remark
	* 2.调用原write_structure函数，传入文件输出流out
	*/
	void write_structure(ofile& out, const std::string& remark) const {
		out << remark;
		write_structure(out);
	}

	/*
	* 函数：重载write_structure，常量函数
	* input：path类对象引用name
	* output: 无
	* function：1.通过path类创建文件输出流out
	* 2.调用原write_structure函数，传入文件输出流out
	*/
	void write_structure(const path& name) const { ofile out(name); write_structure(out); }

	/*
	* 函数：write_model，常量函数
	* input：文件输出流ofile类out，unsigned int数model_number，string类remark
	* output: 无
	* function：向文件输出流输出字符串：MODEL model_number（换行）
	* 2.调用write_structure重载版本函数，传入文件输出流out和remark
	* 3.向文件输出流输出字符串：ENDMDL（换行）
	*/
	void write_model(ofile& out, sz model_number, const std::string& remark) const {
		out << "MODEL " << model_number << '\n';
		write_structure(out, remark);
		out << "ENDMDL\n";
	}

	//声明函数seti，输入一个conf类引用，无返回值
	void seti(const conf& c);
	//声明函数sete，输入一个conf类引用，无返回值
	void sete(const conf& c);
	//声明函数set，输入一个conf类引用，无返回值
	void set (const conf& c);

	//声明常量函数 gyration_radius，输入一个unsigned int数ligand_number，返回double数
	fl gyration_radius(sz ligand_number) const; // uses coords

	/*
	* 函数：movable_atom，常量函数
	* input：unsigned int数i
	* output: atom_base类引用
	* function：判断i的是否小于m_num_movable_atoms，是则返回atoms容器第i个元素，否则报错，终止程序的执行
	*/
	const atom_base& movable_atom  (sz i) const { assert(i < m_num_movable_atoms); return  atoms[i]; }

	/*
	* 函数：movable_coords，常量函数
	* input：unsigned int数i
	* output: 长度为3的矢量vec
	* function：判断i的是否小于m_num_movable_atoms，是则返回coords容器第i个元素，否则报错，终止程序的执行
	*/
	const vec&       movable_coords(sz i) const { assert(i < m_num_movable_atoms); return coords[i]; }


	//声明常量函数atom_coords，输入atom_index类常量引用，返回vec矢量常量引用
	const vec& atom_coords(const atom_index& i) const;

	//声明常量函数distance_sqr_between，输入两个atom_index类常量引用，返回double类型数
	fl distance_sqr_between(const atom_index& a, const atom_index& b) const;

	//声明常量函数atom_exists_between，输入distance_type_matrix类常量引用，两个atom_index类常量引用，unsigned int类型容器常量引用，返回一个布尔值
	bool atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const; // there is an atom closer to both a and b then they are to each other and immobile relative to them

	//声明常量函数distance_type_between，输入distance_type_matrix类常量引用，两个atom_index类常量引用，返回一个枚举类型distance_type
	distance_type distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const;

	// clean up
	//声明常量函数evali，输入precalculate类常量引用，vec矢量常量引用，返回double类型数
	fl evali     (const precalculate& p,                  const vec& v                          ) const;
	//声明常量函数evale，输入precalculate类常量引用，igrid类常量引用，vec矢量常量引用，返回double类型数
	fl evale     (const precalculate& p, const igrid& ig, const vec& v                          ) const;
	//声明函数eval，输入precalculate类常量引用，igrid类常量引用，vec矢量常量引用，conf类常量引用，返回double类型数
	fl eval      (const precalculate& p, const igrid& ig, const vec& v, const conf& c           );
	//声明函数eval_deriv，输入precalculate类常量引用，igrid类常量引用，vec矢量常量引用，conf类常量引用，change类引用，返回double类型数
	fl eval_deriv(const precalculate& p, const igrid& ig, const vec& v, const conf& c, change& g);


	//声明函数eval_intramolecular，输入precalculate类常量引用，vec矢量常量引用，conf类常量引用，返回double类型数
	fl eval_intramolecular(                            const precalculate& p,                  const vec& v, const conf& c);
	//声明函数eval_adjusted，输入scoring_function类常量引用，输入precalculate类常量引用，igrid类常量引用，vec矢量常量引用，conf类常量引用，double类型数，返回double类型数
	fl eval_adjusted      (const scoring_function& sf, const precalculate& p, const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy);

	//声明常量函数rmsd_lower_bound，输入model类对象，返回double类型数
	fl rmsd_lower_bound(const model& m) const; // uses coords
	//声明常量函数rmsd_upper_bound，输入model类对象，返回double类型数
	fl rmsd_upper_bound(const model& m) const; // uses coords
	//声明常量函数rmsd_ligands_upper_bound，输入model类对象，返回double类型数
	fl rmsd_ligands_upper_bound(const model& m) const; // uses coords

	//声明常量函数verify_bond_lengths，无返回
	void verify_bond_lengths() const;
	//声明常量函数about，无返回
	void about() const;


	/*
	* 函数：get_ligand_internal_coords，常量函数
	* input：无
	* output: vec类型容器
	* function：1.判断model成员ligands容器的大小是否等于1，是则往下执行，否则报错，终止程序的执行
	* 2.建立一个临时vecv对象tmp
	* 3.返回model成员ligands容器中起始元素的引用,命名lig
	* 4.i从lig成员begin的大小遍历到lig成员begin的大小，依次往tmp尾部添加model成员internal_coords容器第i个元素
	* 5.返回tmp
	*/
	vecv get_ligand_internal_coords() const { // FIXME rm
		VINA_CHECK(ligands.size() == 1);
		vecv tmp;
		const ligand& lig = ligands.front();
		VINA_RANGE(i, lig.begin, lig.end)
			tmp.push_back(internal_coords[i]);
		return tmp;
	}


	/*
	* 函数：get_ligand_coords，常量函数
	* input：无
	* output: vec类型容器
	* function：1.判断model成员ligands容器的大小是否等于1，是则往下执行，否则报错，终止程序的执行
	* 2.建立一个临时vecv对象tmp
	* 3.返回model成员ligands容器中起始元素的引用,命名lig
	* 4.i从lig成员begin的大小遍历到lig成员begin的大小，依次往tmp尾部添加model成员coords容器第i个元素
	* 5.返回tmp
	*/
	vecv get_ligand_coords() const { // FIXME rm
		VINA_CHECK(ligands.size() == 1);
		vecv tmp;
		const ligand& lig = ligands.front();
		VINA_RANGE(i, lig.begin, lig.end)
			tmp.push_back(coords[i]);
		return tmp;
	}

	/*
	* 函数：get_heavy_atom_movable_coords，常量函数
	* input：无
	* output: vec类型容器
	* function：1.建立一个临时vecv对象tmp
	* 2.i从0遍历到m_num_movable_atoms，依次判断atoms第i个元素的成员el是否等于0，是则往tmp尾部添加model成员coords容器第i个元素
	* 3.返回tmp
	*/
	vecv get_heavy_atom_movable_coords() const { // FIXME mv
		vecv tmp;
		VINA_FOR(i, num_movable_atoms())
			if(atoms[i].el != EL_TYPE_H)
				tmp.push_back(coords[i]);
		return tmp;
	}

	//声明一个常量函数check_internal_pairs，无返回
	void check_internal_pairs() const;
	//声明一个常量函数print_stuff，无返回
	void print_stuff() const; // FIXME rm

	//声明一个常量函数clash_penalty，返回double类型数据
	fl clash_penalty() const;

	//私有成员
private:
	friend struct non_cache;              // 声明友元结构体non_cache，                        此结构体可访问结构体model()的私有成员                                           
	friend struct naive_non_cache;		  // 声明友元结构体naive_non_cache，                  此结构体可访问结构体model()的私有成员 			
	friend struct cache;				  // 声明友元结构体cache，                            此结构体可访问结构体model()的私有成员 					
	friend struct szv_grid;				  // 声明友元结构体szv_grid	，                        此结构体可访问结构体model()的私有成员 			
	friend struct terms;				  // 声明友元结构体terms，                            此结构体可访问结构体model()的私有成员 					
	friend struct conf_independent_inputs;// 声明友元结构体conf_independent_inputs，          此结构体可访问结构体model()的私有成员 	
	friend struct appender_info;		  // 声明友元结构体appender_info，                    此结构体可访问结构体model()的私有成员 			
	friend struct pdbqt_initializer;	  // 声明友元结构体pdbqt_initializer，                此结构体可访问结构体model()的私有成员  - only declared in parse_pdbqt.cpp
	friend struct model_test;			  // 声明友元结构体model_test，                       此结构体可访问结构体model()的私有成员 

	//构造函数model()，初始化m_num_movable_atoms为0，m_atom_typing_used为XS
	model() : m_num_movable_atoms(0), m_atom_typing_used(atom_type::XS) {};


	/*
	* 函数：get_atom，是常量函数
	* input：atom_index类引用i
	* function：判断i的成员in_grid是否为真，是则返回grid_atoms容器的第（对象i的成员i的数值）元素常量引用，否则返回atoms容器的第（对象i的成员i的数值）元素常量引用
	*/
	const atom& get_atom(const atom_index& i) const { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); }
	/*
	* 函数：get_atom
	* input：atom_index类引用i
	* function：判断i的成员in_grid是否为真，是则返回grid_atoms容器的第（对象i的成员i的数值）元素引用，否则返回atoms容器的第（对象i的成员i的数值）元素引用
	*/
	      atom& get_atom(const atom_index& i)       { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); }


	//1.声明一个常量函数write_context，输入context容器对象和一个文件输出流，无返回
	void write_context(const context& c, ofile& out) const;
	


	/*2.重载常量函数write_context，输入context容器对象、一个文件输出流和string类remark，无返回
	* function:向文件输出流输出remark
	*/
	void write_context(const context& c, ofile& out, const std::string& remark) const {
		out << remark;
	}

	/*3.重载常量函数write_context，输入context容器对象和一个文件path类对象，无返回
	* function:通过path类创建一个文件输出流out,再调用write_context函数第1个版本
	*/
	void write_context(const context& c, const path& name) const {
		ofile out(name);
		write_context(c, out);
	}

	/*4.重载常量函数write_context，输入context容器对象、一个文件path类对象和string类remark，无返回
	* function:通过path类创建一个文件输出流out,再调用write_context函数第2个版本
	*/
	void write_context(const context& c, const path& name, const std::string& remark) const {
		ofile out(name);
		write_context(c, out, remark);
	}

	//声明一个常量函数rmsd_lower_bound_asymmetric，输入model类x,y引用，返回一个double类型数
	fl rmsd_lower_bound_asymmetric(const model& x, const model& y) const; // actually static
	

	//声明一个常量函数sz_to_atom_index，输入unsigned int类型数据i，返回一个atom_index类
	atom_index sz_to_atom_index(sz i) const; // grid_atoms, atoms
	//声明一个常量函数bonded_to_HD，输入atom类引用a，返回布尔值
	bool bonded_to_HD(const atom& a) const;
	//声明一个常量函数bonded_to_heteroatom，输入atom类引用a，返回布尔值
	bool bonded_to_heteroatom(const atom& a) const;
	//声明一个常量函数find_ligand，输入unsigned int类型数据a，返回一个unsigned int类型数据
	sz find_ligand(sz a) const;
	//声明一个常量函数bonded_to，输入unsigned int类型数据a，n，和unsigned int类型容器引用szv，无返回
	void bonded_to(sz a, sz n, szv& out) const;
	//重载常量函数bonded_to，输入unsigned int类型数据a，n，返回unsigned int类型容器
	szv bonded_to(sz a, sz n) const;

	//声明一个函数assign_bonds，输入distance_type_matrix类引用，无返回
	void assign_bonds(const distance_type_matrix& mobility); // assign bonds based on relative mobility, distance and covalent length
	//声明一个函数assign_types，无返回
	void assign_types();
	//声明一个函数initialize_pairs，输入distance_type_matrix类引用，无返回
	void initialize_pairs(const distance_type_matrix& mobility);
	//声明一个函数initialize，输入distance_type_matrix类引用，无返回
	void initialize(const distance_type_matrix& mobility);
	//声明一个常量函数clash_penalty_aux，输入interacting_pair类容器引用，返回double数据
	fl clash_penalty_aux(const interacting_pairs& pairs) const;

	vecv internal_coords; //声明一个容器internal_coords，容器里面包含若干结构体vec（长度为3的矢量）
	




	


	vector_mutable<residue> flex; //声明一个vector_mutable容器flex，其元素为residue类
	context flex_context;//声明一个parsed_line类型容器数据

	//声明一个interacting_pairs容器other_pairs
	interacting_pairs other_pairs; // all except internal to one ligand: ligand-other ligands; ligand-flex/inflex; flex-flex/inflex

	sz m_num_movable_atoms; //声明一个unsigned int类型数据
	atom_type::t m_atom_typing_used;//声明一个t类型数据，t是枚举类型
public:
		vector_mutable<ligand> ligands; //声明一个vector_mutable容器ligands，其元素为ligand类
		atomv atoms; // movable, inflex
		vecv coords;          //声明一个容器coords，         容器里面包含若干结构体vec（长度为3的矢量）
		vecv minus_forces;    //声明一个容器minus_forces，   容器里面包含若干结构体vec（长度为3的矢量）	
		atomv grid_atoms;//声明容器grid_atoms，atoms，元素均为atom类
};


#endif

