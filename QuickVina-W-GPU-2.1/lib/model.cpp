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

#include "model.h"
#include "file.h"
#include "curl.h"


/*
* 函数模板get_atom_range,获取原子距离
* input:T类型对象t（如：ligand类）
* output：atom_range类
* 作用：在一个节点上，拷贝一个t的成员node(atom_range类),并将t的成员children容器的所有元素的成员begin的最小值赋值给该拷贝的成员begin，容器所有元素的成员end的最大值赋值给该拷贝的成员end
*        完整说明见图一
*/

template<typename T>
atom_range get_atom_range(const T& t) {
	atom_range tmp = t.node;//定义一个临时atom_range对象tmp，存放t的成员node
	//遍历t的成员children容器,获得容器所有元素的成员begin的最小值并赋值给tmp的成员begin，获得容器所有元素的成员end的最大值并赋值给tmp的成员end
	VINA_FOR_IN(i, t.children) {
		atom_range r = get_atom_range(t.children[i]);
		if(tmp.begin > r.begin) tmp.begin = r.begin;
		if(tmp.end   < r.end  ) tmp.end   = r.end;
	}
	return tmp;
}


//结构体branch_metrics，默认初始化unsigned int类型的数据成员length，corner2corner均为0
struct branch_metrics {
	sz length;
	sz corner2corner;
	branch_metrics() : length(0), corner2corner(0) {}
};


/*
* 函数模板get_branch_metrics
* input:T类型对象t
* output：branch_metrics类
* 作用：见图一解释
*/
template<typename T>
branch_metrics get_branch_metrics(const T& t) {
	//定义一个临时branch_metrics类对象tmp，保存返回值
	branch_metrics tmp;
	//判断t的成员children是否为空，为空直接返回默认初始化的branch_metrics类对象，否则继续执行以下代码，
	if(!t.children.empty()) {
		//定义unsigned int类型数据corner2corner_max，初始化为0
		sz corner2corner_max = 0;
		//定义unsigned int类型数据容器lengths
		szv lengths;

		//遍历t的成员children容器
		VINA_FOR_IN(i, t.children) {
			branch_metrics res = get_branch_metrics(t.children[i]);//递归调用
			if(corner2corner_max < res.corner2corner)//判断corner2corner_max是否小于此子分支的corner2corner，是则赋值corner2corner给corner2corner_max
				corner2corner_max = res.corner2corner;
			lengths.push_back(res.length + 1); //将子分支的长度加一存入容器lengths中   FIXME? weird compiler warning (sz -> unsigned)
		}

		//按升序排序lengths容器的元素
		std::sort(lengths.begin(), lengths.end());

		//返回lengths容器中末尾元素的引用，赋值给tmp的length成员
		tmp.length = lengths.back();

		//tmp的length成员拷贝给tmp的corner2corner成员
		tmp.corner2corner = tmp.length;
		//如果容器lengths的大小大于等于2，tmp的corner2corner成员赋值为lengths最大元素的两倍
		if(lengths.size() >= 2)
			tmp.corner2corner += lengths[lengths.size() - 1];//按语义，我觉得此处应该lengths.size() - 2，不知道是不是程序写错了
		//如果tmp的corner2corner成员小于corner2corner_max，将corner2corner_max的值赋值给corner2corner
		if(tmp.corner2corner < corner2corner_max)
			tmp.corner2corner = corner2corner_max;
	}
	return tmp;
}

/*
* model的成员函数ligand_longest_branch
* input:unsigned int类型ligand_number
* output: model的成员ligands容器的ligand_number个元素对应的branch_metrics类的length数据
*/
sz model::ligand_longest_branch(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).length;
}


/*
* model的成员函数length
* input:unsigned int类型ligand_number
* output: model的成员ligands容器的ligand_number个元素对应的branch_metrics类的corner2corner数据
*/
sz model::ligand_length(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).corner2corner;
}

/*
* ligand的成员函数set_range
* input:无
* output: 无
* 作用：获取ligand的成员node的begin，和其所有后代的成员node的begin，找到最小值赋值给begin
*       获取ligand的成员node的end，  和其所有后代的成员node的end，  找到最大值赋值给end
*/
void ligand::set_range() {
	atom_range tmp = get_atom_range(*this);
	begin = tmp.begin;
	end   = tmp.end;
}

/////////////////// begin MODEL::APPEND /////////////////////////

// FIXME hairy code - needs to be extensively commented, asserted, reviewed and tested

/*
* 结构体：appender_info
* 成员：1.unsigend int类型grid_atoms_size，m_num_movable_atoms，atoms_size
*       2.构造函数函数：
*           初始化时输入一个model对象，
*           model对象的成员grid_atoms容器大小赋值给grid_atoms_size，
*           model对象的成员m_num_movable_atoms赋值给m_num_movable_atoms，
*           model对象的成员atoms容器大小赋值给atoms_size
*/
struct appender_info {
	sz grid_atoms_size;
	sz m_num_movable_atoms;
	sz atoms_size;

	appender_info(const model& m) : grid_atoms_size(m.grid_atoms.size()), m_num_movable_atoms(m.m_num_movable_atoms), atoms_size(m.atoms.size()) {}
};

/*
* 类：appender
* 成员：private:1.appender_info类a_info
*               2.appender_info类b_info
*               3.常量函数new_grid_index：
*							input:unsigned int类型x
*							output:unsigned int类型
*							作用：根据appender类的布尔值is_a，为真返回x，为假返回a_info的成员grid_atoms_size+x
*       
*		public: 1.bool类型is_a
*			    2.构造函数：输入两个model类对象，分别绑定到a_info，b_info，is_a初始化为真
*               3.定义调用操作符（）
*							input:unsigned int类型x
*							output:unsigned int类型
*							作用：根据appender类的布尔值is_a，为真判断x是否小于a_info对象的m_num_movable_atoms，为真返回x
*																												为假返回b_info对象的m_num_movable_atoms+x
*															  为真判断x是否小于b_info对象的m_num_movable_atoms，为真返回a_info对象的m_num_movable_atoms+x															
*                											  												    为假返回a_info对象的atoms_size+x
*               4.重载调用操作符（）
* 							input:atom_index类对象x
* 							output:atom_index类
* 							作用：将x拷贝给tmp，判断tmp的成员in_grid是否为真，是则令tmp的成员i等于将输入tmp的成员inew_grid_index后的返回值
*																			  否则令tmp的成员i等于将输入tmp的成员operator后的返回值
*                                 返回tmp
*				5.常量函数update
*							input:interacting_pair类引用ip
* 							output:无
* 							作用：将ip的成员a输入operator()函数返回新的ip的成员a
*                 				  将ip的成员b输入operator()函数返回新的ip的成员b
*				6.重载update
*							input:矢量vec的引用
* 							output:无
* 							作用：无
*				7.重载update
*							input:ligand类引用lig
*							output:无
*							作用：（1）	1.将此appender类对象(*this)传入transform函数中
*									    2.获取lig的成员begin和end，计算diff = end + begin
*									    3.根据上述定义的定义调用操作符（），对*this使用调用操作符计算（*this)(begin)，并赋值给begin
*									    4.计算end = begin + diff  
* 								  （2） 1.更新lig的成员node的成员begin和end，及其所有后代的成员node的成员begin和end
* 									    2.后代获取方式为例如：后一代lig.children[i]，后二代lig.children[i][j],后三代lig.children[i][j][k]........等等，直到所有的后代没有孩子为止
* 									    3.所有个体更新方式均相同，举一例说明，其他均相同
* 									    4.计算diff= end - begin;将begin输入本appender类调用操作符（），获取新begin;重新计算end = begin + diff
* 								  （3）遍历lig成员pairs，将pairs的元素输入update(interacting_pair& ip)，获取更新后的pairs
* 								  （4）遍历lig成员cont，将cont的元素输入update(parsed_line& p)，获取更新后的cont
*				8.重载update
*							input:residue类引用r
* 							output:无
* 							作用：1.更新r的成员node的成员begin和end，及其所有后代的成员node的成员begin和end
* 							      2.后代获取方式为例如：后一代r.children[i]，后二代r.children[i][j],后三代r.children[i][j][k]........等等，直到所有的后代没有孩子为止
* 								  3.所有个体更新方式均相同，举一例说明，其他均相同
* 							      4.计算diff= end - begin;将begin输入本appender类调用操作符（），获取新begin;重新计算end = begin + diff
*               9.重载update
*							input:parsed_line类引用p
*                           output:无
*							作用：判断p的optional类对象是否为空，不为空则将optional类对象的元素输入函数operator()中，返回值重新赋值给optional类对象
*               10.重载update
*							input:atom类引用a
* 							output:无
* 							作用:遍历a的成员bonds的元素，将元素的成员connected_atom_index输入重载的operator函数，更新connected_atom_index
*				11.函数模板：append
*                  input:T类型的容器引用a，T类型的常量容器引用b         
*				   output:无				
* 				   作用：1.获取a容器的大小a_sz				
* 					     2.将容器b的元素添加到容器a尾部				
* 				         3.将is_a的值赋为真，更新容器a原来的元素，将其每个元素根据T依次传进对应重载版本update函数
* 				         4.将is_a的值赋为假，更新容器a新添加的元素，将其每个元素根据T依次传进对应重载版本update函数
*				12.函数模板：coords_append
* 				   input:T类型的容器引用a，T类型的常量容器引用b         
* 				   output:无				
* 				   作用：1.构造b的拷贝b_copy
*						 2.将is_a的值赋为真，更新容器a的元素，将其每个元素根据T依次传进对应重载版本update函数
*						 3.将is_a的值赋为假，更新容器b_copy的元素，将其每个元素根据T依次传进对应重载版本update函数
*					     4.将b_copy容器头部一部分元素（b_info成员m_num_movable_atoms大小数目）插入a的某一元素（下标为a_info成员m_num_movable_atomsb_info大小）前
* 					     5.将b_copy容器剩下元素插入a的尾部
*/

class appender {
	appender_info a_info;
	appender_info b_info;
	sz new_grid_index(sz x) const {
		return (is_a ? x : (a_info.grid_atoms_size + x)); // a-grid_atoms spliced before b-grid_atoms
	}
public:
	bool is_a;

	appender(const model& a, const model& b) : a_info(a), b_info(b), is_a(true) {}

	sz operator()(sz x) const { // transform coord index
		if(is_a) {
			if(x < a_info.m_num_movable_atoms)  return x; // a-movable unchanged
			else                                return x + b_info.m_num_movable_atoms; // b-movable spliced before a-inflex
		}
		else {
			if(x < b_info.m_num_movable_atoms)  return x + a_info.m_num_movable_atoms; // a-movable spliced before b-movable
			else                                return x + a_info.atoms_size; // all a's spliced before b-inflex
		}
	}
	atom_index operator()(const atom_index& x) const { // transform atom_index
		atom_index tmp(x);
		if(tmp.in_grid) tmp.i = new_grid_index(tmp.i);
		else            tmp.i = operator()(tmp.i);
		return tmp;
	}

	// type-directed old -> new transformations
	void update(interacting_pair& ip) const {
		ip.a = operator()(ip.a);
		ip.b = operator()(ip.b);
	}
	void update(vec& v) const { // coordinates & forces - do nothing
	}
	void update(ligand& lig) const {
		lig.transform(*this); // ligand as an atom_range subclass
		transform_ranges(lig, *this);
		VINA_FOR_IN(i, lig.pairs)
			this->update(lig.pairs[i]);
		VINA_FOR_IN(i, lig.cont)
			this->update(lig.cont[i]); // parsed_line update, below
	}
	void update(residue& r) const {
		transform_ranges(r, *this);
	}
	void update(parsed_line& p) const {
		if(p.second)
			p.second = operator()(p.second.get());
	}
	void update(atom& a) const {
		VINA_FOR_IN(i, a.bonds) {
			bond& b = a.bonds[i];
			b.connected_atom_index = operator()(b.connected_atom_index); // atom_index transformation, above
		}
	}

	// ligands, flex, flex_context, atoms; also used for other_pairs
	template<typename T>
	void append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbbbbbbb
		sz a_sz = a.size();
		vector_append(a, b);

		is_a = true;
		VINA_FOR(i, a_sz)
			update(a[i]);

		is_a = false;
		VINA_RANGE(i, a_sz, a.size())
			update(a[i]);
	}

	// internal_coords, coords, minus_forces, atoms
	template<typename T>
	void coords_append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbaab
		std::vector<T> b_copy(b); // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

		is_a = true;
		VINA_FOR_IN(i, a)
			update(a[i]);

		is_a = false;
		VINA_FOR_IN(i, b_copy)
			update(b_copy[i]);

		// interleave 
		typedef typename std::vector<T>::const_iterator const cci;
		cci b1 = b_copy.begin();
		cci b2 = b_copy.begin() + b_info.m_num_movable_atoms;
		cci b3 = b_copy.end();

		a.insert(a.begin() + a_info.m_num_movable_atoms , b1 , b2);
		a.insert(a.end()                                , b2 , b3);
	}
};


/*
* 定义model的成员函数append
* input:model类常量引用m
* ouput:无
* 作用：1.检查本model类对象的m_atom_typing_used与输入的model类对象的m_atom_typing_used是否相等，不相等则报错，终止程序执行
*       2.用本model类对象与输入model类对象创建一个appender类t
*		3.将输入model类对象成员other_pairs容器的元素添加到本model类对象成员other_pairs容器尾部，并根据model.append函数的更新算法对本model类对象成员other_pairs容器的元素进行更新	
*		4.i从0循环到（本model成员atoms的大小-1）
*				j从0循环到（输入model成员atoms的大小-1）
*						如果i大于等于本model成员m_num_movable_atoms并且j大于等于输入model成员m_num_movable_atoms，则直接执行下一次循环
*						将本model成员atoms的第i个元素绑定到a，输入model成员atoms的第j个元素绑定到b
*						将本model成员m_atom_typing_used输入a的get函数，获得返回值t1，即m_atom_typing_used对应小写形式对应的unsignd int数
*						将输入model成员m_atom_typing_used输入b的get函数，获得返回值t2，即m_atom_typing_used对应小写形式对应的unsignd int数
*						将本model成员m_atom_typing_used输入函数num_atom_types获得返回值n
*						如果ti和t2均小于n，
*									将t的成员is_a设置为真，将i输入t，获取返回值new_i;再将t的成员is_a设置为假，将j输入t，获取返回值new_j
*									如果t1<=t2,得到type_pair_index = t1 + t2*(t2+1)/2，否则type_pair_index = t2 + t1*(t1+1)/2
*									创建interacting_pair对象，用type_pair_index，new_i, new_j初始化，并将该对象添加到other_pairs尾部
*
*		5.检查本model对象成员minus_forces的大小与coords的大小是否相等，不相等则报错，终止程序的执行
*		6.检查输入model对象成员minus_forces的大小与coords的大小是否相等，不相等则报错，终止程序的执行
* 
*		7.将输入model类对象成员internal_coords容器的元素插入到本model类对象成员internal_coordss容器适当位置（根据coords_append函数的更新算法对internal_coords进行更新）
*		8.将输入model类对象成员coords         容器的元素插入到本model类对象成员coords          容器适当位置（根据coords_append函数的更新算法对coords         进行更新）
*		9.将输入model类对象成员minus_forces   容器的元素插入到本model类对象成员minus_forces    容器适当位置（根据coords_append函数的更新算法对minus_forces   进行更新）
*		
*		10.将输入model类对象成员ligands     容器的元素添加到本model类对象成员ligands     容器尾部，并根据model.append函数的更新算法对本model类对象成员ligands     容器的元素进行更新
*		11.将输入model类对象成员flex        容器的元素添加到本model类对象成员flex        容器尾部，并根据model.append函数的更新算法对本model类对象成员flex        容器的元素进行更新
*		12.将输入model类对象成员flex_context容器的元素添加到本model类对象成员flex_context容器尾部，并根据model.append函数的更新算法对本model类对象成员flex_context容器的元素进行更新
*       
*		13.将输入model类对象成员grid_atoms容器的元素添加到本model类对象成员grid_atoms容器尾部，并根据model.append函数的更新算法对本model类对象成员grid_atoms容器的元素进行更新
*		14.将输入model类对象成员atoms容器的元素插入到本model类对象成员atoms容器适当位置（根据coords_append函数的更新算法对atoms进行更新）
*		15.将输入对象的m_num_movable_atoms追加到本对象的m_num_movable_atoms
*/


void model::append(const model& m) {
	VINA_CHECK(atom_typing_used() == m.atom_typing_used());

	appender t(*this, m);

	t.append(other_pairs, m.other_pairs);

	VINA_FOR_IN(i, atoms)
		VINA_FOR_IN(j, m.atoms) {
			if(i >= m_num_movable_atoms && j >= m.m_num_movable_atoms) continue; // no need for inflex-inflex interactions

			const atom& a =   atoms[i];
			const atom& b = m.atoms[j];

			sz t1 = a.get(atom_typing_used());
			sz t2 = b.get(atom_typing_used());
			sz n = num_atom_types(atom_typing_used());

			if(t1 < n && t2 < n) {
				t.is_a =  true;
				sz new_i = t(i);
				t.is_a = false;
				sz new_j = t(j);
				sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
				other_pairs.push_back(interacting_pair(type_pair_index, new_i, new_j));
			}
		}

	VINA_CHECK(  minus_forces.size() ==   coords.size());
	VINA_CHECK(m.minus_forces.size() == m.coords.size());

	t.coords_append(internal_coords, m.internal_coords);
	t.coords_append(         coords, m         .coords);
	t.coords_append(   minus_forces, m   .minus_forces); // for now, minus_forces.size() == coords.size() (includes inflex)

	t.append(ligands,         m.ligands);
	t.append(flex,            m.flex);
	t.append(flex_context,    m.flex_context);

	t       .append(grid_atoms, m.grid_atoms);
	t.coords_append(     atoms, m     .atoms);

	m_num_movable_atoms += m.m_num_movable_atoms;
}

///////////////////  end  MODEL::APPEND /////////////////////////


/////////////////// begin MODEL::INITIALIZE /////////////////////////


/*
* 定义model的成员函数sz_to_atom_index
* input:unsigned int类型i
* ouput:atom_index类
* 作用：判断i是否小于model成员grid_atoms的大小，是则返回一个atom_index类，并且初始化该类的成员i为i，in_grid为true
*                                               否则返回一个atom_index类，并且初始化该类的成员i为(i-model成员grid_atoms的大小)，in_grid为false
*/
atom_index model::sz_to_atom_index(sz i) const {
	if(i < grid_atoms.size()) return atom_index(i                    ,  true);
	else                      return atom_index(i - grid_atoms.size(), false);
}

/*
* 定义model的成员函数distance_type_between
* input:distance_type_matrix类常量引用mobility，atom_index类常量引用i，atom_index类常量引用j
* output:枚举类型distance_type 
*        按以下顺序选择输出
*        1.i的成员in_grid为真并且j的成员in_grid为真返回DISTANCE_FIXED
*        2.只有i的成员in_grid为真，则根据j的成员i是否小于m_num_movable_atoms，是返回DISTANCE_VARIABLE，否返回DISTANCE_FIXED
*		 3.只有j的成员in_grid为真，则根据i的成员i是否小于m_num_movable_atoms，是返回DISTANCE_VARIABLE，否返回DISTANCE_FIXED
*        4.如果i的成员i小于model的成员atoms容器的大小或者j的成员i小于model的成员atoms容器的大小有一个为假，则报错，终止程序执行
*        5.如果i的成员i等于j的成员i，返回DISTANCE_FIXED
*        6.如果i的成员i小于j的成员i，将a,b分别输入mobility的调用运算符，返回mobility成员m_data的第（a + b*(b-1)/2）元素
*        7.如果i的成员i大于j的成员i，将b,a分别输入mobility的调用运算符，返回mobility成员m_data的第（b + a*(a-1)/2）元素
*/
distance_type model::distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const {
	if(i.in_grid && j.in_grid) return DISTANCE_FIXED;
	if(i.in_grid) return (j.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	if(j.in_grid) return (i.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	assert(!i.in_grid);
	assert(!j.in_grid);
	assert(i.i < atoms.size());
	assert(j.i < atoms.size());
	sz a = i.i;
	sz b = j.i;
	if(a == b) return DISTANCE_FIXED;
	return (a < b) ? mobility(a, b) : mobility(b, a);//////////////////////////////////////不知道什么意思
}

/*
* 定义model的成员函数atom_coords
* input:atom_index类引用i
* ouput:vec类常量引用
* 作用：判断i的成员in_grid是否为真，是则返回model的成员grid_atoms容器的第（i的成员i的数值）个元素的成员coords矢量
*                                   否则返回model的成员coords容器的第（i的成员i的数值）个元素
*/
const vec& model::atom_coords(const atom_index& i) const {
	return i.in_grid ? grid_atoms[i.i].coords : coords[i.i];
}

/*
* 定义model的成员函数distance_sqr_between
* input:atom_index类常量引用a,b
* ouput:double类型数
* 作用：将a,b分别传入atom_coords函数，返回两个三维矢量，返回这两个vec中逐个data数组元素差的平方和                                   
*/
fl model::distance_sqr_between(const atom_index& a, const atom_index& b) const {
	return vec_distance_sqr(atom_coords(a), atom_coords(b));
}


/*
* 结构体bond_less
* 定义调用操作符（）
* input:bond类常量引用a,b
* ouput:布尔值
* 作用：比较a的成员connected_atom_index的成员i和b的成员connected_atom_index的成员i,前者大返回false，后者大返回True
*/
struct bond_less { // FIXME rm!?
	bool operator()(const bond& a, const bond& b) const {
		return a.connected_atom_index.i < b.connected_atom_index.i;
	}
};

/*
* 定义model的成员函数atom_exists_between
* input:distance_type_matrix类常量引用mobility，atom_index类常量引用a，atom_index类常量引用b,unsigend int类型vector常量引用relevant_atoms
* output:bool值，寻找满足条件的相关原子是否存在
* 作用：1.根据a，b的in_grid是否为真，有 r2 = (grid_atoms[a.i].coords[0]-grid_atoms[b.i].coords[0])^2+(grid_atoms[a.i].coords[1]-grid_atoms[b.i].coords[1])^2+(grid_atoms[a.i].coords[2]-grid_atoms[b.i].coords[2])^2------------------------a.in_grid = TRUE  , b.in_grid = TRUE
*									    
*									    r2 = (coords[a.i][0]-grid_atoms[b.i].coords[0])^2+(coords[a.i][1]-grid_atoms[b.i].coords[1])^2+(coords[a.i][2]-grid_atoms[b.i].coords[2])^2---------------------------------------------------------a.in_grid = FALSE , b.in_grid = TRUE 
* 
*										r2 = (grid_atoms[a.i].coords[0]-coords[b.i][0])^2+(grid_atoms[a.i].coords[1]-coords[b.i][1])^2+(grid_atoms[a.i].coords[2]-coords[b.i][2])^2---------------------------------------------------------a.in_grid = TRUE  , b.in_grid = FALSE
* 
*										r2 = (coords[a.i][0]-coords[b.i][0])^2+(coords[a.i][1]-coords[b.i][1])^2+(coords[a.i][2]-coords[b.i][2])^2------------------------------------------------------------------------------------------a.in_grid = FALSE , b.in_grid = FALSE
*		2.遍历relevant_atoms内的元素
*				将relevant_atoms的第relevant_atoms_i元素赋值给i
*				判断i是否小于model成员grid_atoms的大小，是则返回一个atom_index类c，并且初始化该类的成员i为i，in_grid为true;否则返回一个atom_index类c，并且初始化该类的成员i为(i-model成员grid_atoms的大小)，in_grid为false
*		 		如果 a等于c或者b等于c 直接进行下一次循环                                        
*				将mobility, a, c输入distance_type_between函数返回得到distance_type枚举类型ac
*				将mobility, b, c输入distance_type_between函数返回得到distance_type枚举类型bc
*				如果ac不等于DISTANCE_VARIABLE，并且bc不等于DISTANCE_VARIABLE，并且distance_sqr_between(a, c)小于r2，并且distance_sqr_between(b, c)小于r2返回真，退出函数
*				其中distance_sqr_between(a, c)根据a，c的in_grid是否为真，有 r2 = (grid_atoms[a.i].coords[0]-grid_atoms[c.i].coords[0])^2+(grid_atoms[a.i].coords[1]-grid_atoms[c.i].coords[1])^2+(grid_atoms[a.i].coords[2]-grid_atoms[c.i].coords[2])^2------------------------a.in_grid = TRUE  , c.in_grid = TRUE
* 											  						    																																																							
* 											  						        r2 = (coords[a.i][0]-grid_atoms[c.i].coords[0])^2+(coords[a.i][1]-grid_atoms[c.i].coords[1])^2+(coords[a.i][2]-grid_atoms[c.i].coords[2])^2---------------------------------------------------------a.in_grid = FALSE , c.in_grid = TRUE 
* 																		    																																																						
* 											  							    r2 = (grid_atoms[a.i].coords[0]-coords[c.i][0])^2+(grid_atoms[a.i].coords[1]-coords[c.i][1])^2+(grid_atoms[a.i].coords[2]-coords[c.i][2])^2---------------------------------------------------------a.in_grid = TRUE  , c.in_grid = FALSE
* 																		    																																																						
* 											  							    r2 = (coords[a.i][0]-coords[c.i][0])^2+(coords[a.i][1]-coords[c.i][1])^2+(coords[a.i][2]-coords[c.i][2])^2------------------------------------------------------------------------------------------a.in_grid = FALSE , c.in_grid = FALSE
* 
*				其中distance_sqr_between(b, c)根据a，c的in_grid是否为真，有 r2 = (grid_atoms[b.i].coords[0]-grid_atoms[c.i].coords[0])^2+(grid_atoms[b.i].coords[1]-grid_atoms[c.i].coords[1])^2+(grid_atoms[b.i].coords[2]-grid_atoms[c.i].coords[2])^2------------------------b.in_grid = TRUE  , c.in_grid = TRUE
* 											  						    																																																							
* 											  						        r2 = (coords[b.i][0]-grid_atoms[c.i].coords[0])^2+(coords[b.i][1]-grid_atoms[c.i].coords[1])^2+(coords[b.i][2]-grid_atoms[c.i].coords[2])^2---------------------------------------------------------b.in_grid = FALSE , c.in_grid = TRUE 
* 																		    																																																						
* 											  							    r2 = (grid_atoms[b.i].coords[0]-coords[c.i][0])^2+(grid_atoms[b.i].coords[1]-coords[c.i][1])^2+(grid_atoms[b.i].coords[2]-coords[c.i][2])^2---------------------------------------------------------b.in_grid = TRUE  , c.in_grid = FALSE
* 																		    																																																						
* 											  							    r2 = (coords[b.i][0]-coords[c.i][0])^2+(coords[a.i][1]-coords[b.i][1])^2+(coords[b.i][2]-coords[c.i][2])^2------------------------------------------------------------------------------------------b.in_grid = FALSE , c.in_grid = FALSE
*		3.以上遍历找不到相关的原子，返回假
* 
* 
* 
* 
* 
*        
*/
bool model::atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const { // there is an atom closer to both a and b then they are to each other and immobile relative to them
	fl r2 = distance_sqr_between(a, b);
	VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
		sz i = relevant_atoms[relevant_atoms_i];
		atom_index c = sz_to_atom_index(i);
		if(a == c || b == c) continue;
		distance_type ac = distance_type_between(mobility, a, c);
		distance_type bc = distance_type_between(mobility, b, c);
		if(ac != DISTANCE_VARIABLE &&
		   bc != DISTANCE_VARIABLE &&
		   distance_sqr_between(a, c) < r2 &&
		   distance_sqr_between(b, c) < r2)
			return true;
	}
	return false;
}

/*
* 结构体beads
* 成员如下：
*			1.double类型radius_sqr
*			2.vector类data，元素为pair类，pair类的first数据成员为vec类，second数据成员为unsignd int类型vector
*			3.构造函数beads：input:unsigned int类型reserve_size,double数radius_sqr
*							 作用：初始化beads类的radius_sqr为输入radius_sqr，data分配能容纳至少reserve_size个元素的内存空间
*			4.函数add：input:unsigned int类型index,vec类常量引用coords
* 					   作用：1.遍历vector类data
*									如果(coords)与(data的第i个元素的first数据成员)的逐个data数组元素差的平方<radius_sqr
*									往data的第i个元素的second数据成员尾部添加元素index
*									退出函数
*							 2.第一步没有找对合适元素的情况下，创建一个pair类tmp,其first数据成员为vec类，second数据成员为unsignd int类型vector
*							   将tmp的first数据成员赋值为coords，second数据成员尾部添加index，将tmp添加到data尾部
*/
struct beads {
	fl radius_sqr;
	std::vector<std::pair<vec, szv> > data;
	beads(sz reserve_size, fl radius_sqr_) : radius_sqr(radius_sqr_) { data.reserve(reserve_size); }
	void add(sz index, const vec& coords) {
		VINA_FOR_IN(i, data) {
			if(vec_distance_sqr(coords, data[i].first) < radius_sqr) {
				data[i].second.push_back(index);
				return;
			}
		}
		// not found
		std::pair<vec, szv> tmp;
		tmp.first = coords;
		tmp.second.push_back(index);
		data.push_back(tmp);
	}
};
/*
* 定义model的成员函数assign_bonds
* input:distance_type_matrix类常量引用mobility
* ouput:空
* 作用：1.定义一个double类型值bond_length_allowance_factor = 1.1 粘结长度允许系数
*		  n = grid_atoms 和 atoms 容器大小之和
*		  bead_radius = 1.5 
*		2.定义一个beads类对象beads_instance，初始化beads_instance的radius_sqr = 225 ,data的预留空间为n
*		3.i从0循环到n
*				定义一个atom_index类i_atom_index，当0 <= i < grid_atoms.size()时，初始化该类的成员i为i，in_grid为true；当i >= grid_atoms.size()时，初始化该类的成员i为(i-model成员grid_atoms的大小)，in_grid为false
*				当0 <= i < grid_atoms.size()时，将i和grid_atoms[i_atom_index.i].coords输入beads_instance成员函数add；当i >= grid_atoms.size()时，将i和coords[i_atom_index.i]输入beads_instance成员函数add
* 
*				当0 <= i < grid_atoms.size()时，在beads_instance.data[i](此i是在遍历data，非前面i从0循环到n)内寻找第一个满足[data[i].first[0] - grid_atoms[i_atom_index.i].coords[0]]^2+[data[i].first[1] - grid_atoms[i_atom_index.i].coords[1]]^2+[data[i].first[2] - grid_atoms[i_atom_index.i].coords[2]]^2<225的beads_instance.data的元素并在data[i].second的尾部存入i
*				如果data不存在满足条件的元素，定义一个pair类tmp，first数据成员存入grid_atoms[i_atom_index.i].coords，second数据成员尾部添加i，将tmp存入beads_instance.data尾部
* 
*				当i >= grid_atoms.size()时，在beads_instance.data[i](此i是在遍历data，非前面i从0循环到n)内寻找第一个满足[data[i].first[0] - coords[i_atom_index.i][0]]^2+[data[i].first[1] - coords[i_atom_index.i][1]]^2+[data[i].first[2] - coords[i_atom_index.i][2]]^2<225的beads_instance.data的元素并在data[i].second的尾部存入i
*				如果data不存在满足条件的元素，定义一个pair类tmp，first数据成员存入coords[i_atom_index.i]，second数据成员尾部添加i，将tmp存入beads_instance.data尾部
* 
*		4.i从0循环到n
*				定义一个atom_index类i_atom_index，当0 <= i < grid_atoms.size()时，初始化该类的成员i为i，in_grid为true；当i >= grid_atoms.size()时，初始化该类的成员i为(i-model成员grid_atoms的大小)，in_grid为false
*				当0 <= i < grid_atoms.size()时，将grid_atoms[i_atom_index.i].coords绑定到i_atom_coords；当i >= grid_atoms.size()时，将coords[i_atom_index.i]绑定到i_atom_coords
*				当0 <= i < grid_atoms.size()时，将grid_atoms[i_atom_index.i]绑定到i_atom；当i >= grid_atoms.size()时，将atoms[i_atom_index.i]绑定到i_atom
*				
*				定义一个double常量类型值max_covalent_r = 1.74----------->atom_kind_data[i].covalent_radius最大值     共价半径
*				定义一个double类型值i_atom_covalent_radius = max_covalent_r = 1.74
*				如果i_atom的ad小于20，i_atom_covalent_radius赋值为atom_kind_data（i_atom.ad）的covalent_radius
* 
*				定义一个unsigned int类型vector -->relevant_atoms
*				定义一个double常量类型值bead_cutoff_sqr = (15 + 1.1*(i_atom_covalent_radius+1.74))^2
* 
*				b从0循环到beads_instance.data的大小-1（遍历元素）
*						当（beads_instance.data[b].first[0] - i_atom_coords[0])^2 + （beads_instance.data[b].first[1] - i_atom_coords[1])^2+（beads_instance.data[b].first[2] - i_atom_coords[2])^2> bead_cutoff_sqr时，直接继续下一个b遍历
*						将beads_instance.data[b].second绑定到bead_elements（unsigned int类型vector）
*						b从bead_elements_i循环到bead_elements的大小-1（遍历元素）
*								定义j = bead_elements的第bead_elements_i个元素
*								定义一个atom_index类j_atom_index，当0 <= j < grid_atoms.size()时，初始化该类的成员i为j，in_grid为true；当i >= grid_atoms.size()时，初始化该类的成员i为(j-model成员grid_atoms的大小)，in_grid为false
*								当0 <= j < grid_atoms.size()时，将grid_atoms[j_atom_index.i]绑定到j_atom；当j >= grid_atoms.size()时，将atoms[j_atom_index.i]绑定到j_atom
*								计算i_atom和j_atom的共价半径和bond_length
*								将mobility, i_atom_index, j_atom_index输入distance_type_between函数，获取i_atom_index和j_atom_index的距离类型dt
*								如果dt != DISTANCE_VARIABLE并且i != j，
*										当0 <= i < grid_atoms.size()并且0 <= j < grid_atoms.size()时，令r2 = [grid_atoms[i_atom_index.i].coords[0]-grid_atoms[j_atom_index.i].coords[0]]^2+[grid_atoms[i_atom_index.i].coords[1]-grid_atoms[j_atom_index.i].coords[1]]^2+[grid_atoms[i_atom_index.i].coords[2]-grid_atoms[j_atom_index.i].coords[2]]^2
*										当0 <= i < grid_atoms.size()并且j > grid_atoms.size()时，令r2 = [grid_atoms[i_atom_index.i].coords[0]-coords[j_atom_index.i][0]]^2+[grid_atoms[i_atom_index.i].coords[1]-coords[j_atom_index.i][1]]^2+[grid_atoms[i_atom_index.i].coords[2]-coords[j_atom_index.i][2]]^2
*										当i > grid_atoms.size()并且0 <= j < grid_atoms.size()时，令r2 = [coords[i_atom_index.i][0]-grid_atoms[j_atom_index.i].coords[0]]^2+[coords[i_atom_index.i][1]-grid_atoms[j_atom_index.i].coords[1]]^2+[coords[i_atom_index.i][2]-grid_atoms[j_atom_index.i].coords[2]]^2
*										当i > grid_atoms.size()并且j > grid_atoms.size()时，令r2 = [coords[i_atom_index.i][0]-coords[j_atom_index.i][0]]^2+[coords[i_atom_index.i][1]-coords[j_atom_index.i][1]]^2+[coords[i_atom_index.i][2]-coords[j_atom_index.i][2]]^2
*										如果r2 < (1.1*(i_atom_covalent_radius+1.74))^2
*										relevant_atoms尾部添加j
*				
*				relevant_atoms_i从0循环到relevant_atoms的大小-1（遍历元素）
*						定义j等于relevant_atoms的relevant_atoms_i个元素
*						如果j <= i,直接继续下一个relevant_atoms_i遍历
*						定义一个atom_index类j_atom_index，当0 <= j < grid_atoms.size()时，初始化该类的成员i为j，in_grid为true；当i >= grid_atoms.size()时，初始化该类的成员i为(j-model成员grid_atoms的大小)，in_grid为false
* 						当0 <= j < grid_atoms.size()时，将grid_atoms[j_atom_index.i]绑定到j_atom；当j >= grid_atoms.size()时，将atoms[j_atom_index.i]绑定到j_atom
*						计算i_atom和j_atom的共价半径和bond_length
*						将mobility, i_atom_index, j_atom_index输入distance_type_between函数，获取i_atom_index和j_atom_index的距离类型dt
*						当0 <= i < grid_atoms.size()并且0 <= j < grid_atoms.size()时，令r2 = [grid_atoms[i_atom_index.i].coords[0]-grid_atoms[j_atom_index.i].coords[0]]^2+[grid_atoms[i_atom_index.i].coords[1]-grid_atoms[j_atom_index.i].coords[1]]^2+[grid_atoms[i_atom_index.i].coords[2]-grid_atoms[j_atom_index.i].coords[2]]^2								
* 						当0 <= i < grid_atoms.size()并且j > grid_atoms.size()时，令r2 = [grid_atoms[i_atom_index.i].coords[0]-coords[j_atom_index.i][0]]^2+[grid_atoms[i_atom_index.i].coords[1]-coords[j_atom_index.i][1]]^2+[grid_atoms[i_atom_index.i].coords[2]-coords[j_atom_index.i][2]]^2
* 						当i > grid_atoms.size()并且0 <= j < grid_atoms.size()时，令r2 = [coords[i_atom_index.i][0]-grid_atoms[j_atom_index.i].coords[0]]^2+[coords[i_atom_index.i][1]-grid_atoms[j_atom_index.i].coords[1]]^2+[coords[i_atom_index.i][2]-grid_atoms[j_atom_index.i].coords[2]]^2
* 						当i > grid_atoms.size()并且j > grid_atoms.size()时，令r2 = [coords[i_atom_index.i][0]-coords[j_atom_index.i][0]]^2+[coords[i_atom_index.i][1]-coords[j_atom_index.i][1]]^2+[coords[i_atom_index.i][2]-coords[j_atom_index.i][2]]^2
*						如果r2 < (1.1*bond_length)^2并且mobility, i_atom_index, j_atom_index, relevant_atoms 代入atom_exists_between不为真（relevant_atoms中不存在与i_atom_index和j_atom_index满足一定距离关系的原子）
*								如果dt 等于 DISTANCE_ROTOR，令rotatable等于1，否则为0
*								length 等于 r2的开方数
*								在i_atom.bonds尾部添加一个bond元素,并用j_atom_index, length, rotatable对其进行初始化
*								在j_atom.bonds尾部添加一个bond元素,并用j_atom_index, length, rotatable对其进行初始化
*/
void model::assign_bonds(const distance_type_matrix& mobility) { // assign bonds based on relative mobility, distance and covalent length
	const fl bond_length_allowance_factor = 1.1;
	sz n = grid_atoms.size() + atoms.size();

	// construct beads
	const fl bead_radius = 15;
	beads beads_instance(n, sqr(bead_radius));
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		beads_instance.add(i, atom_coords(i_atom_index));
	}
	// assign bonds
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		const vec& i_atom_coords = atom_coords(i_atom_index);
		atom& i_atom = get_atom(i_atom_index);

		const fl max_covalent_r = max_covalent_radius(); // FIXME mv to atom_constants
		fl i_atom_covalent_radius = max_covalent_r;
		if(i_atom.ad < AD_TYPE_SIZE)
			i_atom_covalent_radius = ad_type_property(i_atom.ad).covalent_radius;

		//find relevant atoms
		szv relevant_atoms;
		const fl bead_cutoff_sqr = sqr(bead_radius + bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r));
		VINA_FOR_IN(b, beads_instance.data) {
			if(vec_distance_sqr(beads_instance.data[b].first, i_atom_coords) > bead_cutoff_sqr) continue;
			const szv& bead_elements = beads_instance.data[b].second;
			VINA_FOR_IN(bead_elements_i, bead_elements) {
				sz j = bead_elements[bead_elements_i];
				atom_index j_atom_index = sz_to_atom_index(j);
				atom& j_atom = get_atom(j_atom_index);
				const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
				distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
				if(dt != DISTANCE_VARIABLE && i != j) {
					fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
					//if(r2 < sqr(bond_length_allowance_factor * bond_length))
					if(r2 < sqr(bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r)))
						relevant_atoms.push_back(j);
				}
			}
		}
		// find bonded atoms
		VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
			sz j = relevant_atoms[relevant_atoms_i];
			if(j <= i) continue; // already considered
			atom_index j_atom_index = sz_to_atom_index(j);
			atom& j_atom = get_atom(j_atom_index);
			const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
			distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
			fl r2 = distance_sqr_between(i_atom_index, j_atom_index);///////////////////////////////////////////////////////
			if(r2 < sqr(bond_length_allowance_factor * bond_length) && !atom_exists_between(mobility, i_atom_index, j_atom_index, relevant_atoms)) {
				bool rotatable = (dt == DISTANCE_ROTOR);
				fl length = std::sqrt(r2);
				i_atom.bonds.push_back(bond(j_atom_index, length, rotatable));
				j_atom.bonds.push_back(bond(i_atom_index, length, rotatable));
			}

		}
	}
}


/*
* 定义model的成员函数bonded_to_HD
* input:atom类常量引用a
* ouput:bool值
* 作用：1.遍历a的成员bonds的元素
*       2.定义一个临时bond类对象常量引用b，将a的成员bonds的第i个元素绑定到b上
*       3.将b的成员connected_atom_index输入get_atom函数，返回一个atom类对象，判断该对象的成员ad是否等于12，如果是则退出函数返回true，否则继续循环
*       4.循环结束，返回false
*/
bool model::bonded_to_HD(const atom& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).ad == AD_TYPE_HD) 
			return true;
	}
	return false;
}


/*
* 定义model的成员函数bonded_to_heteroatom
* input:atom类常量引用a
* ouput:bool值
* 作用：1.遍历a的成员bonds的元素
*       2.定义一个临时bond类对象常量引用b，将a的成员bonds的第i个元素绑定到a上
*       3.将b的成员connected_atom_index输入get_atom函数，返回一个atom类对象，判断该对象的成员is_heteroatom函数返回值是否为真，如果是则退出函数返回true，否则继续循环
*       4.循环结束，返回false
*/
bool model::bonded_to_heteroatom(const atom& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).is_heteroatom())
			return true;
	}
	return false;
}

/*
* 定义model的成员函数assign_types
* input:无
* ouput:无
* 作用：i从0循环到（model成员grid_atoms的大小和atoms的大小之和-1）
*				将i输入到model成员函数sz_to_atom_index，返回的值赋值到ai上，当0<=i<=(grid_atoms的大小-1)时，ai为atom_index(i,true)，否则为atom_index(i - grid_atoms.size(), false)
*				将ai输入model成员函数get_atom，ai为atom_index(i,true)时，获得model成员grid_atoms的第（ai的成员i）个元素，ai为atom_index(i - grid_atoms.size(), false)时，获得model成员atoms的第i（ai的成员i）个元素，返回的值绑定到a上
*				给a的成员el赋值，若a的成员ad==20且a的成员xs==16 则el=10，否则el=ad输入ad_type_to_el_type函数后的返回值
*				将a的成员xs绑定到x上
*				定义一个bool值acceptor，如果a的成员ad等于9或者10，acceptor为true,否则为false
*				定义一个bool值donor_NorO，如果a的成员el等于10或者将a输入model成员函数bonded_to_HD返回为真，donor_NorO为true,否则为false
*				根据a的成员el大小为x赋值
*							当el =  EL_TYPE_H    =  0时，不执行任何操作
* 							当el =  EL_TYPE_C    =  1时，将a输入model成员函数bonded_to_heteroatom，返回真x赋值为1，否则x赋值为0
* 							当el =  EL_TYPE_N    =  2时，如果acceptor和donor_NorO均为真，x赋值为5，只有acceptor为真x赋值为4，只有donor_NorO为真x赋值为3，都为假x赋值为2
* 							当el =  EL_TYPE_O    =  3时，如果acceptor和donor_NorO均为真，x赋值为9，只有acceptor为真x赋值为8，只有donor_NorO为真x赋值为7，都为假x赋值为6
* 							当el =  EL_TYPE_S    =  4时，x赋值为XS_TYPE_S_P   = 10
* 							当el =  EL_TYPE_P    =  5时，x赋值为XS_TYPE_P_P   = 11
* 							当el =  EL_TYPE_F    =  6时，x赋值为XS_TYPE_F_H   = 12
* 							当el =  EL_TYPE_Cl   =  7时，x赋值为XS_TYPE_Cl_H  = 13
* 							当el =  EL_TYPE_Br   =  8时，x赋值为XS_TYPE_Br_H  = 14
* 							当el =  EL_TYPE_I    =  9时，x赋值为XS_TYPE_I_H   = 15
* 							当el =  EL_TYPE_Met  = 10时，x赋值为XS_TYPE_Met_D = 16
* 							当el =  EL_TYPE_SIZE = 11时，不执行任何操作
*							均不是，报错，终止程序执行
*
*/
void model::assign_types() {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		atom& a = get_atom(ai);
		a.assign_el();
		sz& x = a.xs;

		bool acceptor   = (a.ad == AD_TYPE_OA || a.ad == AD_TYPE_NA); // X-Score forumaltion apparently ignores SA
		bool donor_NorO = (a.el == EL_TYPE_Met || bonded_to_HD(a));

		switch(a.el) {
			case EL_TYPE_H    : break;
			case EL_TYPE_C    : x = bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H; break;
			case EL_TYPE_N    : x = (acceptor && donor_NorO) ? XS_TYPE_N_DA : (acceptor ? XS_TYPE_N_A : (donor_NorO ? XS_TYPE_N_D : XS_TYPE_N_P)); break;
			case EL_TYPE_O    : x = (acceptor && donor_NorO) ? XS_TYPE_O_DA : (acceptor ? XS_TYPE_O_A : (donor_NorO ? XS_TYPE_O_D : XS_TYPE_O_P)); break;
			case EL_TYPE_S    : x = XS_TYPE_S_P; break;
			case EL_TYPE_P    : x = XS_TYPE_P_P; break;
			case EL_TYPE_F    : x = XS_TYPE_F_H; break;
			case EL_TYPE_Cl   : x = XS_TYPE_Cl_H; break;
			case EL_TYPE_Br   : x = XS_TYPE_Br_H; break;
			case EL_TYPE_I    : x = XS_TYPE_I_H; break;
			case EL_TYPE_Met  : x = XS_TYPE_Met_D; break;
			case EL_TYPE_SIZE : break;
			default: VINA_CHECK(false);
		}
	}
}

/*
* 定义model的成员函数find_ligand
* input:unsigned int类型a
* ouput:unsigned int类型值
* 作用：1.遍历model的成员ligands容器
*       2.找到第一个a在ligands容器元素的成员begin和end范围之内的元素，返回其下标
*       3.没有符合条件的则返回ligands容器的大小
*/
sz model::find_ligand(sz a) const {
	VINA_FOR_IN(i, ligands)
		if(a >= ligands[i].begin && a < ligands[i].end)
			return i;
	return ligands.size();
}

/*
* 定义model的成员函数bonded_to
* input:unsigned int类型a，n，unsigned int类型容器out
* ouput:unsigned int类型值
* 作用：1.判断a是否在容器out内，是则函数结束，否则继续执行余下语句
*       2.将a加入到容器out尾部
*       3.判断n是否大于0，是则继续执行余下语句，否则函数结束
*		4.遍历atoms容器第a个元素的成员bonds容器
*       5.将atoms容器第a个元素的成员bonds容器第i个元素绑定到b上
*       6.判断b的成员connected_atom_index类的成员in_grid是否为真，是则函数结束，否则继续执行余下语句
*       7.将b的成员connected_atom_index类的成员i，n-1,out输入bonded_to
* 目的：将符合上述规则的atmos的下标按连接顺序依次存入out尾部，规则见图三
*/

void model::bonded_to(sz a, sz n, szv& out) const {
	if(!has(out, a)) { // not found
		out.push_back(a);
		if(n > 0) 
			VINA_FOR_IN(i, atoms[a].bonds) {
				const bond& b = atoms[a].bonds[i];
				if(!b.connected_atom_index.in_grid)
					bonded_to(b.connected_atom_index.i, n-1, out);
			}
	}
}

/*
* 重载model的成员函数bonded_to
* input:unsigned int类型a，n
* ouput:unsigned int类型vector
* 作用：1.定义一个临时unsigned int类型vector对象tmp
*       2.调用bonded_to函数三输入版本，传入a,n,tmp
*       3.返回tmp
*/
szv model::bonded_to(sz a, sz n) const {
	szv tmp;
	bonded_to(a, n, tmp);
	return tmp;
}

/*
* 定义model的成员函数initialize_pairs
* input:distance_type_matrix类常量引用mobility
* ouput:无
* 作用：遍历model成员atmos内的元素
*			找到第一个i在ligands容器元素的成员begin和end范围之内的元素，返回其下标赋值给i_lig
*			将i,3带入函数bonded_to，返回atoms中与i连接满足一定规则的元素的下标组成的vector--->bonded_atoms，规则见图三
*			j从i+1循环到atoms.size()-1----------->atoms.size()为atmos容器的大小
*				如果i和j均大于model的m_num_movable_atoms（model的移动原子的数量），直接继续下一次循环
*				将i,j带入model中，返回mobility的成员m_data[i + j*(j-1)/2]，判断返回值是否为DISTANCE_VARIABLE并且j不在bonded_atoms，如果都满足，继续执行下面语句，否则继续下一次循环
*				将m_atom_typing_used代入atmos第i个元素的函数get中得到m_atom_typing_used对应的整数值赋给t1
*				将m_atom_typing_used代入atmos第j个元素的函数get中得到m_atom_typing_used对应的整数值赋给t2
*				将m_atom_typing_used代入num_atom_types中得到m_atom_typing_used对应的整数值赋给n
*				如果t1和t2均小于n,继续执行下面语句，否则继续下一次循环
*				将n, t1, t2代入triangular_matrix_index_permissive（三角矩阵索引），当t1<t2时，返回t1 + t2*(t2+1)/2，赋值给type_pair_index
*																				   当t1>t2时，返回t2 + t1*(t1+1)/2，赋值给type_pair_index
*				定义一个interacting_pair类ip，初始化初始化type_pair_index，a，b为type_pair_index, i, j
*				如果i_lig小于ligands（配体容器）大小并且第一个j在ligands容器元素的成员begin和end范围之内的元素，返回其下标等于i_lig
*				将ip存入ligands第i_lig个元素的成员pairs尾部，否则存入model的other_pairs尾部
* 
*/
void model::initialize_pairs(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, atoms) {
		sz i_lig = find_ligand(i);
		szv bonded_atoms = bonded_to(i, 3);
		VINA_RANGE(j, i+1, atoms.size()) {
			if(i >= m_num_movable_atoms && j >= m_num_movable_atoms) continue; // exclude inflex-inflex
			if(mobility(i, j) == DISTANCE_VARIABLE && !has(bonded_atoms, j)) {
				sz t1 = atoms[i].get  (atom_typing_used());
				sz t2 = atoms[j].get  (atom_typing_used());
				sz n  = num_atom_types(atom_typing_used());
				if(t1 < n && t2 < n) { // exclude, say, Hydrogens
					sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
					interacting_pair ip(type_pair_index, i, j);
					if(i_lig < ligands.size() && find_ligand(j) == i_lig)
						ligands[i_lig].pairs.push_back(ip);
					else
						other_pairs.push_back(ip);
				}
			}
		}
	}
}

/*
* 定义model的成员函数initialize
* input:distance_type_matrix类常量引用mobility
* ouput:无
* 作用：遍历model成员ligands内的元素
*			获取ligand[i]的成员node的begin，和其所有后代的成员node的begin，找到最小值赋值给ligand[i].begin
*			获取ligand[i]的成员node的end，  和其所有后代的成员node的end，  找到最大值赋值给ligand[i].end
*		将mobility输入assign_bonds，寻找键合原子
*		给grid_atoms的元素，atoms的元素分配类型
*		将mobilityinitialize_pairs，初始化配体的元素的pair
*/
void model::initialize(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, ligands)
		ligands[i].set_range();
	assign_bonds(mobility);
	assign_types();
	initialize_pairs(mobility);
}

///////////////////  end  MODEL::INITIALIZE /////////////////////////

/*
* 定义model的成员函数num_internal_pairs
* input:无
* ouput:unsigned int类型值
* 作用：1.定义一个临时unsigned int类型值tmp，初始化为0
*       2.求model成员ligands所有元素的成员pairs的大小之和，赋值给tmp并返回
*/
sz model::num_internal_pairs() const {
	sz tmp = 0;
	VINA_FOR_IN(i, ligands)
		tmp += ligands[i].pairs.size();
	return tmp;
}

/*
* 定义model的成员函数get_movable_atom_types，获取可移动的原子类型
* input:atom_type类定义的枚举类型t
* ouput:unsigned int型vector
* 作用：1.定义一个临时unsigned int型vector tmp
*       2.将atom_typing_used_转换成对应的（atom_typing_used_）_TYPE_SIZE赋值给n
*		3.从i=1循环到i=m_num_movable_atoms-1
*               （1）将model成员atoms的第i个元素绑定到a上
*               （2）将atom_typing_used_转换成a成员atom_typing_used_对应小写形式所代表的数，赋值给t
*               （3）如果t小于n并且t不在tmp内，将t存入tmp尾部
*       4.返回tmp
*/
szv model::get_movable_atom_types(atom_type::t atom_typing_used_) const {
	szv tmp;
	sz n = num_atom_types(atom_typing_used_);
	VINA_FOR(i, m_num_movable_atoms) {
		const atom& a = atoms[i];
		sz t = a.get(atom_typing_used_);
		if(t < n && !has(tmp, t))
			tmp.push_back(t);
	}
	return tmp;
}

/*
* 定义model的成员函数get_size
* input:无
* ouput:conf_size类
* 作用：1.获得ligands的扭转计数保存在临时tmp.ligands对象中
*       2.获得flex的扭转计数保存在临时tmp.flex对象中
*		3.返回tmp
*/
conf_size model::get_size() const {
	conf_size tmp;
	tmp.ligands = ligands.count_torsions();
	tmp.flex    = flex   .count_torsions();
	return tmp;
}

/*
* 定义model的成员函数get_initial_conf
* input:无
* ouput:conf类
* 作用：1.获得ligands的扭转计数保存在临时cs.ligands对象中
*       2.获得flex的扭转计数保存在临时cs.ligands对象中
*		3.返回tmp
*		4.定义一个conf对象tmp
*		tmp内		ligands向量的元素个数为cs中的ligands unint向量的元素个数;
*					flex向量的元素个数为cs中的flex uint向量的元素个数;
*
*					ligands[i].torsions double向量中元素个数为s.ligands[ligands.size]个，其值都为0；
*					flex[i].torsions double向量中元素个数为s.flex[flex.size]个，其值都为0；
*		5.	参数设为null
*		    tmp.ligands结构体向量中：
*		    data[0] = 0；data[1] = 0;data[2] = 0;
*		    orientation：一个实部为1，虚部都为0的四元数
*		    torsions向量元素都赋值为0；
*		    tmp.flex结构体向量中：
*		    torsions向量元素都赋值为0
*		6.遍历ligands内元素
*			将tmp.ligands[i].rigid.position设置为ligands[i].node.origin
*		7.返回tmp
*/
conf model::get_initial_conf() const { // torsions = 0, orientations = identity, ligand positions = current
	conf_size cs = get_size();
	conf tmp(cs);
	tmp.set_to_null();
	VINA_FOR_IN(i, ligands)
		tmp.ligands[i].rigid.position = ligands[i].node.get_origin();
	return tmp;
}

/*
* 定义model的成员函数movable_atoms_box
* input:double类型add_to_each_dimension，double类型granularity（粒度）
* ouput:grid_dims数组
* 作用：1.定义两个vec类corner1和corner2，并且均初始化类内的数组元素为0
*       2.i从0循环到m_num_movable_atoms（model的成员）
*				将model成员coords的第i个元素绑定到v上
*				j从0循环到2
*						当i=0时或者v的第j个元素小于corner1的第j个元素，将v的第j个元素赋值给corner1的第j个元素
*						当i=0时或者v的第j个元素大于corner1的第j个元素，将v的第j个元素赋值给corner1的第j个元素
*		 （获得model成员coords所有元素的数组中三个坐标的最小值和最大值，赋给corner1, corner2）		
*		3.corner1[0] = corner1[0]-add_to_each_dimension / 2, corner1[1] = corner1[1]-add_to_each_dimension / 2, corner1[2] = corner1[2]-add_to_each_dimension / 2
* 		  corner2[0] = corner2[0]+add_to_each_dimension / 2, corner2[1] = corner2[1]+add_to_each_dimension / 2, corner2[2] = corner2[2]+add_to_each_dimension / 2
*		4.定义一个grid_dims数组gd
*		  计算corner1和corner2的中心 center[0] = 0.5 * (corner2[0] + corner1[0]) center[1] = 0.5 * (corner2[1] + corner1[1]) center[2] = 0.5 * (corner2[2] + corner1[2])
*		  i从0循环到2
*				gd[i].n = 【(corner2[i] - corner1[i]) / granularity)向上取整】         gd[i].n---------gd[i]的成员n
*				real_span = granularity*gd[i].n，       gd[i].begin = center[i] - real_span/2，          gd[i].end = center[i] + real_span/2             gd[i].begin---------gd[i]的成员begin            gd[i].end---------gd[i]的成员end
*		5.返回gd
* 
* 
*/

grid_dims model::movable_atoms_box(fl add_to_each_dimension, fl granularity) const {
	vec corner1(0, 0, 0), corner2(0, 0, 0);
	VINA_FOR(i, num_movable_atoms()) {
		const vec& v = movable_coords(i);
		VINA_FOR_IN(j, v) {
			if(i == 0 || v[j] < corner1[j]) corner1[j] = v[j];
			if(i == 0 || v[j] > corner2[j]) corner2[j] = v[j];
		}
	}
	corner1 -= add_to_each_dimension / 2;
	corner2 += add_to_each_dimension / 2;

	grid_dims gd;
	{ // always doing this now FIXME ?
		vec center; center = 0.5 * (corner2 + corner1);
		VINA_FOR_IN(i, gd) {
			gd[i].n = sz(std::ceil((corner2[i] - corner1[i]) / granularity));
			fl real_span = granularity * gd[i].n;
			gd[i].begin = center[i] - real_span/2;
			gd[i].end = gd[i].begin + real_span;
		}
	}
	return gd;
}

/*
* 函数string_write_coord
* input:unsigned int类型i，double类型x，string类引用str
* ouput:无
* 作用：1.先检查输入i是否大于0，是则继续下面代码，否则报错，终止程序执行
*       2.i自减1
*		3.定义一个string类输出流out
*       4.对x进行格式设置保存在out中，格式为zzzz.zzz，且有三位有效小数
*		5.将保存在out中的string拷贝到str中，起始处位置为i（即输入i减去1）
*/
void string_write_coord(sz i, fl x, std::string& str) {
	VINA_CHECK(i > 0);
	--i;
	std::ostringstream out;
	//设定点输出格式,小数点后有6位数字,显示浮点数的小数点和小数位数
	out.setf(std::ios::fixed, std::ios::floatfield);
	out.setf(std::ios::showpoint);
	//setw(8)设置输出宽度，8个占位符，右对齐
	//设置浮点值的小数精度为3
	out << std::setw(8) << std::setprecision(3) << x; 
	//检查string输出流中string的大小是否等于8，不等于则报错，终止程序执行
	VINA_CHECK(out.str().size() == 8); 
	//检查str的大小是否大于i+8，小于等于则报错，终止程序执行
	VINA_CHECK(str.size() > i + 8);
	//将string类输出流out中的元素赋值到str中，起始处位置为i
	VINA_FOR(j, 8)
		str[i+j] = out.str()[j];
}

/*
* 函数coords_to_pdbqt_string
* input:vec类常量引用coords，string类常量引用str
* ouput:string类
* 作用：拷贝一个str为tmp，并将coords的三个元素转换成格式为zzzz.zzz，且有三位有效小数，依次将0，1，2位置的元素存入tmp起始位置为30，38，46位置处，返回tmp	
*/
std::string coords_to_pdbqt_string(const vec& coords, const std::string& str) {
	std::string tmp(str);
	string_write_coord(31, coords[0], tmp);
	string_write_coord(39, coords[1], tmp);
	string_write_coord(47, coords[2], tmp);
	return tmp;
}

/*
* 定义model的成员函数write_context
* input:parsed_line类型vector常量引用c，ofile类out
* ouput:void
* 作用：1.调用函数verify_bond_lengths验证连接长度
*		2.遍历c中的元素
*				将c的第i个元素的first数据成员绑定到str上
*				如果c的第i个元素的second数据成员不为空，返回c的第i个元素的second数据成员内保存的值a，获取coords的第a个元素，将其与str输入到coords_to_pdbqt_string函数，打印返回信息并换行，否则打印str并换行
*/
void model::write_context(const context& c, ofile& out) const {
	verify_bond_lengths();
	VINA_FOR_IN(i, c) {
		const std::string& str = c[i].first;
		if(c[i].second) {
			out << coords_to_pdbqt_string(coords[c[i].second.get()], str) << '\n';
		}
		else
			out << str << '\n';
	}
}

/*
* 定义model的成员函数seti
* input:conf类常量引用c
* ouput:void
* 作用：将model成员atoms, internal_coords, 和c成员ligands输入ligands的成员函数set_conf    (*this)[i].set_conf(atoms, coords, c[i]);
*       根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和internal_coords
*/
void model::seti(const conf& c) {
	ligands.set_conf(atoms, internal_coords, c.ligands);
}


/*
* 定义model的成员函数sete
* input:conf类常量引用c
* ouput:void
* 作用：遍历c.ligands的元素，将model成员internal_coords, coords, ligands[i].begin（ligand成员）, ligands[i].end（ligand成员）输入函数c.ligands[i].rigid.apply，获得坐标coords
* 根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords				
*/
void model::sete(const conf& c) {
	VINA_FOR_IN(i, ligands)
		c.ligands[i].rigid.apply(internal_coords, coords, ligands[i].begin, ligands[i].end);
	flex.set_conf(atoms, coords, c.flex);
}

/*
* 定义model的成员函数set
* input:conf类常量引用c
* ouput:void
* 根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
*/
void model::set         (const conf& c) {
	ligands.set_conf(atoms, coords, c.ligands);
	flex   .set_conf(atoms, coords, c.flex);
}


/*
* model的成员函数gyration_radius
* input:unsigned int类型ligand_number
* output:double类型数
* 作用：根据输入ligand_number，得到ligand第ligand_number个元素，获取它的成员begin,end，根据begin,end-1范围查找atoms的元素的成员el是否等于EL_TYPE_H
* 求解公式1，返回公式1结果的值
*/
fl model::gyration_radius(sz ligand_number) const {
	//ligand_number
	VINA_CHECK(ligand_number < ligands.size()); //检查ligand_number是否小于ligands容器的大小，不小于则报错，终止程序执行


	const ligand& lig = ligands[ligand_number]; //将ligands容器第ligand_number个元素绑定到lig上

	//定义内部变量double类型变量acc = 0，unsigned类型counter = 0
	fl acc = 0;
	unsigned counter = 0;
	//从lig的成员begin大小循环到lig的成员end大小-1
	VINA_RANGE(i, lig.begin, lig.end) {
		//如果model成员atoms的第i个元素的成员el不等于0
		//acc = acc+运算coords第i个元素和lig的成员node的成员origin逐个data数组元素差的平方
		//counter自加1
		if(atoms[i].el != EL_TYPE_H) { // only heavy atoms are used
			acc += vec_distance_sqr(coords[i], lig.node.get_origin()); // FIXME? check!
			++counter;
		}
	}
	//如果model成员atoms没有一个元素它的成员el等于0，返回0，否则返回见公式1
	return (counter > 0) ? std::sqrt(acc/counter) : 0;
}


/*
* model的成员函数eval_interacting_pairs   评估相互作用对
* input:precalculate类常量引用p，interacting_pairs类常量引用pairs，vecv类(vec类型vector)常量引用coords
* output:double类型数
* 作用：1.返回p的m_cutoff_sqr赋值给cutoff_sqr（截止面积），定义一个double类型e=0
*		2.遍历pairs的元素
*				将第i个元素绑定到ip上
*				计算coords[ip.a]（ip.a-->ip的成员a）和coords[ip.b]（ip.b-->ip的成员b）的平方距离(coords[ip.a][0]-coords[ip.b][0])^2+(coords[ip.a][1]-coords[ip.b][1])^2+(coords[ip.a][2]-coords[ip.b][2])^2，赋值给r2
*				判断r2是否小于cutoff_sqr，是则继续执行下面语句，否则继续下一次循环
*				将ip.type_pair_index和r2代入p的eval_fast，返回p.data.m_data[ip.type_pair_index].fast[sz(factor * r2)](其中sz为强制类型转换，factor为p.data.m_data[ip.type_pair_index]的数据成员)赋值给tmp
*				当v < epsilon_fl（运行编译程序的计算机所能识别的最小非零浮点数）时，tmp = 0；否则tmp = tmp*(v / (v + tmp))
*				e累加tmp
*		3.返回e
* 
*/
fl eval_interacting_pairs(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords) { // clean up
	const fl cutoff_sqr = p.cutoff_sqr();
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		fl r2 = vec_distance_sqr(coords[ip.a], coords[ip.b]);
		if(r2 < cutoff_sqr) {
			fl tmp = p.eval_fast(ip.type_pair_index, r2);
			curl(tmp, v);
			e += tmp;
		}
	}
	return e;
}

/*
* model的成员函数eval_interacting_pairs_deriv   评估相互作用对导数
* input:precalculate类常量引用p，double类型数v，interacting_pairs类常量引用pairs，vecv类(vec类型vector)常量引用coords，vecv类(vec类型vector)引用coords
* output:double类型数
* 作用：1.返回p的m_cutoff_sqr赋值给cutoff_sqr（截止面积），定义一个double类型e=0
*		2.遍历pairs的元素
*				将第i个元素绑定到ip上
*				计算coords[ip.a]（ip.a-->ip的成员a）和coords[ip.b]（ip.b-->ip的成员b）的差r[0] = coords[ip.a][0]-coords[ip.b][0]，r[1] = coords[ip.a][1]-coords[ip.b][1]，r[2] = coords[ip.a][2]-coords[ip.b][2]
*				计算r2[0] = r[0]^2，r2[1] = r[1]^2，r2[2] = r[2]^2
*				判断r2是否小于cutoff_sqr，是则继续执行下面语句，否则继续下一次循环
* 
*				将ip.type_pair_index和r2代入p的eval_deriv，返回返回pair类（firs和second数据成员均为double）pr，
*				pr(p1.first  + （factor * r2 -  sz(factor * r2)） * (p2.first  - p1.first), p1.second + （factor * r2 -  sz(factor * r2)）* (p2.second - p1.second))，其中p1 = smooth[sz(factor * r2)]，p2 = smooth[sz(factor * r2)+1]均为pair类
*				(其中sz为强制类型转换，factor为p.data.m_data[ip.type_pair_index]的数据成员，smooth为p.data.m_data[ip.type_pair_index]的数据成员（pair类型vector）)赋值给pr
* 
*				force[0] = tmp.second * r[0]，force[1] = tmp.second * r[1]，force[2] = tmp.second * r[2]         tmp.second------------>tmp的second数据成员
*				当tmp.first > 0 并且 v <  epsilon_fl（运行编译程序的计算机所能识别的最小非零浮点数）时，tmp.first = 0，force[0] = 0，force[1] = 0，force[2] = 0；
*				当tmp.first > 0 并且 v >= epsilon_fl（运行编译程序的计算机所能识别的最小非零浮点数）时，force[0] = force[0]*((v / (v + tmp.first))^2)，force[1] = force[1]*((v / (v + tmp.first))^2)，force[2] = force[2]*((v / (v + tmp.first))^2),
*																									   tmp.first = tmp.first*(v / (v + tmp.first))	
*				e累加tmp.first
*				forces[ip.a][0] = forces[ip.a][0] - force[0],forces[ip.a][1] = forces[ip.a][1] - force[1],forces[ip.a][2] = forces[ip.a][2] - force[2]
*				forces[ip.b][0] = forces[ip.b][0] + force[0],forces[ip.b][1] = forces[ip.b][1] + force[1],forces[ip.b][2] = forces[ip.b][2] + force[2]
*		3.返回e
* 
*/
fl eval_interacting_pairs_deriv(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces) { // adds to forces  // clean up
	const fl cutoff_sqr = p.cutoff_sqr();
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		vec r; r = coords[ip.b] - coords[ip.a]; // a -> b
		fl r2 = sqr(r);
		if(r2 < cutoff_sqr) {
			pr tmp = p.eval_deriv(ip.type_pair_index, r2);
			vec force; force = tmp.second * r;
			curl(tmp.first, force, v);
			e += tmp.first;
			// FIXME inefficient, if using hard curl
			forces[ip.a] -= force; // we could omit forces on inflex here
			forces[ip.b] += force;
		}
	}
	return e;
}

/*
* model的成员函数evali
* input:precalculate类常量引用p，vec类常量引用v
* output:double类型数
* 作用：1.定义一个double类型e=0
*		2.遍历ligands的元素				
*				将p,v[0],ligands第i个元素的成员pairs,internal_coords输入eval_interacting_pairs函数，返回值累加到e上
*		3.返回e
*/
fl model::evali(const precalculate& p,                                  const vec& v                          ) const { // clean up
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, internal_coords); // probably might was well use coords here
	return e;
}

/*
* model的成员函数evale
* input:precalculate类常量引用p，igrid类常量引用ig，vec类常量引用v
* output:double类型数
* 作用：1.定义一个double类型e,具体查看输入ig中的eval定义,ig应该是igrid派生类对象
*		2.遍历ligands的元素
*				将p,v[2],other_pairs, coords输入eval_interacting_pairs函数，返回值累加到e上
*		3.返回e
*/
fl model::evale(const precalculate& p, const igrid& ig, const vec& v                          ) const { // clean up
	fl e = ig.eval(*this, v[1]);
	e += eval_interacting_pairs(p, v[2], other_pairs, coords);
	return e;
}


/*
* model的成员函数eval
* input:precalculate类常量引用p，igrid类常量引用ig，vec类常量引用v，conf类常量引用c
* output:double类型数
* 作用：1.根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
*		2.将p,ig,v输入evale函数，返回值赋值到e上
*		3.遍历ligands的元素
*				将p,v[0],ligands第i个元素的成员pairs,coords输入eval_interacting_pairs函数，返回值累加到e上
*		4.返回e
*/
fl model::eval         (const precalculate& p, const igrid& ig, const vec& v, const conf& c           ) { // clean up
	set(c);
	fl e = evale(p, ig, v);
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords
	return e;
}

/*
* model的成员函数eval_deriv
* input:precalculate类常量引用p，igrid类常量引用ig，vec类常量引用v，conf类常量引用c，change类引用g
* output:double类型数
* 作用：1.根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
*		2.定义一个double类型e,具体查看输入ig中的eval_deriv定义,ig应该是igrid派生类对象
*		3.将p,v[2], other_pairs, coords,minus_forces输入eval_interacting_pairs_deriv函数，返回值累加到e上
*		4.遍历ligands的元素
*				将p,v[0],ligands第i个元素的成员pairs,coord,minus_forces输入eval_interacting_pairs_deriv函数，返回值累加到e上
*		5.将coords, minus_forces, g.ligands输入ligands.derivative中，分别对ligands的所有元素求导数，ligands[i].derivative(coords, minus_forces, g.ligands[i])
*		6.将coords, minus_forces, g.flex输入flex.derivative中，分别对flex的所有元素求导数，flex[i].derivative(coords, minus_forces, g.flex[i])
*		7.返回e
*/
fl model::eval_deriv  (const precalculate& p, const igrid& ig, const vec& v, const conf& c, change& g) { // clean up
	set(c);
	fl e = ig.eval_deriv(*this, v[1]); // sets minus_forces, except inflex
	e += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords, minus_forces); // adds to minus_forces
	VINA_FOR_IN(i, ligands)
		e += eval_interacting_pairs_deriv(p, v[0], ligands[i].pairs, coords, minus_forces); // adds to minus_forces
	// calculate derivatives
	ligands.derivative(coords, minus_forces, g.ligands);
	flex   .derivative(coords, minus_forces, g.flex); // inflex forces are ignored
	return e;
}




/*
* model的成员函数eval_intramolecular
* input:precalculate类常量引用p，vec类常量引用v，conf类常量引用c
* output:double类型数
* 作用：1.根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
*		2.定义一个double类型e = 0
*		3.遍历ligands的元素
*				将p,v[0],ligands第i个元素的成员pairs,coords输入eval_interacting_pairs_deriv函数，返回值累加到e上
*		4.将m_atom_typing_used输入num_atom_types得到m_atom_typing_used对应的TYPE_SIZE，赋值给nat
*		  返回p.m_cutoff_sqr赋值给cutoff_sqr
*		5.i从0循环到m_num_movable_atoms-1
*				判断i是否在ligands容器任意一个元素的成员begin和end范围之内，是则直接下一个i循环，否则继续执行以下语句		
*				将atmos的第i个元素绑定在a上
*				返回m_atom_typing_used对应的a的unit变量（el, ad, xs, sy之一）赋值给t1
*				当t1>=nat时直接下一个i循环，否则继续执行以下语句
*				j从0循环到grid_atoms.size()-1（遍历grid_atoms的元素）
*						将atmos的第j个元素绑定在b上
*						返回m_atom_typing_used对应的b的unit变量（el, ad, xs, sy之一）赋值给t2
*						当t2>=nat时直接下一个j循环，否则继续执行以下语句
*						r2 = (coords[i][0]-b.coords[0])^2 +(coords[i][1]-b.coords[1])^2 + (coords[i][2]-b.coords[2])^2
*						如果r2 < cutoff_sqr
*								如果t1<=t2，type_pair_index = t1 + t2*(t2+1)/2；如果t1>t2，type_pair_index = t2 + t1*(t1+1)/2
*								this_e = p.data.m_data[type_pair_index].fast[sz(factor * r2)]，factor为p.data.m_data[type_pair_index]的成员，sz为强制类型转换为unsigned int		
*								当this_e>0并且v[1]<0.1*max_fl(当前处理器下double型变量最大值)时，如果v[1] < epsilon_fl=2.22045e-16，this_e=0，否则this_e=this_e*(v[1] / (v[1] + this_e))
*								将this_e累加到e上
* 
*		6.i从0循环到other_pairs.size-1
*				将other_pairs的第interacting_pair个元素绑定在a上
*				如果pair.a在ligands容器任意一个元素的成员begin和end范围之内或者pair.b在ligands容器任意一个元素的成员begin和end范围之内，则直接下一个i循环，否则继续执行以下语句	
*				令r2 = （coords[pair.a][0] - coords[pair.b][0]）^2 + （coords[pair.a][1] - coords[pair.b][1]）^2 + （coords[pair.a][2] - coords[pair.b][2]）^2
*				如果r2 < cutoff_sqr
*						this_e = p.data.m_data[pair.type_pair_index].fast[sz(factor * r2)]，factor为p.data.m_data[type_pair_index]的成员，sz为强制类型转换为unsigned int		
*						当this_e>0并且v[2]<0.1*max_fl(当前处理器下double型变量最大值)时，如果v[2] < epsilon_fl=2.22045e-16，this_e=0，否则this_e=this_e*(v[2] / (v[2] + this_e))
*						将this_e累加到e上
*		7.返回e
*/
fl model::eval_intramolecular(const precalculate& p, const vec& v, const conf& c) {
	set(c);
	fl e = 0;

	// internal for each ligand
	VINA_FOR_IN(i, ligands)
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords

	sz nat = num_atom_types(atom_typing_used());
	const fl cutoff_sqr = p.cutoff_sqr();

	// flex-rigid
	VINA_FOR(i, num_movable_atoms()) {
		if(find_ligand(i) < ligands.size()) continue; // we only want flex-rigid interaction
		const atom& a = atoms[i];
		sz t1 = a.get(atom_typing_used());
		if(t1 >= nat) continue;
		VINA_FOR_IN(j, grid_atoms) {
			const atom& b = grid_atoms[j];
			sz t2 = b.get(atom_typing_used());
			if(t2 >= nat) continue;
			fl r2 = vec_distance_sqr(coords[i], b.coords);
			if(r2 < cutoff_sqr) {
				sz type_pair_index = triangular_matrix_index_permissive(nat, t1, t2);
				fl this_e = p.eval_fast(type_pair_index, r2);
				curl(this_e, v[1]);
				e += this_e;
			}
		}
	}

	// flex-flex
	VINA_FOR_IN(i, other_pairs) {
		const interacting_pair& pair = other_pairs[i];
		if(find_ligand(pair.a) < ligands.size() || find_ligand(pair.b) < ligands.size()) continue; // we only need flex-flex
		fl r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
		if(r2 < cutoff_sqr) {
			fl this_e = p.eval_fast(pair.type_pair_index, r2);
			curl(this_e, v[2]);
			e += this_e;
		}
	}
	return e;
}

/*
* 定义model的成员函数eval_adjusted 
* input:scoring_function类常量引用sf，precalculate类常量引用p，igrid类常量引用ig，vec类常量引用v，conf类常量引用c，double类型数intramolecular_energy（分子内能量）
* ouput:double类型数
* 作用：1.将p, ig, v, c输入eval函数中获得估计值e
*		2.将本model和e-intramolecular_energy输入sf成员函数conf_independent，返回独立配置评估值
*		
*/
fl model::eval_adjusted      (const scoring_function& sf, const precalculate& p, const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy) {
	fl e = eval(p, ig, v, c); // sets c
	return sf.conf_independent(*this, e - intramolecular_energy);
}

/*
* 定义model的成员函数rmsd_lower_bound_asymmetric 下界非对称均方根误差
* input:model类常量引用x,y
* ouput:double类型数
* 作用：1.获取x的成员m_num_movable_atoms赋值给n，判断n是否等于y的成员m_num_movable_atoms，不相等则报错，终止程序执行
*		2.定义一个double类型数sum = 0,unsigned类型数counter = 0
*		3.i从0循环到n-1
*			将x的成员atoms的第i个元素绑定到a上
*			如果a的成员el不等于0
*				（1）定义一个double类型数r2等于当前处理器下double型变量最大值
*				（2）j从0循环到n-1
*						将y的成员atoms的第j个元素绑定到b上
*						如果a的成员el与b的成员el相等并且b不是氢原子（!(b.ad==6||b.ad==12)）
*							[1]定义一个double类型数this_r2等于x成员coords第i个元素与y成员coords第j个元素逐个data数组元素差的平方
*							[2]如果this_r2 < r2，将this_r2赋值给r2
*				判断r2<0.1*max_fl(当前处理旗下double型变量最大值),如果不小于报错，终止程序执行
*				将r2累加到sum上，counter自增一
*		4.如果counter等于0，返回0，否则返回sum / counter的开方数	
*/

fl model::rmsd_lower_bound_asymmetric(const model& x, const model& y) const { // actually static
	sz n = x.m_num_movable_atoms; 
	VINA_CHECK(n == y.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, n) {
		const atom& a =   x.atoms[i];
		if(a.el != EL_TYPE_H) {
			fl r2 = max_fl;
			VINA_FOR(j, n) {
				const atom& b = y.atoms[j];
				if(a.same_element(b) && !b.is_hydrogen()) {
					fl this_r2 = vec_distance_sqr(x.coords[i], 
					                              y.coords[j]);
					if(this_r2 < r2)
						r2 = this_r2;
				}
			}
			assert(not_max(r2));
			sum += r2;
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

/*
* 定义model的成员函数rmsd_lower_bound  下限均方根误差
* input:model类常量引用m
* ouput:double类型数
* 作用：1.将*this（本model）和m按顺序输入函数rmsd_lower_bound_asymmetric，将m和*this（本model）按顺序输入函数rmsd_lower_bound_asymmetric，返回这两个返回值较大的
*/
fl model::rmsd_lower_bound(const model& m) const {
	return (std::max)(rmsd_lower_bound_asymmetric(*this,     m),
		            rmsd_lower_bound_asymmetric(    m, *this));
}

/*
* 定义model的成员函数rmsd_upper_bound  上限均方根误差
* input:model类常量引用m
* ouput:double类型数
* 作用：1.检查本model成员m_num_movable_atoms与m的成员m_num_movable_atoms是否相等，不相等则报错，终止程序执行
*		2.定义一个double类型数sum = 0,unsigned类型数counter = 0
*		3.i从0循环到（本model成员m_num_movable_atoms-1）
*			将本model的成员atoms的第i个元素绑定到a上
*			将m的成员atoms的第i个元素绑定到b上
*			检查a成员ad与b的成员ad是否相等，不相等则报错，终止程序执行
*			检查a成员xs与b的成员xs是否相等，不相等则报错，终止程序执行
*			如果a成员el不等于0
*				（1）将本model成员coords第i个元素与m成员coords第i个元素逐个data数组元素差的平方累加到sum上
*				（2）counter自增一
*		4.如果counter等于0，返回0，否则返回sum / counter的开方数
*/
fl model::rmsd_upper_bound(const model& m) const {
	VINA_CHECK(m_num_movable_atoms == m.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, m_num_movable_atoms) {
		const atom& a =   atoms[i];
		const atom& b = m.atoms[i];
		assert(a.ad == b.ad);
		assert(a.xs == b.xs);
		if(a.el != EL_TYPE_H) {
			sum += vec_distance_sqr(coords[i], m.coords[i]);
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

/*
* 定义model的成员函数rmsd_ligands_upper_bound 配体上限均方根误差
* input:model类常量引用m
* ouput:double类型数
* 作用：1.检查本model成员ligands容器的大小与m的成员ligands容器的大小是否相等，不相等则报错，终止程序执行
*		2.定义一个double类型数sum = 0,unsigned类型数counter = 0
*		3.ligand_i从0循环到（本model成员ligands的大小-1）
*			将本model的成员ligands的第ligand_i个元素绑定到lig  上
*			将  m    的成员ligands的第ligand_i个元素绑定到m_lig上
*			检查lig成员begin与m_lig的成员begin 是否相等，不相等则报错，终止程序执行
*			检查lig成员end  与m_lig的成员end   是否相等，不相等则报错，终止程序执行
*			i从lig成员begin循环到（lig成员end-1）
*				 将本model的成员atmos的第i个元素绑定到a上
* 				 将  m    的成员atmos的第i个元素绑定到b上
*				 检查a成员ad与b的成员ad是否相等，不相等则报错，终止程序执行
* 				 检查a成员xs与b的成员xs是否相等，不相等则报错，终止程序执行
*				 如果a成员el不等于0
*				 	（1）将本model成员coords第i个元素与m成员coords第i个元素逐个data数组元素差的平方累加到sum上
* 				 	（2）counter自增一
*		4.如果counter等于0，返回0，否则返回sum / counter的开方数
*/
fl model::rmsd_ligands_upper_bound(const model& m) const {
	VINA_CHECK(ligands.size() == m.ligands.size());
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR_IN(ligand_i, ligands) {
		const ligand&   lig =   ligands[ligand_i];
		const ligand& m_lig = m.ligands[ligand_i];
		VINA_CHECK(lig.begin == m_lig.begin);
		VINA_CHECK(lig.end   == m_lig.end);
		VINA_RANGE(i, lig.begin, lig.end) {
			const atom& a =   atoms[i];
			const atom& b = m.atoms[i];
			assert(a.ad == b.ad);
			assert(a.xs == b.xs);
			if(a.el != EL_TYPE_H) {
				sum += vec_distance_sqr(coords[i], m.coords[i]);
				++counter;
			}
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

/*
* 定义model的成员函数verify_bond_lengths 验证连接长度
* input:无
* ouput:无
* 作用：i从0循环到（model成员grid_atoms的大小和atoms的大小之和-1）
*				将i输入到model成员函数sz_to_atom_index，返回的值赋值到ai上，当0<=i<=(grid_atoms的大小-1)时，ai为atom_index(i,true)，否则为atom_index(i - grid_atoms.size(), false)
* 				将ai输入model成员函数get_atom，ai为atom_index(i,true)时，获得model成员grid_atoms的第i（ai的成员i）个元素，ai为atom_index(i - grid_atoms.size(), false)时，获得model成员atoms的第i（ai的成员i）个元素，返回的值绑定到a上
* 				j从0循环到a成员bonds的大小-1
*					a成员bonds的第j个元素绑定到b上
*					将ai与b的成员connected_atom_index输入距离平方函数，得到返回值并计算开方数，得到d（见公式3）
* 					比较d与b的成员length，相等则定义ok为真，否则为假
*					如果ok为假
*						在控制台显示d=d的数值
*						在控制台显示b.length=b.length的数值
*					检查ok是否为真，不为真则报错，终止程序执行
*						
*		
*/
void model::verify_bond_lengths() const {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		const atom& a = get_atom(ai);
		VINA_FOR_IN(j, a.bonds) {
			const bond& b = a.bonds[j];
			fl d = std::sqrt(distance_sqr_between(ai, b.connected_atom_index));
			bool ok = eq(d, b.length);
			if(!ok) {
				VINA_SHOW(d);
				VINA_SHOW(b.length);
			}
			VINA_CHECK(ok);
		}
	}
}

/*
* 定义model的成员函数check_internal_pairs 验证连接长度
* input:无
* ouput:无
* 作用：遍历model成员容器ligands的元素
*			将ligands的第i个元素绑定到lig上
* 			将lig的成员pairs绑定到pairs上
*			遍历容器pairs的元素
*				将pairs的第j个元素绑定到ip上
*				检查ip的成员a是否大于等于lig的成员begin，小于则报错，终止程序执行
* 				检查ip的成员b是否小于lig的成员end，不小于则报错，终止程序执行
*				
*
* 检查配体的begin是否小于等于相互作用对的a，检查配体的end是否大于相互作用对的b
* 
*/
void model::check_internal_pairs() const {
	VINA_FOR_IN(i, ligands) {
		const ligand& lig = ligands[i];
		const interacting_pairs& pairs = lig.pairs;
		VINA_FOR_IN(j, pairs) {
			const interacting_pair& ip = pairs[j];
			VINA_CHECK(ip.a >= lig.begin);
			VINA_CHECK(ip.b  < lig.end);
		}
	}
}

/*
* model的成员函数:about
* input:无
* output:无
* 作用：在控制台显示atom_typing_used()   = model成员m_atom_typing_used的值
* 		在控制台显示num_movable_atoms()  = model成员m_num_movable_atoms的值
*		在控制台显示num_internal_pairs() = model成员ligands所有元素的成员pairs的大小之和
* 		在控制台显示num_other_pairs()	 = model成员other_pairs的大小
*		在控制台显示num_ligands()		 = model成员ligands的大小
* 		在控制台显示num_flex()		     = model成员flex的大小
*
*/
void model::about() const {
	VINA_SHOW(atom_typing_used());
	VINA_SHOW(num_movable_atoms());
	VINA_SHOW(num_internal_pairs());
	VINA_SHOW(num_other_pairs());
	VINA_SHOW(num_ligands());
	VINA_SHOW(num_flex());
}

/*
* model的成员函数:about
* input:无
* output:无
* 作用：
*		遍历coords的元素，在控制台打印a的成员，即
* 		打印
* 			coords:
* 			[coords[i][0],coords[i][1],coords[i][2]]
*			........................................
*		遍历internal_coords的元素，在控制台打印a的成员，即
* 		打印
* 			internal_coords:
* 			[internal_coords[i][0],internal_coords[i][1],internal_coords[i][2]]
* 			...................................................................
*		遍历atoms的元素，分别绑定在a上，在控制台打印a的成员，即
* 		打印
* 			atoms:
* 			a.el a.ad a.xs a.sy a.charge
* 			a.bonds的大小 [a.coords[0],a.coords[1],a.coords[2]]       
*			...................................................
*		遍历grid_atoms的元素，分别绑定在a上，在控制台打印a的成员，即
*		打印
*			grid_atoms:
*			a.el a.ad a.xs a.sy a.charge
*			a.bonds的大小 [a.coords[0],a.coords[1],a.coords[2]]        
* 		    ...................................................
*/
void model::print_stuff() const {
	std::cout << "coords:\n";
	VINA_FOR_IN(i, coords)
		printnl(coords[i]);

	std::cout << "internal_coords:\n";
	VINA_FOR_IN(i, internal_coords)
		printnl(internal_coords[i]);

	std::cout << "atoms:\n";
	VINA_FOR_IN(i, atoms) {
		const atom& a = atoms[i];
		std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
		std::cout << a.bonds.size() << "  "; printnl(a.coords);
	}

	std::cout << "grid_atoms:\n";
	VINA_FOR_IN(i, grid_atoms) {
		const atom& a = grid_atoms[i];
		std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
		std::cout << a.bonds.size() << "  "; printnl(a.coords);
	}
	about();
}


/*
*input:double型r和covalent_r
*output：若 r < 0 或者 covalent < 运行编译程序的计算机所能识别的最小非零浮点数 时报错
*否则r / covalent_r >2 时返回0
*    r / covalent_r <= 2 时返回1-（r / covalent_r）*（r / covalent_r）/4
*/
fl pairwise_clash_penalty(fl r, fl covalent_r) {
	// r = 0          -> max_penalty 
	// r = covalent_r -> 1
	// elsewhere      -> hyperbolic function
	assert(r >= 0);//r<0报错
	assert(covalent_r > epsilon_fl);//covalent<运行编译程序的计算机所能识别的最小非零浮点数时报错
	const fl x = r / covalent_r;
	if(x > 2) return 0;
	return 1-x*x/4;
}

/*
* model成员函数clash_penalty_aux，冲突惩罚
*input:interacting_pair类vector
*output:doule类数
*作用：
*	1.定义一个double变量e = 0
*	2.遍历vector
*		(1)将vector的第i个元素绑定到ip上
*       (2)计算公式2，得到r
*       (3)计算atoms第ip.a（ip的成员a）和ip.b（ip的成员b）的共价半径和，得到covalent_r
*       (4)计算成对冲突惩罚，将r和covalent_r代入函数pairwise_clash_penalty，并将其累加到e上
*   3.返回e
*/
fl model::clash_penalty_aux(const interacting_pairs& pairs) const {
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		const fl r = std::sqrt(vec_distance_sqr(coords[ip.a], coords[ip.b]));
		const fl covalent_r = atoms[ip.a].covalent_radius() + atoms[ip.b].covalent_radius();
		e += pairwise_clash_penalty(r, covalent_r);
	}
	return e;
}

/*
* model成员函数clash_penalty，冲突惩罚
*input:无
*output:doule类数
*作用：
*	1.定义一个double变量e = 0
*	2.遍历model成员ligands，
*         将ligands第i个元素的成员pairs输入clash_penalty_aux，得到结果累加到e上
*   3.将model成员other_pairs输入clash_penalty_aux，得到结果累加到e上
*   4.返回e
*/
fl model::clash_penalty() const {
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += clash_penalty_aux(ligands[i].pairs);
	e += clash_penalty_aux(other_pairs);
	return e;
}
