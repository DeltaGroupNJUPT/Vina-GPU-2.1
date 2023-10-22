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

#include "terms.h"
#include "brick.h"

/*		outout——tmp:double 向量；
		定义了一个double向量，然后将结构体conf_independent_inputs中的
		num_tors、num_rotors、num_heavy_atoms、num_hydrophobic_atoms、ligand_max_num_h_bonds、num_ligands、ligand_lengths_sum
		按顺序添加到tmp的向量的最末尾一个元素后面。
*/

conf_independent_inputs::operator flv() const {
	flv tmp;                                     //double向量
	tmp.push_back(num_tors);                    
	tmp.push_back(num_rotors);
	tmp.push_back(num_heavy_atoms);
	tmp.push_back(num_hydrophobic_atoms);
	tmp.push_back(ligand_max_num_h_bonds);
	tmp.push_back(num_ligands);
	tmp.push_back(ligand_lengths_sum);
	return tmp;
}

/*	
	input——m：model结构体，i：atom_index；
	output——acc：unsigned
	bonds：bond结构体向量，等于atom结构体中的 bonds；
	b：bond结构体；
	a：atom 结构体；

	a是以b输出原子索引获得grid_atoms向量中的内容
	若a不是氢则acc加1
*/
unsigned conf_independent_inputs::num_bonded_heavy_atoms(const model& m, const atom_index& i) const { // FIXME? - could be static, but I don't feel like declaring function friends
	unsigned acc = 0;
	const std::vector<bond>& bonds = m.get_atom(i).bonds;
	VINA_FOR_IN(j, bonds) {                   //for(j=0;j<bonds.size;j++)
		const bond& b = bonds[j];
		const atom& a = m.get_atom(b.connected_atom_index);
		if(!a.is_hydrogen())    //是否为氢(ad==6||ad==12);
			++acc;
	}
	return acc;
}

/*  重配体原子的可旋转键数
	input——m：model结构体；i：atom_index结构体
	output——acc：unsigned
	bonds：bond结构体向量，等于atom结构体中的 bonds；
	b：bond结构体；
	a：atom 结构体；

	a是以b输出原子索引获得grid_atoms向量中的内容
	若a不是氢、b.rotatable为true并且num_bonded_heavy_atoms结构体中acc》1
	则acc加1
*/
// FIXME? - could be static, but I don't feel like declaring function friends
unsigned conf_independent_inputs::atom_rotors(const model& m, const atom_index& i) const { // the number of rotatable bonds to heavy ligand atoms
	unsigned acc = 0;
	const std::vector<bond>& bonds = m.get_atom(i).bonds;
	VINA_FOR_IN(j, bonds) {                          //for(i=0;i<bonds.size;i++)
		const bond& b = bonds[j];
		const atom& a = m.get_atom(b.connected_atom_index);
		if(b.rotatable && !a.is_hydrogen() && num_bonded_heavy_atoms(m, b.connected_atom_index) > 1) // not counting CH_3, etc
			++acc;
	}
	return acc;
}


/*  重配体原子的可旋转键数
input——m：model结构体；
	num_tors、num_rotors、num_heavy_atoms、num_hydrophobic_atoms、
	ligand_max_num_h_bonds、ligand_lengths_sum这些double类型数据都初始化为0
	num_ligands为ligands的元素个数（vector_mutable继承vector）
*/
conf_independent_inputs::conf_independent_inputs(const model& m) {
	num_tors = 0;
	num_rotors = 0;
	num_heavy_atoms = 0;
	num_hydrophobic_atoms = 0;
	ligand_max_num_h_bonds = 0;
	num_ligands = m.ligands.size();
	ligand_lengths_sum = 0;

	VINA_FOR_IN(ligand_i, m.ligands) {                     //for(ligand_i=0;ligand_i< m.ligands.size;ligand_i++)
		const ligand& lig = m.ligands[ligand_i];
		ligand_lengths_sum += m.ligand_length(ligand_i);
		VINA_RANGE(i, lig.begin, lig.end) {                //for(i=lig.begin;i<lig.end;i++)
			const atom& a = m.atoms[i];
			if(a.el != EL_TYPE_H) {
				unsigned ar = atom_rotors(m, atom_index(i, false));

				num_tors += 0.5 * ar;

				if(ar > 2) num_rotors += 0.5;
				else num_rotors += 0.5 * ar;

				++num_heavy_atoms;
				if(xs_is_hydrophobic(a.xs))
					++num_hydrophobic_atoms;

				if(xs_is_acceptor(a.xs) || xs_is_donor(a.xs))
					++ligand_max_num_h_bonds;
			}
		}
	}
}


/*		
		outout——tmp: string 向量；
		定义了一个string向量，然后将结构体conf_independent_inputs中的
		“num_tors”、“num_rotors”、“num_heavy_atoms”、“num_hydrophobic_atoms”、
		“ligand_max_num_h_bonds”、“num_ligands”、“ligand_lengths_sum”字符串
		按顺序添加到tmp的向量的最末尾一个元素后面。
		最后转换了tmp使其成为double向量并判断交换后的tmp double向量元素个数是否与转化前的tmp string向量元素个数相等，
		是的话返回tmp（转化前的）
*/

std::vector<std::string> conf_independent_inputs::get_names() const { // FIXME should probably be static
	std::vector<std::string> tmp;
	tmp.push_back("num_tors");
	tmp.push_back("num_rotors");
	tmp.push_back("num_heavy_atoms");
	tmp.push_back("num_hydrophobic_atoms");
	tmp.push_back("ligand_max_num_h_bonds");
	tmp.push_back("num_ligands");
	tmp.push_back("ligand_lengths_sum");
	VINA_CHECK(static_cast<flv>(*this).size() == tmp.size()); // FIXME?  转换之后操作的并非原有对象，
															   //而是一个副本，是一个临时变量。
	return tmp;
}

/*
num_tors、num_rotors、num_heavy_atoms、num_hydrophobic_atoms、
ligand_max_num_h_bonds、num_ligands、ligand_lengths_sum这些double类型数据都初始化为0
*/
conf_independent_inputs::conf_independent_inputs() : 
	num_tors(0), num_rotors(0), num_heavy_atoms(0), 
	num_hydrophobic_atoms(0), ligand_max_num_h_bonds(0), num_ligands(0), 
	ligand_lengths_sum(0) {}

/*
	input——a、b：double向量；
	inter——n：unint
	output——acc：double；
	n为a、b向量元素个数的最小值

	demo：a{2，3，4，5}，a{1，1，2}，则n=3
		  acc=2*1+3*1+2*4=13

*/
inline fl inner_product_shortest(const flv& a, const flv& b) {
	sz n = (std::min)(a.size(), b.size());
	fl acc = 0;
	VINA_FOR(i, n)                             //for(i=0;i<n;i++)
		acc += a[i] * b[i];
	return acc;
}


/*
       input——weights：double向量，include_internal：bool；
	   inter——e为factors结构体中的e（double向量）；i为factors结构体中的i（double向量）；
	   output—— tmp：double

	   demo：
	   e{2，3，4，5}，weights{1，1，2}，i{2，3}
	   if（include_internal=false）（e*weights）
			 n=3
			 tmp=2*1+3*1+2*4=13；
	   else     （i*weights）
	         n=2
	         tmp=2*1+3*1=5

*/
fl factors::eval(const flv& weights, bool include_internal) const {
	fl tmp = inner_product_shortest(e, weights);     //e为factors结构体中的e（double向量）
	if(include_internal)
		tmp += inner_product_shortest(i, weights);   //i为factors结构体中的i（double向量）
	return tmp;
}


// terms
/*
 input——enabled_only：bool；
 inter——tmp：字符串向量；
	该函数通过输入的bool类型enabled_only和term_set结构体中的enabled的bool向量来判断是否在tmp向量末尾添加term_set结构体中fun[i].name
		（enabled_only为0或者enabled[i]大于等于1）
*/
std::vector<std::string> terms::get_names(bool enabled_only) const { // does not include conf-independent
	std::vector<std::string> tmp;
	/*先确定enabled向量的元素个数是否和fun向量元素相同，是的话继续否则终止；
		最后如果enabled_only不为0或者enabled[i]大于等于1则在tmp向量末尾添加fun[i].name（term结构体中name）*/
	distance_additive_terms.get_names(enabled_only, tmp);//一个继承了distance_additive_terms结构体的term_set结构体
	           usable_terms.get_names(enabled_only, tmp);//一个继承了usable结构体的term_set结构体，在上面处理完的tmp字符串向量末尾添加
	         additive_terms.get_names(enabled_only, tmp);//一个继承了additive结构体的term_set结构体，在上面处理完的tmp字符串向量末尾添加
	   intermolecular_terms.get_names(enabled_only, tmp);//一个继承了intermolecular结构体的term_set结构体，在上面处理完的tmp字符串向量末尾添加
	return tmp;
}


/*
	该函数返回distance_additive_terms结构体中fun向量元素个数、usable_terms结构体中fun向量元素个数、
	additive_terms中fun向量元素个数的总和

*/
sz terms::size_internal() const {
	return distance_additive_terms.size() + usable_terms.size() + additive_terms.size();
	//size()：该函数返回fun向量元素个数，term_set结构体中fun（）任意类型的指针向量，这个和vector类似不过效率更高
	//
}



/*
	input——enabled_only：bool类型
	inter——enabled：bool类型向量
	在enabled_only=false或者term_set结构体中的enabled[i]>0情况下对fun[i]的元素个数进行累加（i=0-fun向量元素个数）
*/
sz terms::size_conf_independent(bool enabled_only) const { // number of parameters does not necessarily equal the number of operators
	sz acc = 0;
	VINA_FOR_IN(i, conf_independent_terms)    //for(i=0;i<conf_independent_terms.size(fun向量元素个数);i++)
		if(!enabled_only || conf_independent_terms.enabled[i])//enabled_only为0或者enabled[i]大于等于1
			acc += conf_independent_terms[i].size();         //conf_independent_terms[i]=fun[i]
	return acc;
}


/*

	这个函数是在寻找 distance_additive_terms.fun[i].cutoff、usable_terms.fun[i].cutoff、
	additive_terms.fun[i].cutoff中的最大值

*/

fl terms::max_r_cutoff() const {
	fl tmp = 0;
	tmp = (std::max)(tmp, distance_additive_terms.max_cutoff());//max_cutoff():寻找fun[i].cutoff中的最大值
	tmp = (std::max)(tmp,            usable_terms.max_cutoff());
	tmp = (std::max)(tmp,          additive_terms.max_cutoff());
	return tmp;
}

/*
	input——m：model结构体，i、j：atom_index结构体，r：double，out：double向量
	inter——a、b：atom结构体
	a、b：判断i、j的成员in_grid是否为真，是则返回grid_atoms容器的第（对象i、j的成员i、j的数值）元素常量引用，否则返回atoms容器的第（对象i、j的成员i、j的数值）元素常量引用
*/
void terms::eval_additive_aux(const model& m, const atom_index& i, const atom_index& j, fl r, flv& out) const { // out is added to
	const atom& a = m.get_atom(i);
	const atom& b = m.get_atom(j);

	sz offset = 0;
	VINA_FOR_IN(k, distance_additive_terms)                      //for(k=0;k<term_set.fun.size().size;k++)
		if(r < distance_additive_terms[k].cutoff)                //(r<term_set<distance_additive_terms>.fun[i].cutoff)
			out[k] += distance_additive_terms[k].eval(a, b, r);  //get（）函数是结构体atom_base成员，
															     //因为atom_typing_use=XS，所以返回XS_TYPE_SIZE

	offset += distance_additive_terms.size();                    
	VINA_FOR_IN(k, usable_terms)								 //for(k=0;k<term_set<usable_terms>.fun.size().size;k++)
		if(r < usable_terms[k].cutoff)
			out[offset + k] += usable_terms[k].eval(a, b, r);    //在上面的out[k]后面添加

	offset += usable_terms.size();
	VINA_FOR_IN(k, additive_terms)							      //for(k=0;k<term_set<additive_terms>.fun.size().size;k++)                            
		if(r < additive_terms[k].cutoff)
			out[offset + k] += additive_terms[k].eval(m, i, j);   //在上面的out[k]后面添加
	
	VINA_CHECK(offset + additive_terms.size() == size_internal());//size_internal（）该函数返回distance_additive_terms结构体中fun向量元素个数、usable_terms结构体中fun向量元素个数、
																//additive_terms中fun向量元素个数的总和
}


/*外部因素（可能吧）
	input——m：model结构体
	output——tmp：double向量

*/
flv terms::evale(const model& m) const {
	VINA_CHECK(m.ligands.size() == 1);      //ligands为模板参数为ligand结构体的vector_mutable结构体
	VINA_CHECK(m.flex.size() == 0);
	VINA_CHECK(m.atoms.size() == m.num_movable_atoms()); // no inflex

	flv tmp(size(), 0);            //size()该函数返回distance_additive_terms结构体中fun向量元素个数、
								   //usable_terms结构体中fun向量元素个数、
								   //additive_terms中fun向量元素个数、
								   //intermolecular_terms的向量元素个数总和
	fl max_r_cutoff_sqr = sqr(max_r_cutoff()); //在寻找 distance_additive_terms.fun[i].cutoff、usable_terms.fun[i].cutoff、
											   //additive_terms.fun[i].cutoff中的最大值的平方

	grid_dims box = m.movable_atoms_box(0); // add nothing，有3个元素的grid_dim结构体数组
	vec box_begin = grid_dims_begin(box);   //将三个grid_dim结构体中的begin存入向量
	vec box_end   = grid_dims_end  (box);   //将三个grid_dim结构体中的end存入向量

	szv relevant_atoms;                     //unint向量
	VINA_FOR_IN(j, m.grid_atoms)           //for(j=0;j< m.grid_atoms.size;j++)
		if(brick_distance_sqr(box_begin, box_end, m.grid_atoms[j].coords) < max_r_cutoff_sqr)
			/*
			demo：
			box_begin（1，2，3），	box_end（7，8，9），m.grid_atoms[j].coords（0，9，7）
			所以
			closest(1，8，7)（closest是m.grid_atoms[j].coords最接近box_begin、box_end中的值）
			返回
			(1 - 0) ^ 2 + (8 - 9) ^ 2 + (7 - 7) ^ 2 = 2
			max_r_cutoff_sq为distance_additive_terms.fun[i].cutoff、usable_terms.fun[i].cutoff、
			additive_terms.fun[i].cutoff中的最大值的平方
			*/
			relevant_atoms.push_back(j);    //如果if成立，将j加到relevant_atoms向量的最后面

	VINA_FOR(i, m.num_movable_atoms()) {                 //for(i=0;i<m.num_movable_atoms;i++)
		const vec& coords = m.coords[i];                 
		VINA_FOR_IN(relevant_j, relevant_atoms) {        //for(relevant_j=0;relevant_j<relevant_atoms.size;relevant_j++)
			const sz j = relevant_atoms[relevant_j];     //将relevant_atoms中的元素别名为j
			const atom& b = m.grid_atoms[j];             //将grid_atoms中的元素别名为b
			fl d2 = vec_distance_sqr(coords, b.coords); //计算两个vec中逐个data数组元素差的平方

			if(d2 > max_r_cutoff_sqr) continue; // most likely scenario   if成立直接下一次VINA_FOR_IN
			fl d = std::sqrt(d2);               //对d2开根号
			eval_additive_aux(m, atom_index(i, false), atom_index(j, true), d, tmp);
		}
	}
	sz offset = size_internal();       //该函数返回distance_additive_terms结构体中fun向量元素个数、usable_terms结构体中fun向量元素个数、
										//additive_terms中fun向量元素个数的总和
	VINA_FOR_IN(k, intermolecular_terms)
		tmp[offset + k] += intermolecular_terms[k].eval(m);

	VINA_CHECK(offset + intermolecular_terms.size() == tmp.size()); 
	return tmp;
}


/*内部因素（可能吧）
	input——m：model结构体
	output——tmp：double向量
*/

flv terms::evali(const model& m) const {
	VINA_CHECK(m.ligands.size() == 1);
	VINA_CHECK(m.flex.size() == 0);
	VINA_CHECK(m.atoms.size() == m.num_movable_atoms());

	flv tmp(size_internal(), 0);                         //size_internal（）该函数返回distance_additive_terms结构体中fun向量元素个数、usable_terms结构体中fun向量元素个数、
	                                                     //	additive_terms中fun向量元素个数的总和
	fl max_r_cutoff_sqr = sqr(max_r_cutoff());          //在寻找 distance_additive_terms.fun[i].cutoff、usable_terms.fun[i].cutoff、
											            //additive_terms.fun[i].cutoff中的最大值的平方
	const interacting_pairs& pairs = m.ligands.front().pairs;    //别名interacting_pair结构体向量
	VINA_FOR_IN(i, pairs) {                                      //for(i=0;i<pairs.size;i++)
		const interacting_pair& ip = pairs[i];                   //别名interacting_pair结构体
		fl d2 = vec_distance_sqr(m.coords[ip.a], m.coords[ip.b]);//计算两个vec中逐个data数组元素差的平方
		if(d2 > max_r_cutoff_sqr) continue; // most likely scenario    if成立直接下一个循环
		fl d = std::sqrt(d2);                                    //对d2开根号
		eval_additive_aux(m, atom_index(ip.a, false), atom_index(ip.b, false), d, tmp);
	}
	return tmp;
}


/*外部因素增强，只有单个配体系统需要用到这个程序
	input——m：model结构体
	output——tmp：double向量
*/
flv terms::evale_robust(const model& m) const {
	VINA_CHECK(m.ligands.size() == 1); // only single-ligand systems are supported by this procedure

	flv tmp(size(), 0);                             

	fl max_r_cutoff_sqr = sqr(max_r_cutoff());

	grid_dims box = m.movable_atoms_box(0);            //add nothing，有3个元素的grid_dim结构体
	vec box_begin = grid_dims_begin(box);	           //将三个grid_dim结构体中的begin存入向量
	vec box_end = grid_dims_end(box);                  //将三个grid_dim结构体中的end存入向量

	const sz n  = num_atom_types(m.atom_typing_used());/*n=
														原子数类型
														EL型返回11;
														AD型返回20;
														XS型返回17;
														SY型返回18
														*/

	std::vector<atom_index> relevant_atoms;              //atom_index结构体向量

	VINA_FOR_IN(j, m.grid_atoms) {
		const atom& a = m.grid_atoms[j];
		const sz t = a.get(m.atom_typing_used());        //枚举型t对应返回值
		if(brick_distance_sqr(box_begin, box_end, a.coords) < max_r_cutoff_sqr && t < n) // exclude, say, Hydrogens
	          /*
	          demo：
	          box_begin（1，2，3），	box_end（7，8，9），m.grid_atoms[j].coords（0，9，7）
	          所以
	          closest(1，8，7)（closest是m.grid_atoms[j].coords最接近box_begin、box_end中的值）
	          返回
	          (1 - 0) ^ 2 + (8 - 9) ^ 2 + (7 - 7) ^ 2 = 2
	          max_r_cutoff_sq为distance_additive_terms.fun[i].cutoff、usable_terms.fun[i].cutoff、
	          additive_terms.fun[i].cutoff中的最大值的平方
	          */
			relevant_atoms.push_back(atom_index(j, true));//将(atom_index(j, true)存入relevant_atoms结构体向量的最后面
	}

	VINA_FOR_IN(j, m.atoms) {
		const atom& a = m.atoms[j];
		const vec& a_coords = m.coords[j];
		if(m.find_ligand(j) < m.ligands.size()) continue; // skip ligand atoms, add only flex/inflex  成立直接进入下一个循环
		const sz t = a.get(m.atom_typing_used());
		if(brick_distance_sqr(box_begin, box_end, a_coords) < max_r_cutoff_sqr && t < n) // exclude, say, Hydrogens
			relevant_atoms.push_back(atom_index(j, false));   //将(atom_index(j, false)存入上面的relevant_atoms结构体向量的最后面
	}

	VINA_FOR_IN(lig_i, m.ligands) {											
		const ligand& lig = m.ligands[lig_i];								
		VINA_RANGE(i, lig.begin, lig.end) {									
			const vec& coords = m.coords[i];								//将m结构体向量中的coords向量别名
			const atom& a = m.atoms[i];										//将m结构体向量中的atoms结构体别名
			const sz t = a.get(m.atom_typing_used());						//将m结构体向量中的m_atom_typing_used变量别名

			if(t < n) { // exclude, say, Hydrogens
				VINA_FOR_IN(relevant_j, relevant_atoms) {
					const atom_index& j = relevant_atoms[relevant_j];
					fl d2 = vec_distance_sqr(coords, m.atom_coords(j));    //计算两个vec中逐个data数组元素差的平方
					if(d2 > max_r_cutoff_sqr) continue; // most likely scenario
					fl d = std::sqrt(d2);
					eval_additive_aux(m, atom_index(i, false), j, d, tmp);
				}
			}
		}
	}

	sz offset = size_internal();
	VINA_CHECK(intermolecular_terms.size() == 0);
	VINA_CHECK(offset + intermolecular_terms.size() == tmp.size());

	return tmp;
}


/*	获取外部内部元素
	input——m：model结构体
	output——tmp：factors结构体
*/
factors terms::eval(const model& m) const {
	factors tmp;
	tmp.e = evale(m);     //外部因素（可能吧）
	tmp.i = evali(m);     //内部因素（可能吧）
	return tmp;
}



/*
	input——in：conf_independent_inputs；it：迭代器：
	inter——conf_independent_terms：模板形参为conf_independent的term_set结构体
	output——x：double；
	if(conf_independent结构体中enabled[i]的值大于0)（i=0-conf_independent结构体中fun向量元素个数）
		x为  x + read_iterator(i) * in.******   // read_iterator(i)为迭代器指向下一个元素 ，******应该是取决于调用这个函数的前一个everthing.cpp中的函数决定的
*/
fl terms::eval_conf_independent(const conf_independent_inputs& in, fl x, flv::const_iterator& it) const { // evaluates enabled only
	VINA_FOR_IN(i, conf_independent_terms)   //for(i=0;i<conf_independent_terms.size(fun向量元素个数);i++ ）
		if(conf_independent_terms.enabled[i]) 
			x = conf_independent_terms[i].eval(in, x, it);     // x + read_iterator(i) * in.ligand_lengths_sum    
	return x;
}

/*
	input——v:double向量
	当distance_additive_terms、usable_terms、additive_terms、intermolecular_terms结构体中的enabled的值大于等于1时
	此时将迭代器里元素依次（即v中的元素）添加到tmp向量最末尾，否则不添加；
	
*/
flv terms::filter_external(const flv& v) const {
	flv tmp;
	flv::const_iterator i = v.begin();      //将迭代器i初始化为指向输入v向量的第一个元素
	distance_additive_terms.filter(i, tmp); //当term_set结构体中enabled的值大于等于1时将此时迭代器里元素添加到tmp向量最末尾
	           usable_terms.filter(i, tmp);
	         additive_terms.filter(i, tmp);
	   intermolecular_terms.filter(i, tmp);
	VINA_CHECK(i == v.end());              //判断迭代器i是否遍历输入v向量
	return tmp;
}


/*
         input——v:double向量
         当distance_additive_terms、usable_terms、additive_terms、intermolecular_terms结构体中的enabled的值大于等于1时
         此时将迭代器里元素依次（即v中的元素）添加到tmp向量最末尾，否则不添加
         同filter_external（）函数功能一样

*/
flv terms::filter_internal(const flv& v) const {
	flv tmp;
	flv::const_iterator i = v.begin();
	distance_additive_terms.filter(i, tmp);
	           usable_terms.filter(i, tmp);
	         additive_terms.filter(i, tmp);
	VINA_CHECK(i == v.end());
	return tmp;
}



/*		
		input——f：factors结构体；
		inter——tmp：factors结构体；
		tmp.e得到f.e过滤后的向量，tmp.i得到f.i过滤后的向量，最后返回factors结构体的函数。
*/
factors terms::filter(const factors& f) const {
	factors tmp;
	tmp.e = filter_external(f.e);
	tmp.i = filter_internal(f.i);
	return tmp;
}


/*
	这个函数用于打印term_set结构体中fun[i].name和继承了conf_independent_terms结构体
	的term_set结构体的中fun[i].name，以及fun向量元素个数总和和cutoff中的最大值
*/

void terms::display_info() const {
	std::vector<std::string> enabled_names = get_names(true);//返回包含term_set结构体中fun[i].name的string向量，
	std::cout << "Enabled terms: \n";  
	VINA_FOR_IN(i, enabled_names)                 //for(i=0;i<enabled_names.size;i++)
		std::cout << enabled_names[i] << '\n';   //将fun[i].name打印出来
	std::cout << '\n';

	std::vector<std::string> enabled_operators;
	conf_independent_terms.get_names(true, enabled_operators);   //out返回包含继承了conf_independent_terms结构体
																 //的term_set结构体的中fun[i].name的string向量，
	std::cout << "Enabled conf-independent operators: \n"; 
	VINA_FOR_IN(i, enabled_operators)						//for(i=0;i<enabled_operators.size;i++)
		std::cout << enabled_operators[i] << '\n';				//将fun[i].name打印出来
	std::cout << '\n';

	VINA_SHOW(size());                          //打印一次size（）=向量元素个数总和
												//该函数返回distance_additive_terms结构体中fun向量元素个数、
								                //usable_terms结构体中fun向量元素个数、
								                //additive_terms中fun向量元素个数、
									            //intermolecular_terms的fun向量元素个数总和

	VINA_SHOW(size_internal());               //打印一次size_internal()=fun向量元素个数的总和
											  //该函数返回distance_additive_terms结构体中fun向量元素个数、
										      //usable_terms结构体中fun向量元素个数、
										      //additive_terms中fun向量元素个数的总和

	VINA_SHOW(max_r_cutoff());              //打印一次max_r_cutoff()=cutoff的最大值
											//distance_additive_terms.fun[i].cutoff、usable_terms.fun[i].cutoff、
										    //additive_terms.fun[i].cutoff中的最大值
}
