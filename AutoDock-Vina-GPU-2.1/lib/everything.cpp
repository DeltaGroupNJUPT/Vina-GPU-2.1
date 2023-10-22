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

#include "everything.h"
#include "int_pow.h"



// distance_additive terms   输入double类型的x   width  输出e的（-（x/width）的平方） double
inline fl gaussian(fl x, fl width) {
	return std::exp(-sqr(x/width));
}


/*结构体electrostatic继承distance_additive    定义模板参数UINT类型i
* 定义double   cap
* 函数electrostatic给父类中的string参数name赋值
* 定义double类型函数eval  输入结构体a、b  double r  输出double      
	返回charge1*charge2 * cap或者r的i次方分之一和cap的较小值
*/
template<unsigned i>
struct electrostatic : public distance_additive {
	fl cap;
	electrostatic(fl cap_, fl cutoff_) : distance_additive(cutoff_), cap(cap_) {
		name = std::string("electrostatic(i=") + to_string(i) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom_base& a, const atom_base& b, fl r) const {
		fl tmp = int_pow<i>(r);   //r的i次方
		fl q1q2 = a.charge * b.charge;//charge1*charge2
		if(tmp < epsilon_fl) return q1q2 * cap;//epsilon_fl   2.2204460492503131e-016   
		else                 return q1q2 * (std::min)(cap, 1/int_pow<i>(r));//返回  r的i次方分之一和cap的较小值
	}
};        
/*
function：
输入：结构体  原子类型  输出  double类型数据
* 判断  如果结构体的参数ad<16 则返回数组atom_kind_data[ad]中的solvation  数组是结构体类型的
 如果结构体的参数xs<16则返回-0.00110

*/
fl solvation_parameter(const atom_type& a) {
	if(a.ad < AD_TYPE_SIZE) return ad_type_property(a.ad).solvation;//ad<20    返回数组atom_kind_data[ad]中的solvation  数组是结构体类型的
	else if(a.xs == XS_TYPE_Met_D) return metal_solvation_parameter;//metal_solvation_parameter = -0.00110
	VINA_CHECK(false); 
	return 0; // placating the compiler
}                
/*
function：
输入：结构体  原子类型  输出  double类型数据
* 判断  如果结构体的参数ad<16 则返回数组atom_kind_data[ad]中的volume  数组是结构体类型的
 如果结构体的参数xs<17则返回 （xs_vdw_radii[xs]）的三次方*4*pi / 3  

*/
fl volume(const atom_type& a) {
	if(a.ad < AD_TYPE_SIZE) return ad_type_property(a.ad).volume;
	else if(a.xs < XS_TYPE_SIZE) return 4*pi / 3 * int_pow<3>(xs_radius(a.xs));
	VINA_CHECK(false);
	return 0; // placating the compiler
}                

/*结构体ad4_solvation 继承distance_additive    
* 定义double    desolvation_sigma;solvation_q  bool类型  charge_dependent
* 函数ad4_solvation给父类中的string参数name赋值
* 定义double类型函数eval  输入结构体a、b  double r  
	返回double类型的数tmp
*/
struct ad4_solvation : public distance_additive {
	fl desolvation_sigma;
	fl solvation_q;
	bool charge_dependent;
	ad4_solvation(fl desolvation_sigma_, fl solvation_q_, bool charge_dependent_, fl cutoff_) : distance_additive(cutoff_), solvation_q(solvation_q_), charge_dependent(charge_dependent_), desolvation_sigma(desolvation_sigma_) {
		name = std::string("ad4_solvation(d-sigma=") + to_string(desolvation_sigma) + ", s/q=" + to_string(solvation_q) + ", q=" + to_string(charge_dependent) + ", c=" + to_string(cutoff) + ")";
	}//   name被赋值 ad4_solvation(d-sigma=    ,s/q=   ,q=     ,c=    )”
	fl eval(const atom_base& a, const atom_base& b, fl r) const {
		fl q1 = a.charge;
		fl q2 = b.charge;

		VINA_CHECK(not_max(q1));//q1<0.17976931348623158e+308
		VINA_CHECK(not_max(q2));//q2<0.17976931348623158e+308

		sz t1 = a.ad;
		sz t2 = b.ad;

		fl solv1 = solvation_parameter(a);
		fl solv2 = solvation_parameter(b);

		fl volume1 = volume(a);
		fl volume2 = volume(b);

		fl my_solv = charge_dependent ? solvation_q : 0;

		fl tmp = ((solv1 + my_solv * std::abs(q1)) * volume2 + 
			    (solv2 + my_solv * std::abs(q2)) * volume1) * std::exp(-sqr(r/(2*desolvation_sigma)));

		VINA_CHECK(not_max(tmp));//tmp<0.17976931348623158e+308
		return tmp;
	}
};
/*
定义double类型的函数    输入  uint xs_t1、xs_t2   输入double ：数组中的数据 xs_vdw_radii[xs_t1]+xs_vdw_radii[xs_t2]
*/
inline fl optimal_distance(sz xs_t1, sz xs_t2) {
	return xs_radius(xs_t1) + xs_radius(xs_t2);    
}
/*结构体gauss继承usable    
* 定义double   offset   width
* 函数gauss给父类中的string参数name赋值   name被赋值 gauss(o=     ,w=    ,c=     )
* 定义double类型函数eval  输入double类型t1   t2  r     输出double
	返回e的-{(r-xs_vdw_radii[xs_t1]-xs_vdw_radii[xs_t2]-offset)/width}的平方
*/
struct gauss : public usable {
	fl offset; // added to optimal distance
	fl width;
	gauss(fl offset_, fl width_, fl cutoff_) : usable(cutoff_), offset(offset_), width(width_) {
		name = std::string("gauss(o=") + to_string(offset) + ", w=" + to_string(width) + ", c=" + to_string(cutoff) + ")";
	}//name被赋值 gauss(o=, w=, c = )
	fl eval(sz t1, sz t2, fl r) const {
		return gaussian(r - (optimal_distance(t1, t2) + offset), width);
	}// e的-{(r-xs_vdw_radii[xs_t1]-xs_vdw_radii[xs_t2]-offset)/width}的平方
};

/*结构体repulsion继承usable
* 定义double   offset   
* 函数gauss给父类中的string参数name赋值   name被赋值 repulsion(o=     )
* 定义double类型函数eval  输入double类型t1   t2  r     输出double
	返回（r-xs_vdw_radii[t1]-xs_vdw_radii[t2]-offset）的平方
*/


struct repulsion : public usable {
	fl offset; // added to vdw
	repulsion(fl offset_, fl cutoff_) : usable(cutoff_), offset(offset_) {
		name = std::string("repulsion(o=") + to_string(offset) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const {
		fl d = r - (optimal_distance(t1, t2) + offset);//d=r-xs_vdw_radii[t1]-xs_vdw_radii[t2]-offset
		if(d > 0) //判断d>0
			return 0;
		return d*d;//返回d的平方
	}
};
//输入double   a、b、c    判断：c<a<b或者c<=b<=a返回0      a<b<c或者b<=a<=c返回1   其他的 返回double（c-a）/(b-a)
inline fl slope_step(fl x_bad, fl x_good, fl x) {   
	if(x_bad < x_good) {
		if(x <= x_bad) return 0;
		if(x >= x_good) return 1;
	}
	else {
		if(x >= x_bad) return 0;
		if(x <= x_good) return 1;
	}
	return (x - x_bad) / (x_good - x_bad);
}


/*结构体hydrophobic继承usable
* 定义double   good  bad
* 函数hydrophobic给父类中的string参数name赋值   name被赋值 hydrophobic(g=     ,b=     ,c=     )
* 定义double类型函数eval  输入double类型t1   t2  r     输出double  判断 t1=0或12或13或14或15  并且  t2=0或12或13或14或15 判断slope_step函数的条件
	返回  (（r-xs_vdw_radii[t1]-xs_vdw_radii[t2])-bad）/(good-bad)
*/
struct hydrophobic : public usable {
	fl good;
	fl bad;
	hydrophobic(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
		name = "hydrophobic(g=" + to_string(good) + ", b=" + to_string(bad) + ", c=" + to_string(cutoff) + ")";
	}//name被赋值 hydrophobic(g=     ,b=     ,c=     )
	fl eval(sz t1, sz t2, fl r) const {
		if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))//判断 t1=0或12或13或14或15  并且  t2=0或12或13或14或15
			return slope_step(bad, good, r - optimal_distance(t1, t2));//   (（r-xs_vdw_radii[t1]-xs_vdw_radii[t2])-bad）/(good-bad)
		else return 0;
	}
};


/*结构体non_hydrophobic继承usable
* 定义double   good  bad
* 函数hydrophobic给父类中的string参数name赋值   name被赋值 non_hydrophobic(g=     ,b=     ,c=     )
* 定义double类型函数eval  输入double类型t1   t2  r     输出double  判断 t1!=0或12或13或14或15  并且  t2!=0或12或13或14或15 判断slope_step函数的条件 
	返回  (（r-xs_vdw_radii[t1]-xs_vdw_radii[t2])-bad）/(good-bad)
*/
struct non_hydrophobic : public usable {
	fl good;
	fl bad;
	non_hydrophobic(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
		name = "non_hydrophobic(g=" + to_string(good) + ", b=" + to_string(bad) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const {
		if(!xs_is_hydrophobic(t1) && !xs_is_hydrophobic(t2))
			return slope_step(bad, good, r - optimal_distance(t1, t2));
		else return 0;
	}
};
/*定义模板参数 m n uint   输入  double    position  depth c_n   c_m  输出  改变参数c_n   c_m
//c_n = position的n次方  * depth * m / (fl(n)-fl(m));
//c_m = position的m次方  * depth * n/ (fl(m)-fl(n));
*/
template<unsigned n, unsigned m>
void find_vdw_coefficients(fl position, fl depth, fl& c_n, fl& c_m) {
	BOOST_STATIC_ASSERT(n != m); //确保n不等于m
	c_n = int_pow<n>(position) * depth * m / (fl(n)-fl(m));
	c_m = int_pow<m>(position) * depth * n / (fl(m)-fl(n));
}

/*定义模板参数 i  j uint 
结构体vdw继承usable
* 定义double   smoothing  cap
* 函数vdw给父类中的string参数name赋值   name被赋值 vdw(i=     ,j=     ,s=     ,^=     ,c=      )
* 定义double类型函数eval  输入double类型t1   t2  r    定义inter  double d0  depth  c_i  c_j  
输出double  min(cap, c_i / r的i次方 + c_j / r的j次方)或者cap
*/
template<unsigned i, unsigned j>
struct vdw : public usable {
	fl smoothing;
	fl cap;
	vdw(fl smoothing_, fl cap_, fl cutoff_) 
		: usable(cutoff_), smoothing(smoothing_), cap(cap_) {
		name = "vdw(i=" + to_string(i) + ", j=" + to_string(j) + ", s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const {
		fl d0 = optimal_distance(t1, t2);//xs_vdw_radii[t1]+xs_vdw_radii[t2]
		fl depth = 1; 
		fl c_i = 0;
		fl c_j = 0;
		find_vdw_coefficients<i, j>(d0, depth, c_i, c_j);//改变c_i  c_j
		if     (r > d0 + smoothing) r -= smoothing;//改变r
		else if(r < d0 - smoothing) r += smoothing;
		else r = d0;

		fl r_i = int_pow<i>(r);//r的i次方
		fl r_j = int_pow<j>(r);//r的j次方
		if(r_i > epsilon_fl && r_j > epsilon_fl)//r的i次方和j次方同时大于2.2204460492503131e-016
			return (std::min)(cap, c_i / r_i + c_j / r_j);//返回较小值
		else 
			return cap;
	}
};



/*结构体non_dir_h_bond继承usable
* 定义double   good  bad
* 函数non_dir_h_bond给父类中的string参数name赋值   name被赋值 non_dir_h_bond(g=     ,b=    )
* 定义double类型函数eval  输入double类型t1   t2  r     输出double  判(t1=3或5或7或9或16  并且t2=4或5或8或9)  或者  (t2=3或5或7或9或16  并且t1=4或5或8或9) 
判断slope_step函数的条件   返回  (（r-xs_vdw_radii[t1]-xs_vdw_radii[t2])-bad）/(good-bad)
*/
struct non_dir_h_bond : public usable {
	fl good;
	fl bad;
	non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
		name = std::string("non_dir_h_bond(g=") + to_string(good) + ", b=" + to_string(bad) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const {
		if(xs_h_bond_possible(t1, t2))//t1=3或5或7或9或16  并且t2=4或5或8或9  或者  t2=3或5或7或9或16  并且t1=4或5或8或9
			return slope_step(bad, good, r - optimal_distance(t1, t2));  //(（r-xs_vdw_radii[t1]-xs_vdw_radii[t2])-bad）/(good-bad)
		return 0;
	}
};
/*定义 函数read_iterator   输入   i  迭代器     输出double
* function   返回迭代器指向的容器的元素  迭代器指向下一个元素
*/
inline fl read_iterator(flv::const_iterator& i) {
	fl x = *i; 
	++i;
	return x;
}
/*定义 函数smooth_div   输入   x   y   double      输出double
* function   x的绝对值<2.2204460492503131e-016  返回0   y的绝对值<2.2204460492503131e-016   并且xy同号返回 1.7976931348623158e+308         异号返回- 1.7976931348623158e+308
否则返回  x/y
*/
fl smooth_div(fl x, fl y) {
	if(std::abs(x) < epsilon_fl) return 0;
	if(std::abs(y) < epsilon_fl) return ((x*y > 0) ? max_fl : -max_fl); // FIXME I hope -max_fl does not become NaN
	return x / y;
}//x的绝对值<2.2204460492503131e-016  返回0   y的绝对值<2.2204460492503131e-016   并且xy同号返回 1.7976931348623158e+308         异号返回- 1.7976931348623158e+308


/*结构体num_tors_add继承conf_independent
* 函数non_dir_h_bond给父类中的string参数name赋值   name被赋值 num_tors_add
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i 
* w=（迭代器指向的值）迭代器指向下一个元素
输出double  x + w * num_tors(结构体的参数);
*/

struct num_tors_add : public conf_independent {
	num_tors_add() { name = "num_tors_add"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		//fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		fl w = read_iterator(i); // FIXME?
		return x + w * in.num_tors;
	}
};


/*结构体num_tors_sqr继承conf_independent
* 函数num_tors_sqr给父类中的string参数name赋值   name被赋值 num_tors_sqr
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i
* * w=0.1*（迭代器指向的值）迭代器指向下一个元素
输出double  x + w * （num_tors的平方) / 5;
*/
struct num_tors_sqr : public conf_independent {
	num_tors_sqr() { name = "num_tors_sqr"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		return x + w * sqr(fl(in.num_tors)) / 5;
	}
};

/*结构体num_tors_sqrt继承conf_independent
* 函数num_tors_sqrt给父类中的string参数name赋值   name被赋值 num_tors_sqrt
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i
* w=0.1*（迭代器指向的值）迭代器指向下一个元素
输出double  x + w * （num_tors的平方根) /（5.0的平方根) ;
*/
struct num_tors_sqrt : public conf_independent {
	num_tors_sqrt() { name = "num_tors_sqrt"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		return x + w * std::sqrt(fl(in.num_tors)) / sqrt(5.0);
	}
};
/*结构体num_tors_div继承conf_independent
* 函数num_tors_div给父类中的string参数name赋值   name被赋值 num_tors_div
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i
w=0.1*（迭代器指向的值+1） 迭代器指向下一个元素
若x的绝对值<2.2204460492503131e-016  返回0   若(1+w*num_tors/5.0)的绝对值<2.2204460492503131e-016   并且xy同号返回 1.7976931348623158e+308         异号返回- 1.7976931348623158e+308
否则返回  x/ (1+w*num_tors/5.0)
*/
struct num_tors_div : public conf_independent {
	num_tors_div() { name = "num_tors_div"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.1 * (read_iterator(i) + 1); // w is in [0..0.2]
		return smooth_div(x, 1 + w * in.num_tors/5.0);
	}
};
/*结构体ligand_length继承conf_independent
* 函数ligand_length给父类中的string参数name赋值   name被赋值 ligand_length
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i
* w=(迭代器指向的值） 迭代器指向下一个元素
返回double    x + w * in.ligand_lengths_sum
*/
struct ligand_length : public conf_independent {
	ligand_length() { name = "ligand_length"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = read_iterator(i);
		return x + w * in.ligand_lengths_sum;
	}
};
/*结构体num_ligands继承conf_independent
* 函数num_ligands给父类中的string参数name赋值   name被赋值 num_ligands
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i
*  w=1*（迭代器指向的值） 迭代器指向下一个元素
* 返回x + w * in.num_ligands
*/
struct num_ligands : public conf_independent {
	num_ligands() { name = "num_ligands"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 1 * read_iterator(i); // w is in [-1.. 1]
		return x + w * in.num_ligands;
	}
};
/*结构体num_heavy_atoms_div继承conf_independent
* 函数num_heavy_atoms_div给父类中的string参数name赋值   name被赋值 num_heavy_atoms_div
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i
*   w=0.05*（迭代器指向的值） 迭代器指向下一个元素
若x的绝对值<2.2204460492503131e-016  返回0   若(1+w*num_heavy_atoms)的绝对值<2.2204460492503131e-016   并且xy同号返回 1.7976931348623158e+308         异号返回- 1.7976931348623158e+308
否则返回  x/ (1+w*num_heavy_atoms)
*/
struct num_heavy_atoms_div : public conf_independent {
	num_heavy_atoms_div() { name = "num_heavy_atoms_div"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.05 * read_iterator(i); 
		return smooth_div(x, 1 + w * in.num_heavy_atoms); 
	}
};
/*结构体num_heavy_atoms继承conf_independent
* 函数num_heavy_atoms给父类中的string参数name赋值   name被赋值 num_heavy_atoms
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i
* *   w=0.05*（迭代器指向的值） 迭代器指向下一个元素
返回double   x + w * num_heavy_atoms;
*/
struct num_heavy_atoms : public conf_independent {
	num_heavy_atoms() { name = "num_heavy_atoms"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.05 * read_iterator(i); 
		return x + w * in.num_heavy_atoms;
	}
};
/*结构体num_hydrophobic_atoms继承conf_independent
* 函数num_hydrophobic_atoms给父类中的string参数name赋值   name被赋值 num_hydrophobic_atoms
* 定义uint函数     返回1  只读
* 定义double类型函数eval  输入结构体  in   double类型  x    迭代器  i
* w=0.05*（迭代器指向的值） 迭代器指向下一个元素
* 返回double  x + w * in.num_hydrophobic_atoms
*/
struct num_hydrophobic_atoms : public conf_independent {
	num_hydrophobic_atoms() { name = "num_hydrophobic_atoms"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.05 * read_iterator(i); 
		return x + w * in.num_hydrophobic_atoms;
	}
};
/*定义 结构体everything的 构造函数everything  
* 定义 uint  d    double   cutoff
*/
everything::everything() { // enabled according to design.out227
	const unsigned d = 0; // default
	const fl cutoff = 8; //6;

	// FIXME? enable some?
	//// distance_additive
	//add(d, new ad4_solvation(3.6, 0.01097,  true, cutoff)); // desolvation_sigma, solvation_q, charge_dependent, cutoff
	//add(d, new ad4_solvation(3.6, 0.01097, false, cutoff)); // desolvation_sigma, solvation_q, charge_dependent, cutoff

	//add(d, new electrostatic<1>(100, cutoff)); // cap, cutoff
	//add(d, new electrostatic<2>(100, cutoff)); // cap, cutoff

	//add(d, new gauss(0,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(0.5, 0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1.5, 0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2.5, 0.3, cutoff)); // offset, width, cutoff

	add(1, new gauss(0, 0.5, cutoff));//在fun最末尾添加gauss，在enable最末尾添加0，1
	// offset, width, cutoff // WEIGHT: -0.035579
	//add(d, new gauss(1, 0.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.5, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 0.7, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 0.7, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.7, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 0.9, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(4, 1.5, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 2.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 2.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 2.0, cutoff)); // offset, width, cutoff
	add(1, new gauss(3, 2.0, cutoff)); //在fun最末尾添加gauss，在enable最末尾添加0，1
	// offset, width, cutoff // WEIGHT: -0.005156
	//add(d, new gauss(4, 2.0, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(4, 3.0, cutoff)); // offset, width, cutoff

	//add(d, new repulsion( 0.4, cutoff)); // offset, cutoff
	//add(d, new repulsion( 0.2, cutoff)); // offset, cutoff
	add(1, new repulsion( 0.0, cutoff));//在fun最末尾添加repulsion，在enable最末尾添加0，1
	// offset, cutoff // WEIGHT:  0.840245
	//add(d, new repulsion(-0.2, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.4, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.6, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.8, cutoff)); // offset, cutoff
	//add(d, new repulsion(-1.0, cutoff)); // offset, cutoff

	//add(d, new hydrophobic(0.5, 1, cutoff)); // good, bad, cutoff
	add(1, new hydrophobic(0.5, 1.5, cutoff)); //在fun最末尾添加hydrophobic，在enable最末尾添加0，1
	// good, bad, cutoff // WEIGHT:  -0.035069
	//add(d, new hydrophobic(0.5, 2, cutoff)); // good, bad, cutoff
	//add(d, new hydrophobic(0.5, 3, cutoff)); // good, bad, cutoff

	//add(1, new non_hydrophobic(0.5, 1.5, cutoff));

	//add(d, new vdw<4,  8>(   0, 100, cutoff)); // smoothing, cap, cutoff

	add(1, new non_dir_h_bond(-0.7, 0, cutoff)); //在fun最末尾添加non_dir_h_bond，在enable最末尾添加0，1
	// good, bad, cutoff // WEIGHT:  -0.587439
	//add(d, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff
	//add(d, new non_dir_h_bond(-0.7, 0.2, cutoff)); // good, bad, cutoff
	//add(d, new non_dir_h_bond(-0.7, 0.4, cutoff)); // good, bad, cutoff
	// additive

	// conf-independent
	//add(d, new num_ligands());

	add(1, new num_tors_div()); //在fun最末尾添加num_tors_div，在enable最末尾添加0，1
	// WEIGHT: 1.923 -- FIXME too close to limit?
	//add(d, new num_heavy_atoms_div());
	//add(d, new num_heavy_atoms());
	//add(1, new num_tors_add());
	//add(d, new num_tors_sqr());
	//add(d, new num_tors_sqrt());
	//add(d, new num_hydrophobic_atoms());
	///add(1, new ligand_length());

	//add(d, new num_tors(100, 100, false)); // cap, past_cap, heavy_only
	//add(1, new num_tors(100, 100,  true)); // cap, past_cap, heavy_only
	//add(d, new num_tors(  2,   1,  true)); // cap, past_cap, heavy_only
	//add(d, new num_heavy_atoms());
	//add(d, new ligand_max_num_h_bonds());
	//add(1, new num_ligands());
}
