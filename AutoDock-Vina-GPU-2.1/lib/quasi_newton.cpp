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

#include "quasi_newton.h"
#include "bfgs.h"



/*
      m——model结构体；p——precalculate结构体指针；ig——igrid结构体指针；
	  v——vec结构体；c——conf结构体；g——change结构体；
	 ·eval_deriv（）函数作用:
			1.根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
			2.定义一个double类型e,具体查看输入ig中的eval_deriv定义,ig应该是igrid派生类对象
			3.将p,v[2], other_pairs, coords,minus_forces输入eval_interacting_pairs_deriv函数，返回值累加到e上
			4.遍历ligands的元素
					将p,v[0],ligands第i个元素的成员pairs,coord,minus_forces输入eval_interacting_pairs_deriv函数，返回值累加到e上
			5.将coords, minus_forces, g.ligands输入ligands.derivative中，分别对ligands的所有元素求导数，ligands[i].derivative(coords, minus_forces, g.ligands[i])
			6.将coords, minus_forces, g.flex输入flex.derivative中，分别对flex的所有元素求导数，flex[i].derivative(coords, minus_forces, g.flex[i])
			7.返回e（double）
	  ·operator（）函数
			input——c：conf结构体，g：change结构体
			用于返回eval_deriv（）函数返回值


*/
struct quasi_newton_aux {
	model* m;
	const precalculate* p;
	const igrid* ig;           
	const vec v;
	quasi_newton_aux(model* m_, const precalculate* p_, const igrid* ig_, const vec& v_) : m(m_), p(p_), ig(ig_), v(v_) {}
	fl operator()(const conf& c, change& g) {
		const fl tmp = m->eval_deriv(*p, *ig, v, c, g);
		return tmp;                                    //e
	}
};



/*
		m——model结构体；p——precalculate结构体；ig——igrid结构体；out——output_type；
		g——change结构体；v——vec结构体
		out.e=res=bfgs函数返回值f0
*/
void quasi_newton::operator()(model& m, const precalculate& p, const igrid& ig, output_type& out, change& g, const vec& v) const { // g must have correct size
	quasi_newton_aux aux(&m, &p, &ig, v);                                      //aux为重构quasi_newton_aux结构体
	fl res = bfgs(aux, out.c, g, max_steps, average_required_improvement, 10); //f0
	out.e = res;                                                               
}

