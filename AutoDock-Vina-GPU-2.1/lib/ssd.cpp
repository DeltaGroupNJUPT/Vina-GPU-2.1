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

#include "ssd.h"

 /*  input：
			m——model结构体；p——precalculate结构体；ig——igrid结构体；out——output_type结构体；
			g——change结构体；v——vec结构体

			·eval_deriv（）函数作用:
			1.根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
			2.定义一个double类型e,具体查看输入ig中的eval_deriv定义,ig应该是igrid派生类对象
			3.将p,v[2], other_pairs, coords,minus_forces输入eval_interacting_pairs_deriv函数，返回值累加到e上
			4.遍历ligands的元素
					将p,v[0],ligands第i个元素的成员pairs,coord,minus_forces输入eval_interacting_pairs_deriv函数，返回值累加到e上
			5.将coords, minus_forces, g.ligands输入ligands.derivative中，分别对ligands的所有元素求导数，ligands[i].derivative(coords, minus_forces, g.ligands[i])
			6.将coords, minus_forces, g.flex输入flex.derivative中，分别对flex的所有元素求导数，flex[i].derivative(coords, minus_forces, g.flex[i])
			7.返回e（double）

 */

// clean up
void ssd::operator()(model& m, const precalculate& p, const igrid& ig, output_type& out, change& g, const vec& v) const { // g must have correct size
	out.e = m.eval_deriv(p, ig, v, out.c, g);			 //out结构体中的double类型e
	fl factor = initial_factor;							 //ssd结构体中的initial_factor变量
	VINA_U_FOR(i, evals) {								 //unsigned=unsigned int
		if(factor < min_factor) break;					 //如果initial_factor<min_factor
		output_type candidate(out);						 //out结构体复制给candidate，candidate做改变不会使out变化
		candidate.c.increment(g, -factor);				 //配体刚性变化中位置、方向、扭矩以及柔性中扭矩进行标准化
		change candidate_g(g);							 //g结构体复制给candidate_g，candidate_g做改变不会使g变化
		candidate.e = m.eval_deriv(p, ig, v, candidate.c, candidate_g);
		if(candidate.e <= out.e) {                
			out = candidate;                              //candidate赋给out
			g = candidate_g;				              //candidate_g赋给g	
			factor *= up;                                 //factor=factor* up
		}
		else {
			factor *= down;                              //factor=factor* down
		}
	}
	out.coords = m.get_heavy_atom_movable_coords();      //从i=0遍历到m_num_movable_atoms(model结构体中),依次判断atoms（atom结构体向量）成员变量el是否等于0，
													     //是则往out.coords尾部添加model成员coords容器第i个元素
}
