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

#include "szv_grid.h"
#include "brick.h"
/*
定义构造函数szv_grid  初始化赋值  为结构体里的参数 m_init    m_range   m_data赋值
*/
szv_grid::szv_grid(const model& m, const grid_dims& gd, fl cutoff_sqr) : m_data(gd[0].n, gd[1].n, gd[2].n) {
	vec end;
	VINA_FOR_IN(i, gd) {
		m_init[i] = gd[i].begin;
		end   [i] = gd[i].end;
	}//将grid_dim的double参数begin  和end   赋值
	m_range = end - m_init;//为m_range赋值

	const sz nat = num_atom_types(m.atom_typing_used());	//atom_typing_used枚举类型参数 //num_atom_types 得到枚举对应返回值
	

 
//若get函数得到的t对应返回值 uint  <   num_atom_types函数得到的t对应返回值 uint  并且
//三个数组经过函数brick_distance_sqr得到的数 double < cutoff_sqr
//那么在 relevant_indexes后边加元素i 
	szv relevant_indexes;  //uint型的向量
	VINA_FOR_IN(i, m.grid_atoms) {//grid_atoms  结构体atom类型的向量   
		const atom& a = m.grid_atoms[i];
		if(a.get(m.atom_typing_used()) < nat && brick_distance_sqr(m_init, end, a.coords) < cutoff_sqr)//get   得到t对应返回值uint      atom_typing_used返回枚举型t类型数据
			relevant_indexes.push_back(i);
	}
	/*循环执行 
	x y z为    000   【  a=grid_atoms【relevant_indexes0】  判断  三个数组经过函数brick_distance_sqr后是否<cutoff_sqr
	                                                    m_data[i + m_i*(j + m_j*k)]后加入元素relevant_indexes0
			   000   【  a=grid_atoms【relevant_indexes1】     判断  三个数组经过函数brick_distance_sqr后是否<cutoff_sqr
														m_data[i + m_i*(j + m_j*k)]后加入元素relevant_indexes1        】
			   000   【  a=grid_atoms【relevant_indexes2】     判断  三个数组经过函数brick_distance_sqr后是否<cutoff_sqr
														m_data[i + m_i*(j + m_j*k)]后加入元素relevant_indexes2        】
	x y z为	   001  .....
			   002  .....
			   00m_k  ....
			   m_im_jm_k....

   */
	VINA_FOR(x, m_data.dim0())
	VINA_FOR(y, m_data.dim1())
	VINA_FOR(z, m_data.dim2()) {
		VINA_FOR_IN(ri, relevant_indexes) {
			const sz i = relevant_indexes[ri];
			const atom& a = m.grid_atoms[i];
			if(brick_distance_sqr(index_to_coord(x, y, z), index_to_coord(x+1, y+1, z+1), a.coords) < cutoff_sqr)
				m_data(x, y, z).push_back(i);
		}
	}
}
/*定义double 函数average_num_possibilities   
* 
循环执行  xyz从000  到m_im_jm_k     counter += m_data(x, y, z).size()
返回double类型    counter /（m_i*m_j*m_k）

*/
fl szv_grid::average_num_possibilities() const {
	sz counter = 0;
	VINA_FOR(x, m_data.dim0())
	VINA_FOR(y, m_data.dim1())
	VINA_FOR(z, m_data.dim2()) {
		counter += m_data(x, y, z).size();
	}
	return fl(counter) / (m_data.dim0() * m_data.dim1() * m_data.dim2());
}
/*定义只读函数possibilities  返回uint型向量    输入长度3的数组  coords
* inter     长度3uint类型数组  index
* function  ：给index[i]赋值    返回   m_data[index[0] + m_i*(index[1] + m_j*index[2])]      m_data是array3d<szv>类型
*/
const szv& szv_grid::possibilities(const vec& coords) const {
	boost::array<sz, 3> index;
	VINA_FOR_IN(i, index) {//循环3次    
		assert(coords[i] + epsilon_fl >= m_init[i]);//判断条件  epsilon_fl  ：2.2204460492503131e-016
		assert(coords[i] <= m_init[i] + m_range[i] + epsilon_fl);
		const fl tmp = (coords[i] - m_init[i]) * m_data.dim(i) / m_range[i];
		index[i] = fl_to_sz(tmp, m_data.dim(i) - 1);//若0<tmp<m_data.dim(i) - 1则 强制整形并输出
	}
	return m_data(index[0], index[1], index[2]);//m_data是array3d<szv>类型  返回m_data[index[0] + m_i*(index[1] + m_j*index[2])]
}

/*input     uint  i j  k
inter   index  tmp   长度3的数组
output  长度3的数组tmp
function  tmp[i] = m_init[i] + m_range[i] * index[i] / m_data.dim(i);    i= 0、1、2
其中   index(i)= i  j   k   dim(i)=m_i;m_j;m_k;
*/
vec szv_grid::index_to_coord(sz i, sz j, sz k) const {
	vec index(i, j, k);
	vec tmp;
	VINA_FOR_IN(n, tmp) 
		tmp[n] = m_init[n] + m_range[n] * index[n] / m_data.dim(n);
	return tmp;
}
/*定义grid_dims类型的函数  szv_grid_dims    输入  只读grid_dims类型的参数  gd   输出  grid_dims类型的tmp
inter： 定义grid_dims类型参数tmp   里边包含三个grid_dim类型结构体
function：给tmp[i]的参数begin   end    n赋值  然后返回tmp
*/
grid_dims szv_grid_dims(const grid_dims& gd) {
	grid_dims tmp;    
	VINA_FOR_IN(i, tmp) {
		tmp[i].begin = gd[i].begin;//grid_dims中的double类型参数begin和end传递
		tmp[i].end   = gd[i].end;
		fl n_fl = (gd[i].end - gd[i].begin) / 3; // 3A preferred size
		int n_int = int(n_fl);
		tmp[i].n     = (n_int < 1) ?  1 : sz(n_int);//给tmp的uint类型参数n赋值
	}
	return tmp;
}

