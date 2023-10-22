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

#ifndef VINA_SZV_GRID_H
#define VINA_SZV_GRID_H

#include "model.h"
#include "grid_dim.h"
#include "array3d.h"
/*结构体szv_grid
*声明函数szv_grid  输入  类 只读  m   gd  ， double   cutoff_sqr
* 声明uint向量类型 函数 possibilities    只读      输入为  向量coords
* 声明函数average_num_possibilities只读
* 定义类型参数 模板参数uint的向量   m_data
* 定义长度为3的数组m_init   m_range  
* 定义函数index_to_coord  输入参数为uint  i，j，k   只读  返回长度为3的数组
* 声明一个结构体组（3个结构体）型的函数szv_grid_dims
*/
struct szv_grid {
	szv_grid(const model& m, const grid_dims& gd, fl cutoff_sqr);
	const szv& possibilities(const vec& coords) const;
	fl average_num_possibilities() const;
public:
	array3d<szv> m_data;
	vec m_init;
	vec m_range;
	vec index_to_coord(sz i, sz j, sz k) const;
};

grid_dims szv_grid_dims(const grid_dims& gd);


#endif
