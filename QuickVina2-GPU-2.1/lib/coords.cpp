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

#include "coords.h"
/*
	input——a、b：vec结构体向量；
	inter——acc:double;
	demo:前提a、b元素个数一样，并且不为空，否则返回0
	a{（7，8，9），（5，6，7）}，b{（1，2，3），（3，4，5）} size=2
	acc=（(7-1)^2+(8-2)^2+(9-3)^2）+（(5-3)^2+(6-4)^2+(7-5)^2）=120
	返回 （120/2）^(1/2)=7.74597
*/
fl rmsd_upper_bound(const vecv& a, const vecv& b) {
	VINA_CHECK(a.size() == b.size());
	fl acc = 0;
	VINA_FOR_IN(i, a)       //for(i=0;i<a.size;i++)
		acc += vec_distance_sqr(a[i], b[i]);
	return (a.size() > 0) ? std::sqrt(acc / a.size()) : 0;
}


/*
	input——a、b.coords：vec结构体向量，b：output_type结构体向量
	output——tmp：<unint,1.79769e+308>,可以用first、second调用
	demo:前提a、b.coords元素个数一样，并且不为空，否则返回res=0
	 a{（7，8，9），（5，6，7）}，b.coords{（1，2，3），（3，4，5）} size=2
	res=[(（(7-1)^2+(8-2)^2+(9-3)^2）+（(5-3)^2+(6-4)^2+(7-5)^2）)/2]^(1/2)=7.74597
	返回tmp<2,7.74597>
*/

std::pair<sz, fl> find_closest(const vecv& a, const output_container& b) {
	std::pair<sz, fl> tmp(b.size(), max_fl);          //max_fl=1.79769e+308
	VINA_FOR_IN(i, b) { 
		fl res = rmsd_upper_bound(a, b[i].coords);
		if(i == 0 || res < tmp.second)
			tmp = std::pair<sz, fl>(i, res);
	}
	return tmp;
}


/*
	input——out：output_type结构体向量，t：output_type结构体，min_rmsd：double，max_size：unint
	根据find_closest函数的demo，closest_rmsd<2,7.74597>
	·如果out向量元素>2，min_rmsd>7.74597（表示有一个非常小的值）并且t.e < out[2].e
		out[2]=t;
	·否则（表示没有小的值）
	如果out向量元素<max_size,在out向量末尾加上t（output_type）结构体
	否则如果out向量不为空并且t.e < out末尾元素的e值，将out末尾元素赋值为t
	·最后将out向量以out结构体中e进行升序排序

*/
void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size) {
	std::pair<sz, fl> closest_rmsd = find_closest(t.coords, out);
	if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
		if(t.e < out[closest_rmsd.first].e) { // the new one is better, apparently
			out[closest_rmsd.first] = t; // FIXME? slow
		}
	}
	else { // nothing similar
		if(out.size() < max_size)   //new开辟的空间在堆上，而一般声明的变量存放在栈上
			out.push_back(new output_type(t)); // the last one had the worst energy - replacing 
		else
			if(!out.empty() && t.e < out.back().e) // FIXME? - just changed
				out.back() = t; // FIXME? slow
	}
	out.sort();
}
