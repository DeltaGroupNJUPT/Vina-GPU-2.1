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

#ifndef VINA_MONTE_CARLO_H
#define VINA_MONTE_CARLO_H

#include "ssd.h"
#include "incrementable.h"
#include "commonMacros.h"
/*
* 结构体monte_carlo
* 成员：1.unsigned类型num_steps
*		2.double类型temperature
*		3.vec矢量hunt_cap
*		4.double类型min_rmsd
*		5.unsigned类型num_saved_mins
*		6.double类型mutation_amplitude
*		7.ssd类对象ssd_par
*		8.构造函数：初始化num_steps=2500，temperature=1.2，hunt_cap=[10, 1.5, 10], min_rmsd=0.5, num_saved_mins=50, mutation_amplitude=2
*       9.声明定义调用操作符（），为常量函数
*					input:model类引用,precalculate类常量引用,igrid类常量引用,precalculate类常量引用,igrid类常量引用,vec矢量常量引用,vec矢量常量引用,incrementable类指针,generator类引用
*					output:output_type类
*		9.声明many_runs函数，为常量函数
* 					input:model类引用,precalculate类常量引用,igrid类常量引用,vec矢量常量引用,vec矢量常量引用,unsigned int类型数,generator类引用
* 					output:output_type类
*		10.声明single_run函数，为常量函数
* 					input:model类引用,output_type类，precalculate类常量引用,igrid类常量引用,generator类引用
* 					output:无
*		11.重载声明定义调用操作符（），为常量函数
* 					input:model类引用,output_container类引用，precalculate类常量引用,igrid类常量引用,precalculate类常量引用,igrid类常量引用,vec矢量常量引用,vec矢量常量引用,incrementable类指针,generator类引用
* 					output:无
*					相比原版本输入多了output_container类引用，无返回值
* 		9.声明many_runs函数，为常量函数
* 					input:model类引用,output_container类引用，precalculate类常量引用,igrid类常量引用,vec矢量常量引用,vec矢量常量引用,unsigned int类型数,generator类引用
* 					output:无
*					相比原版本输入多了output_container类引用，无返回值
*/
struct monte_carlo {
	unsigned num_steps;
	fl temperature;
	vec hunt_cap;
	fl min_rmsd;
	sz num_saved_mins;
	fl mutation_amplitude;
	ssd ssd_par;

	int thread; //parallelism 20210917 Glinttsd
	std::vector<int> search_depth; // 20210813 Glinttsd



	monte_carlo() : num_steps(2500), temperature(1.2), hunt_cap(10, 1.5, 10), min_rmsd(0.5), num_saved_mins(50), mutation_amplitude(2) {} // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  num_steps = 50*lig_atoms = 2500

	output_type operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
	output_type many_runs(model& m, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const;

	void single_run(model& m, output_type& out, const precalculate& p, const igrid& ig, rng& generator) const;
	// out is sorted
	void operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
	void many_runs(model& m, output_container& out, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const;
	std::vector<output_type> cl_to_vina(output_type_cl result_ptr[], int exhaus) const;
	void generate_uniform_position(const vec corner1, const vec corner2, std::vector<vec>& uniform_data, int exhaustiveness) const;
	//void print_process(boost::progress_display* p_d);
};

#endif

