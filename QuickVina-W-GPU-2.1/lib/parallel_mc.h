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

#ifndef VINA_PARALLEL_MC_H
#define VINA_PARALLEL_MC_H

#include "monte_carlo.h"


/*
* 结构体parallel_mc
* 成员：1.monte_carlo类mc
*		2.unsigned int类型num_tasks
*		3.unsigned int类型num_threads
*		4.bool类型display_progress
*		5.构造函数：初始化num_tasks=8, num_threads=1, display_progress为真 	
*		6.定义调用操作符（），为常量函数
* 					input:model类常量引用,output_container类引用，precalculate类常量引用,igrid类常量引用,precalculate类常量引用,igrid类常量引用,vec矢量常量引用,vec矢量常量引用,generator类引用
* 				    output:无					
*/
struct parallel_mc {
	monte_carlo mc;
	sz num_tasks;
	sz num_threads;
	bool display_progress;
	parallel_mc() : num_tasks(8), num_threads(1), display_progress(true) {}
	void operator()(const model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const;
};

#endif
