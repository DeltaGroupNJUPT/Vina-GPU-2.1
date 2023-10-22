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

#include "manifold.h"
#include "recent_history.h"
#include "coords.h"
#include "quasi_newton.h"

const vec authentic_v(1000, 1000, 1000); 
/*定义操作符重载    
input：m  model    p：precalculate结构体  ig、ig_widened：igrid结构体   p_widened：p_widened结构体  corner1、corner2：向量   generator:随机数产生器
out：output_type
function：定义tmp   output_container    然后当前对象执行parallel_mc::operator()  给tmp赋值    判断tmp不为空    并且返回tmp的首个元素
*/
output_type manifold::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, generator); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}
/*定义结构体指针类型函数  select_rs     只读
* input   output_container  结构体向量  mf  只读     rstart：uint      随机数产生器 generator
*判断结构体向量不空并且元素个数大于等于rstart     则  返回这个结构体向量中的第i个结构体地址     i为在0-（元素个数-1）之间的随机数 i 
* 否则 返回空指针
*/
const conf* select_rs(const output_container& mf, sz rstart, rng& generator) { // clean up
	if(!mf.empty() && mf.size() >= rstart) { //向量不空并且元素个数大于等于rstart
		sz i = random_sz(0, mf.size()-1, generator);//产生一个在0-（元素个数-1）之间的随机数 i uint类型
		return &(mf[i].c);//返回这个结构体向量中的第i个结构体地址
	}
	return NULL;//否则指针为空
}
/*
* 判断结构体向量中的内部配体变化是否相似   并且返回最后一次结果的相反值
input：conf类型结构体C   output_container类型结构体向量mf    布尔类型向量internal_too_close    double  ：exclusion（用于判断相似函数internal_too_close中作为判断标准）
out： tmp  bool
判断结构体向量mf[i]中的结构体conf中的结构体向量torsions是否和ligands结构体向量中的torsions相似，并且判断结果给internal_too_close赋值
最后一次判断的结果是相似的话   返回fause   否则返回ture

*/
bool determine_iu(const conf& c, const output_container& mf, std::vector<bool>& internal_too_close, fl exclusion) { // clean up
	bool tmp = true;
	VINA_FOR_IN(i, mf) {
		bool t = c.internal_too_close(mf[i].c, exclusion);/* 判断ligands结构体向量中的torsions是否和mf[i].conf.结构体向量中的torsions相似是的话返回true，否则返回false */
		internal_too_close[i] = t;//判断结果给向量赋值
		if(t) 
			tmp = false;//最后一次判断的结果是相似的话   返回fause   否则返回ture
	}
	return tmp;
}
/*

判断结构体向量中的内部配体变化和外部配体变化是否都相似  是的话返回false  否则返回ture
input  c  conf结构体   mf  output_container  ：结构体向量   internal_too_close  bool向量  （内部配体相似的判断结果） exclusion：scale结构体
out    bool 

external_too_close：判断外部匹配变化是否相似   ：先判断ligands[i].rigid中的position、orientation是否和c.ligands[i].rigid中的position、orientation相似，不是的话返回false；
否则判断flex[i]中的torsions向量和c.flex[i]中的torsions向量是否相似，不是的话返回false；
要不最后返回true
*/

bool conf_is_legal(const conf& c, const output_container& mf, const std::vector<bool>& internal_too_close, const scale& exclusion) {
	assert(mf.size() == internal_too_close.size());// 判断mf结构体向量的元素个数是否和internal_too_close向量的元素个数相等
	VINA_FOR_IN(i, mf)
		if(internal_too_close[i] && c.external_too_close(mf[i].c, exclusion))//判断内部变化和外部变化是否都相似
			return false;
	return true;
}
/*input   internal_conf：conf结构体只读     mf： 结构体向量  只读    internal_too_close：布尔类型向量 只读   uniq： bool    spread、exclusion ： scale结构体 只读     rp：double  rs： conf结构体指针只读  num_attempts：uint   failed ：bool    generator：随机数产生器
out： tmp  conf结构体 
function：  把输入internal_conf过渡给tmp    
			循环执行             num_attempts次  
						1、对tmp生成外部变化：刚性（位置、方向、扭矩）、柔性（扭矩）
						2、如果uniq为0且tmp内部配体外部配体都相似  返回tmp   ，
						3、如果执行到(attempt + 1 >= num_attempts)  则返回tmp  且failed变成ture
			若接到错误，则返回internal_conf
*/
conf generate_external(const conf& internal_conf, const output_container& mf, const std::vector<bool>& internal_too_close, bool uniq, const scale& spread, const scale& exclusion, fl rp, const conf* rs, unsigned num_attempts, bool& failed, rng& generator) {
	// FIXME if one ligand with no side chains, don't use the same (pos, orient) more than once 
	failed = false;
	VINA_U_FOR(attempt, num_attempts) {//for  （attempt=0；attempt<num_attempts;attempt++）
		conf tmp(internal_conf);
		tmp.generate_external(spread, rp, rs, generator); //生成外部变化：刚性（位置、方向、扭矩）、柔性（扭矩）。
		if(uniq || conf_is_legal(tmp, mf, internal_too_close, exclusion))//如果uniq为0且tmp内部配体外部配体都相似  返回tmp
			return tmp;
		if(attempt + 1 >= num_attempts) { 
			failed = true;
			return tmp;
		}
	}
	VINA_CHECK(false); 
	return internal_conf; // placating the compiler
}
/*input   double   from   to    progress
* out    double
* function：判断 from 是否< 2.2204460492503131e-016  
是   return  from   
否则   return    from * e的(progress * log(to/from))次方
*/
fl extrapolate_cap(fl from, fl to, fl progress) {
	if(from < epsilon_fl) return from;
	return from * std::exp(progress * std::log(to/from));
}


/*input   from to： 只读 double向量   progress   ：double
out   double向量
function：  跟tmp[i]赋值并返回  
判断 from[i] 是否< 2.2204460492503131e-016 
是  tmp[i]=from[i] 
否  tmp[i]=from[i] * e的(progress * log(to[i]/from[i]))次方
*/
vec extrapolate_cap(const vec& from, const vec& to, fl progress) {
	vec tmp;
	VINA_FOR_IN(i, tmp)
		tmp[i] = extrapolate_cap(from[i], to[i], progress);
	return tmp;
}
/*
* 输入:   phase  :  uint , out:引用output_type结构体 , m：引用类model, mf：引用output_container结构体向量 只读,  p：引用结构体precalculate 只读,  ig：引用结构体igrid 只读, corner2corner、init_manifold_factor：double , e_internal_stats：引用结构体manifold, par：引用结构体manifold只读, rng& generator：随机数产生器)
  输出：   void

function: 循环执行  step从0-num_steps
  ｛
		1、复制输入的函数out生成candidate，然后对candidate的配置conf执行generate_internal函数用来生成配体内部变化：位置、方向、扭矩。
		2、判断结构体向量中的内部配体变化是否相似结果存在向量internal_too_close中   并且返回最后一次结果的相反值
		3、对配体执行seti函数：遍历c.ligands的元素，将model成员internal_coords, coords, ligands[i].begin（ligand成员）, ligands[i].end（ligand成员）输入函数c.ligands[i].rigid.apply，获得坐标coords   根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords

	                       其中：generate_internal函数的输入参数为：spread.torsion, rp, rs, generator。

								spread（scale类型结构体）：参数为（orner2corner（输入） * manifold_factor,pi * manifold_factor,pi * manifold_factor）
								                         manifold_factor=init_manifold_factor(输入) * std::exp(- par.manifold_lambda * step / par.num_steps)
								rp                     ： 1 - std::exp(std::log(1-par.max_prob) * phase/par.num_phases)

								rs                     ： mf结构体向量中的 第随机个结构体地址

         4、判断rigid_trials>0
		 5、循环执行   i从0-rigid_trials
										 ｛

											1、定义bool类型  failed =false
											2、定义candidate2（candidate一样类型的结构体）并且给其中的配体conf赋值并改变failed
																	赋值方法：


											3、判断： 如果failed=ture  且 （candidate_is_set=ture   或者 i+1 < rigid_trials ）  跳过本次循环
												  否则执行函数sete：遍历c.ligands的元素，将model成员internal_coords, coords, ligands[i].begin（ligand成员）, ligands[i].end（ligand成员）输入函数c.ligands[i].rigid.apply，获得坐标coords   根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords

											4、对candidate2的e赋值，赋值方法：1.定义一个double类型e,具体查看输入ig中的eval定义,ig应该是igrid派生类对象
																		 2.遍历ligands的元素 将p,v[2],other_pairs, coords输入eval_interacting_pairs函数，返回值累加到e上

											5、判断： 如果candidate_is_set为1并且  candidate2.e < candidate.e  那么不执行
												   否则 执行  ｛candidate_is_set=ture  candidate = candidate2 ｝
										   ｝

		6、检测candidate_is_set    if（ture）
					                    ｛    
					                      1、定义 	fl e_internal = (3 * step < 2 * num_steps) ? 0 : m.evali(p, hunt_cap);
																	其中  evali作用：1.定义一个double类型e=0
*																	                2.遍历ligands的元素   将p,v[0],ligands第i个元素的成员pairs,internal_coords输入eval_interacting_pairs函数，返回值累加到e上
*																	                3.返回e
					
					                       2、在add函数中输入参数e_internal    来改变e_internal_stats（输入）里的error_estimate_sqr和x_estimate 

																	this_error_sqr=（e_internal - x_estimate）^2;
																	
																	error_estimate_sqr=weight * this_error_sqr + (1 - weight) * error_estimate_sqr
																	x_estimate         = weight * x              + (1 - weight) * x_estimate;

											3、改变candidate中的参数e    candidate.e += e_internal;
					
											4、判断：   if(step == 0 || candidate.e < out.e)   out = candidate;

						                     ｝

	｝
		
		
		*/



void manifold_phase(sz phase, output_type& out, model& m, const output_container& mf, const precalculate& p, const igrid& ig, fl corner2corner, fl init_manifold_factor, recent_history& e_internal_stats, const manifold& par, rng& generator) { // out.first is starting conf on input
	out.e = max_fl; // 结构体里的e赋值

	const fl rp = 1 - std::exp(std::log(1-par.max_prob) * phase/par.num_phases); // 
	//rp=1-e的(log(1-max_prob)*phase/num_phases)次方
	fl e_internal_cost = par.relative_pair_cost * m.num_internal_pairs();//e_internal_cost=relative_pair_cost*pairs向量的元素个数
	fl e_external_cost = par.relative_pair_cost * m.num_other_pairs() + m.num_movable_atoms();
	const sz rigid_trials = 1 +  fl_to_sz(1 + par.cost_factor * e_internal_cost / e_external_cost, 100);
	//前者小于100大于0 则强制整形   用于下边的判断执行
	unsigned num_steps = (std::max)(unsigned(1), par.num_steps);

	VINA_U_FOR(step, num_steps) {//循环num_steps
		const fl manifold_factor = init_manifold_factor * std::exp(- par.manifold_lambda * step / par.num_steps);
		//用于下边的scale结构体的初始化
		vec hunt_cap; hunt_cap = extrapolate_cap(par.hunt_cap, authentic_v, fl(step) / num_steps);

		const scale spread(corner2corner * manifold_factor,
			               pi            * manifold_factor,
						   pi            * manifold_factor);//定义scale结构体  并初始化   用到了for循环的数据manifold_factor

		sz rstart = fl_to_sz(par.rstart_fraction * par.num_phases, par.num_phases);
		//用于构成结构体指针rs    用于函数判断执行用的
		const conf* rs = select_rs(mf, rstart, generator);
		//返回mf结构体向量中的 第随机个结构体地址
		output_type candidate = out;
		candidate.c.generate_internal(spread.torsion, rp, rs, generator); // 生成配体内部变化：位置、方向、扭矩。

		std::vector<bool> internal_too_close(mf.size());//和mf一样长度的向量
		bool uniq = determine_iu(candidate.c, mf, internal_too_close, par.exclusion.torsion); // internal_too_close is output
		//判断结构体向量中的内部配体变化是否相似   并且返回最后一次结果的相反值
		m.seti(candidate.c);
/*SETI
* 定义model的成员函数seti
* input:conf类常量引用c
* ouput:void
* 作用：将model成员atoms, internal_coords, 和c成员ligands输入ligands的成员函数set_conf    (*this)[i].set_conf(atoms, coords, c[i]);
*       根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和internal_coords
*/



		VINA_CHECK(rigid_trials > 0);
		bool candidate_is_set = false;//
		VINA_FOR(i, rigid_trials) {
			bool failed = false;
			output_type candidate2(generate_external(candidate.c, mf, internal_too_close, uniq, spread, par.exclusion, rp, rs, par.num_attempts, failed, generator), 0);
			//可以改变  generate_external： failed   并且返回conf结构体        而本句是  初始化   candidate2
			if(failed && 
				(candidate_is_set || i+1 < rigid_trials)) // candidate is set or a chance to set it remains
				continue;//如果failed=ture  且 （candidate_is_set=ture   或者 i+1 < rigid_trials ）  跳过本次循环
			m.sete(candidate2.c);

			/*SETE
* 定义model的成员函数sete
* input:conf类常量引用c
* ouput:void
* 作用：遍历c.ligands的元素，将model成员internal_coords, coords, ligands[i].begin（ligand成员）, ligands[i].end（ligand成员）输入函数c.ligands[i].rigid.apply，获得坐标coords
* 根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords				
*/
			
			candidate2.e = m.evale(p, ig, hunt_cap);  

			/*
* model的成员函数evale
* input:precalculate类常量引用p，igrid类常量引用ig，vec类常量引用v
* output:double类型数
* 作用：1.定义一个double类型e,具体查看输入ig中的eval定义,ig应该是igrid派生类对象
*		2.遍历ligands的元素
*				将p,v[2],other_pairs, coords输入eval_interacting_pairs函数，返回值累加到e上
*		3.返回e 
*/
			if(!candidate_is_set || candidate2.e < candidate.e) {
				candidate_is_set = true;
				candidate = candidate2; // this works because internal conf is the same in both // FIXME inefficient, all this
			}//如果candidate_is_set为1并且  candidate2.e < candidate.e  那么不执行   否则执行  ｛candidate_is_set=ture  candidate = candidate2 ｝
		}

		VINA_CHECK(candidate_is_set);
		if(true) {
			fl e_internal = (3 * step < 2 * num_steps) ? 0 : m.evali(p, hunt_cap);
			/*EVALI
* model的成员函数evali
* input:precalculate类常量引用p，vec类常量引用v
* output:double类型数
* 作用：1.定义一个double类型e=0
*		2.遍历ligands的元素
*				将p,v[0],ligands第i个元素的成员pairs,internal_coords输入eval_interacting_pairs函数，返回值累加到e上
*		3.返回e
*/
			e_internal_stats.add(e_internal);//改变e_internal_stats里的error_estimate_sqr和x_estimate
			/*add：
	input——x:double；
	inter——this_error_sqr：double；
		this_error_sqr=（x - x_estimate）^2;
		error_estimate_sqr=weight * this_error_sqr + (1 - weight) * error_estimate_sqr
		x_estimate         = weight * x              + (1 - weight) * x_estimate;
*/
			candidate.e += e_internal;//改变e
			if(step == 0 || candidate.e < out.e)
				out = candidate;
		}
	}
}
/*input：	 m,					引用model类型结构体
			out,				引用结构体向量output_container
			p,                  引用结构体precalculate   只读
			ig_widened、 ig,    引用结构体igrid   只读
			p_widened,          引用结构体precalculate   只读
			corner2 corner1,    引用向量   只读
			rng& generator      随机数产生器
	out：  void

	function：
				1、定义manifold类型结构体tuning_par   并把authentic_v赋值给器其参数hunt_cap
				2、定义corner2corner  用于给manifold_phase的输入
				3、定义conf_size类型s和change类型结构体gS   s赋值给g
				4、定义recent_history类型e_internal_stats  并赋值
	            5、执行循环  phase从0-np
							｛      1、定义output_type类型tmp并初始化
									2、tmp里的配体conf执行函数randomize，参数(corner1, corner2, generator)，生成刚性的随机位置、方向、扭矩和柔性中的随机扭矩
									3、执行manifold_phase函数
									4、判断：   如果use_ssd（结构体manifold的参数）为1，那么对结构体ssd_par执行重载函数   改变tmp的e参数  tmp.e=res=bfgs函数返回值f0
											   否则把结构体ssd_par的evals参数赋值给 quasi_newton_par（quasi_newton类型的结构体）里面的max_steps参数
											   然后对quasi_newton_par执行重载函数   改变tmp的e参数  tmp.e=res=bfgs函数返回值f0
									5、tmp的配体conf执行函数set          根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
									6、给tmp的coords赋值    执行get_heavy_atom_movable_coords，常量函数
									7、复制结构体find_closest给closest_rmsd
									8、判断：//若 closest_rmsd第一个参数<out的元素个数且 第二个参数<min_rmsd   且tmp.e < out[closest_rmsd.first].e
											//则out[closest_rmsd.first] = tmp
											//否则在out的后面加入tmp
							｝
						6、	定义output_container类型 null_array只读
						7、判断out不为空    并对其参数排序   判断排序正确

	*/


void manifold::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const {
	manifold tuning_par;
	tuning_par.hunt_cap = authentic_v;//给tuning_par的向量赋值

	fl corner2corner = 0.5 * (corner2 - corner1).norm();//  norm：sqrt（sqr(data[0]) + sqr(data[1]) + sqr(data[2])）
	//用于给manifold_phase赋值
	conf_size s = m.get_size();//定义S   用于给  g赋值
/*get_size
* 作用： 1.获得ligands的扭转计数保存在临时tmp.ligands对象中
*       2.获得flex的扭转计数保存在临时tmp.flex对象中
*		3.返回tmp赋值给  s
*/
	change g(s);//定义change类型结构体并初始化

	sz np = (std::max)(sz(1), num_phases);//用于循环  取1和num_phases最大的
	recent_history e_internal_stats(0, 10, 10);//定义结构体recent_history并且初始化    用于 给manifold_phase赋值



	VINA_FOR(phase, np) {
		output_type tmp(s, 0);
		tmp.c.randomize(corner1, corner2, generator);
/*生成刚性的随机位置、方向、扭矩和柔性中的随机扭矩；
ligands结构体向量中生成随机的position结构体和orientation四元数以及[-pi,pi]的随机向量torsions，
flex结构体向量中生成[-pi,pi]的随机向量torsions
*/
		manifold_phase(phase, tmp, m, out, p_widened, ig_widened, corner2corner, 1, e_internal_stats, *this, generator); // 1 - initial manifold
		
		
		if(use_ssd)
			ssd_par(m, p, ig, tmp, g, authentic_v);//（）重载 改变tmp的e参数  tmp.e=res=bfgs函数返回值f0
		else {
			quasi_newton quasi_newton_par; 
			quasi_newton_par.max_steps = ssd_par.evals;
			quasi_newton_par(m, p, ig, tmp, g, authentic_v);//（）重载 改变tmp的e参数  tmp.e=res=bfgs函数返回值f0
		}


		m.set(tmp.c);
/*set
* 根据配体配置c，获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
*/
		tmp.coords = m.get_heavy_atom_movable_coords();//给tmp的coords赋值
/*
* 函数：get_heavy_atom_movable_coords，常量函数
* function：1.建立一个临时vecv对象tmp
* 2.i从0遍历到m_num_movable_atoms，依次判断atoms第i个元素的成员el是否等于0，是则往tmp尾部添加model成员coords容器第i个元素
* 3.返回tmp
*/
		VINA_CHECK(tmp.coords.size() > 0);//判断coords是否为空
		std::pair<sz, fl> closest_rmsd = find_closest(tmp.coords, out);
		//复制find_closest
		

		//若 closest_rmsd第一个参数<out的元素个数且 第二个参数<min_rmsd   且tmp.e < out[closest_rmsd.first].e
		//则out[closest_rmsd.first] = tmp
		//否则在out的后面加入tmp
		if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
			if(tmp.e < out[closest_rmsd.first].e) { // the new one is better, apparently
				out[closest_rmsd.first] = tmp; // FIXME? slow
			}
		}
		else { // nothing similar
			out.push_back(new output_type(tmp));
		}
	}
	// final tunings
	const output_container null_array;//结构体向量
	VINA_CHECK(!out.empty());
	out.sort();//将out中结构体的参数按照升序排列
	VINA_CHECK(out.front().e <= out.back().e); // //验证排序   make sure the sorting worked in the correct order
}
