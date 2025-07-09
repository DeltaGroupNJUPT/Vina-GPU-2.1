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

#include <algorithm> // fill, etc

#if 0 // use binary cache ��ע�͵�
	// for some reason, binary archive gives four huge warnings in VC2008
	#include <boost/archive/binary_oarchive.hpp>
	#include <boost/archive/binary_iarchive.hpp>
	typedef boost::archive::binary_iarchive iarchive;
	typedef boost::archive::binary_oarchive oarchive;
#else // use text cache
	#include <boost/archive/text_oarchive.hpp>
	#include <boost/archive/text_iarchive.hpp>
	typedef boost::archive::text_iarchive iarchive;
	typedef boost::archive::text_oarchive oarchive;
#endif 

#include <boost/serialization/split_member.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/static_assert.hpp>
#include "cache.h"
#include "file.h"
#include "szv_grid.h"

//OpenCL Related
#include "wrapcl.h"
//#define DISPLAT_RUNTIME
//end OpenCL Related
/*
	scoring_function_version=scoring_function_version_��string��
	��gd=gd_����3��Ԫ�ص�grid_dim�ṹ�����飻
	��slope=slope_��double��
	��atu=atom_typing_used_��atom_type�ṹ����ö���ͱ���t��
	��grids=num_atom_types(atom_typing_used_)��
	  atom_typing_used_��atom_type�ṹ���е�ö���ͱ���t��ΪEL�ͷ���11��AD�ͷ���20��XS�ͷ���17��SY�ͷ���18 


*/
cache::cache(const std::string& scoring_function_version_, const grid_dims& gd_, fl slope_, atom_type::t atom_typing_used_) 
: scoring_function_version(scoring_function_version_), gd(gd_), slope(slope_), atu(atom_typing_used_), grids(num_atom_types(atom_typing_used_)) {}



/*
	input����m��model�ṹ�壬v��double
	output����e��double
	inter����atu��atom_type�ṹ���е�ö���ͱ���t ��
	��nat��unint������atuΪEL�ͷ���11��AD�ͷ���20��XS�ͷ���17��SY�ͷ���18 ��
	��t��unint������atuΪEL�ͷ���el��AD�ͷ���ad��XS�ͷ���xs��SY�ͷ���sy,��ʼ�� el=11; ad=20; xs=17; sy=18 ��
	�� m.atoms��model�ṹ���е�atom�ṹ��������m.coords��vec�ṹ������;
	�� m.num_movable_atoms()������model�ṹ���еģ�unint��m_num_movable_atoms����;
	��slope��double
	��������Ҫ���ܣ�
		��i=0��m_num_movable_atoms��
		if(t >= nat)   ��i++��ֱ�ӽ�����һ��ѭ����
		else 
			g=grid�ṹ��������t��Ԫ�ر�����
			m_data��ģ���β�Ϊdouble��arrat3d �ࣩ�Ĳ���sz m_i, m_j, m_kͬʱ>0�������������ֹ����          
			e=e+g.evaluate(m.coords[i], slope, v)��
	    ѭ����������e��
*/
fl cache::eval      (const model& m, fl v) const { // needs m.coords
	fl e = 0;
	sz nat = num_atom_types(atu);			//atu��atom_type�ṹ���е�ö���ͱ���t��ΪEL�ͷ���11��AD�ͷ���20��XS�ͷ���17��SY�ͷ���18 

	VINA_FOR(i, m.num_movable_atoms()) {	  //for(i=0;i<m_num_movable_atoms;i++)��
		const atom& a = m.atoms[i];           //atoms��atom�ṹ������
		sz t = a.get(atu);					//atu��atom_type�ṹ���е�ö���ͱ���t����ʼ��el=11; ad=20; xs=17; sy=18
		if(t >= nat) continue;              //��������ֱ�ӽ�����һ��ѭ�����������������
		const grid& g = grids[t];           //grid�ṹ��������t��Ԫ�ر���
		assert(g.initialized());            //��m_data��ģ���β�Ϊdouble��arrat3d �ࣩ�Ĳ���sz m_i, m_j, m_kͬʱ>0�򷵻�1  ����0
		e += g.evaluate(m.coords[i], slope, v);   //slope��double�����double���͵�����derivָ��Ϊnull 
	}
	return e;
}


/*
	input����m��model�ṹ�壬v��double
	output����e��double
	inter����atu��atom_type�ṹ���е�ö���ͱ���t ��
	��nat��unint������atuΪEL�ͷ���11��AD�ͷ���20��XS�ͷ���17��SY�ͷ���18 ��
	��t��unint������atuΪEL�ͷ���el��AD�ͷ���ad��XS�ͷ���xs��SY�ͷ���sy,��ʼ�� el=11; ad=20; xs=17; sy=18 ��
	�� m.atoms��model�ṹ���е�atom�ṹ��������m.coords�� m.minus_forces��vec�ṹ������;
	�� m.num_movable_atoms()������model�ṹ���еģ�unint��m_num_movable_atoms����;
	��slope��double
	��������Ҫ���ܣ�
	��i=0��m_num_movable_atoms��
	if(t >= nat)   ��i++���ҽṹ��minus_forces�е�data[0] = data[1] = data[2] = 0��Ȼ��ֱ�ӽ�����һ��ѭ����
	else
	g=grid�ṹ��������t��Ԫ�ر�����
	������deriv��vec�ṹ��
	m_data��ģ���β�Ϊdouble��arrat3d �ࣩ�Ĳ���sz m_i, m_j, m_kͬʱ>0�������������ֹ����
	e=e+g.evaluate(m.coords[i], slope, v��deriv)��
	m.minus_forces[i] = deriv;   //�ṹ�帳ֵ�ṹ��
	ѭ����������e��

*/
fl cache::eval_deriv(      model& m, fl v) const { // needs m.coords, sets m.minus_forces
	fl e = 0;
	sz nat = num_atom_types(atu);

	VINA_FOR(i, m.num_movable_atoms()) {
		const atom& a = m.atoms[i];
		sz t = a.get(atu);
		if(t >= nat) { m.minus_forces[i].assign(0); continue; }//�ṹ��minus_forces�е�data[0] = data[1] = data[2] = 0
		const grid& g = grids[t];
		assert(g.initialized());
		vec deriv;
		e += g.evaluate(m.coords[i], slope, v, deriv);
		m.minus_forces[i] = deriv;
	}
	return e;
}
/*
* Added by Glinttsd
* Function to get grids
*/
std::vector<grid> cache::get_grids() const{
	return this->grids;
}
/*
* Added by Glinttsd
* Function to get atu
*/
int cache::get_atu() const {
	return this->atu;
}
/*
* Added by Glinttsd
* Function to get atu
*/
double cache::get_slope() const {
	return this->slope;
}


#if 0 // No longer doing I/O of the cache  ��ע�͵�
void cache::read(const path& p) {
	ifile in(p, std::ios::binary);
	iarchive ar(in);
	ar >> *this;
}

void cache::write(const path& p) const {
	ofile out(p, std::ios::binary);
	oarchive ar(out);
	ar << *this;
}
#endif


/*
	���л��������л��Ĳ�����scoring_function_version(string)��gd(grid_dims�ṹ��)��
	atu��atom_type�ṹ���е�ö���ͱ���t ����grids��grid�ṹ��������
*/
template<class Archive>
void cache::save(Archive& ar, const unsigned version) const {
	ar & scoring_function_version;
	ar & gd;
	ar & atu;
	ar & grids;
}


/*
	�����л�
	��name_tmp��string��
	��gd_tmp��grid_dims�ṹ�壻
	��atu_tmp��atom_type�ṹ���е�ö���ͱ���t��
	��grids��grid�ṹ������
	������ܣ�load���ĺͷ��ͣ�save���Ĳ���ͬ�׳��쳣
*/
template<class Archive>
void cache::load(Archive& ar, const unsigned version) {
	std::string name_tmp;       ar & name_tmp;       if(name_tmp != scoring_function_version) throw energy_mismatch();  //�׳��쳣
	grid_dims   gd_tmp;         ar &   gd_tmp;       if(!eq(gd_tmp, gd))                    throw grid_dims_mismatch();
	atom_type::t atu_tmp;       ar &  atu_tmp;       if(atu_tmp != atu)                       throw cache_mismatch();

	ar & grids;
}


/*
	input����m:model�ṹ�壬p��precalculate�ṹ��
			 atom_types_needed��unint������display_progress��bool;
    inter����grids��grid�ṹ��������gd��grid_dims�ṹ��
*/

void cache::populate(const model& m, const precalculate& p, const szv& atom_types_needed, bool display_progress) {
	szv needed;                                     //unint����
	VINA_FOR_IN(i, atom_types_needed) {            //for(i=0;i<atom_types_needed.size;i++)
		sz t = atom_types_needed[i];               //��atom_types_needed�����е�Ԫ�����θ�ֵ��t
		if(!grids[t].initialized()) {              //grid�ṹ������άdouble����m_data�Ĳ���sz m_i, m_j, m_kͬʱ>0�򷵻�1  ����0
			needed.push_back(t);                   // needed����ĩβ����t
			grids[t].init(gd);                     //��grid�ṹ���е����vec�ṹ�帳ֵ
		}
	} 


	if(needed.empty())                             //���neededΪ�ռ�����grid�ṹ������άdouble����m_data�Ĳ���sz m_i, m_j, m_kͬʱ>0
		return;                                    //�������ؿ�ֵ 
												   //����needed��Ϊ��
	flv affinities(needed.size());                 //double��������needed����Ԫ��һ����

	sz nat = num_atom_types(atu);				   //atu��atom_type�ṹ���е�ö���ͱ���t��ΪEL�ͷ���11��AD�ͷ���20��XS�ͷ���17��SY�ͷ���18

	grid& g = grids[needed.front()];               //needed.front()������needed��һ��Ԫ�ص�����

	const fl cutoff_sqr = p.cutoff_sqr();          //����precalculate�ṹ���е�m_cutoff_sqr

	grid_dims gd_reduced = szv_grid_dims(gd);      //gd����3��Ԫ�ص�grid_dim�ṹ�����飬����grid_dims��tmp[i]�Ĳ���begin   end    n��ֵ  Ȼ�󷵻�tmp
	szv_grid ig(m, gd_reduced, cutoff_sqr);        //��index[i]��ֵ    ����   m_data[index[0] + m_i*(index[1] + m_j*index[2])]      m_data��array3d<szv>����       

	VINA_FOR(x, g.m_data.dim0()) {                  //for(x=0;x<m_i;x++)��m_i��array3d���еĳ�Ա����
		VINA_FOR(y, g.m_data.dim1()) {              //for(y=0;y<m_j;y++)��m_j��array3d���еĳ�Ա����
			VINA_FOR(z, g.m_data.dim2()) {          //for(z=0;z<m_k;z++)��m_k��array3d���еĳ�Ա����
				std::fill(affinities.begin(), affinities.end(), 0); //��affinities�������0
				vec probe_coords;                  //����һ��vec�ṹ��
				probe_coords = g.index_to_argument(x, y, z);
				/*m_init��m_factor_inv��grid�ṹ���е�vec�ṹ�壬�䶼����������
				vec��m_init[0] + m_factor_inv[0] * x,
					 m_init[1] + m_factor_inv[1] * y,
					 m_init[2] + m_factor_inv[2] * z ��
				*/
				const szv& possibilities = ig.possibilities(probe_coords);	//��index[i]��ֵ    ����   m_data[index[0] + m_i*(index[1] + m_j*index[2])]��index�� m_data��array3d<szv>����
				VINA_FOR_IN(possibilities_i, possibilities) {				//for(possibilities_i=0;possibilities_i<possibilities.size;possibilities_i++)
					const sz i = possibilities[possibilities_i];			//��possibilities����Ԫ�����θ�ֵ��i
					const atom& a = m.grid_atoms[i];						//grid_atoms������atom�ṹ����������Ϊa
					const sz t1 = a.get(atu);								//atu��atom_type�ṹ���е�ö���ͱ���t����ʼ��el=11; ad=20; xs=17; sy=18
					
					if(t1 >= nat) continue;									// ��������ֱ�ӽ�����һ��ѭ�����������������
					const fl r2 = vec_distance_sqr(a.coords, probe_coords);
					/*
						demo:vec a.coords(4,4,4)
						     vec probe_coords(1,2,3)
						  r2 = vec_distance_sqr(a,b)
						//r2 = (4-1)^2+(4-2)^2+(4-3)^2=14 
					*/

					if(r2 <= cutoff_sqr) {                                 
						VINA_FOR_IN(j, needed) {                           //for(j=0;j<needed.size;j++)
							const sz t2 = needed[j];                       //��needed����Ԫ�����θ�ֵ��t2
							assert(t2 < nat);
							const sz type_pair_index = triangular_matrix_index_permissive(num_atom_types(atu), t1, t2);
							/*atu��EL�ͷ���11;AD�ͷ���20;XS�ͷ���17;SY�ͷ���18
							��t2 < num_atom_types(atu)��t1< num_atom_types(atu)����£�
		                    ��t1 <= t2������t1 + t2*(t2+1)/2��
		                    ���򷵻�t2 + t1*(t1+1)/2	
							*/
							affinities[j] += p.eval_fast(type_pair_index, r2);
						}
					}           //if(r2 <= cutoff_sqr)����
				}               //for(possibilities_i=0;possibilities_i<possibilities.size;possibilities_i++)����
				VINA_FOR_IN(j, needed) {						//for(j=0;j<needed.size;j++)
					sz t = needed[j];							//��needed����Ԫ�����θ�ֵ��t
					assert(t < nat);
					grids[t].m_data(x, y, z) = affinities[j];   // m_data[x+ m_i*(y + m_j*z)]=affinities[j]
				}
			}               //for(z=0;z<m_k;z++)����
		}					//for(y=0;y<m_j;y++)����
	}						//for(x=0;x<m_i;x++)����
}

bool cache::m_data_check(const szv_grid ig)const {
	//Check if m_data is too large
	size_t max_m_data_size = 0;
	for (int i = 0; i < ig.m_data.m_data.size(); i++) {
		size_t m_data_size = ig.m_data.m_data[i].size();
		if (m_data_size > max_m_data_size)max_m_data_size = m_data_size;
	}
	return(		ig.m_data.dim0() <= MAX_M_DATA_MI
			||	ig.m_data.dim1() <= MAX_M_DATA_MJ
			||	ig.m_data.dim2() <= MAX_M_DATA_MK
			||	max_m_data_size <= MAX_NUM_OF_EVERY_M_DATA_ELEMENT);
	//end Check
}
#ifdef OPENCL_PART_1
void cache::populate_cl(const model& m, const precalculate& p, const szv& atom_types_needed, bool display_progress) {
/**************************************************************************/
/***************************    OpenCL Init    ****************************/
/**************************************************************************/

	cl_int err;
	cl_platform_id* platforms;
	cl_device_id* devices;
	cl_context context;
	cl_command_queue queue;
	cl_int gpu_platform_id = 0;
	SetupPlatform(&platforms, &gpu_platform_id);
	SetupDevice(platforms, &devices, gpu_platform_id);
	SetupContext(platforms, devices, &context, 1, gpu_platform_id);
	SetupQueue(&queue, context, devices);
	char* program_file;
	cl_program program_cl;
	cl_program program;
	size_t program_size;
	//Read kernel source code
	/*
	const std::string default_work_path = "D:/VScode_Project/AutoDock_vina_Opencl_float";
	const std::string include_path = default_work_path + "/OpenCL/inc"; //FIX it
	const std::string addtion = "";

	//read_file(&program_file, &program_size, default_work_path + "/OpenCL/src/kernels/kernel2.cl");
	//program_cl = clCreateProgramWithSource(context, 1, (const char**)&program_file, &program_size, &err); checkErr(err);

	char* program_file_n[1];
	size_t program_size_n[1];
	std::string file_paths[1] = {default_work_path + "/OpenCL/src/kernels/kernel1.cl" }; // The order of files is important!
	read_n_file(program_file_n, program_size_n, file_paths, 1);
	std::string final_file;
	size_t final_size = 1 - 1; // count '\n'
	for (int i = 0; i < 1; i++) {
		if (i == 0) final_file = program_file_n[0];
		else final_file = final_file + '\n' + (std::string)program_file_n[i];
		final_size += program_size_n[i];
	}
	const char* final_files_char = final_file.data();

	program_cl = clCreateProgramWithSource(context, 1, (const char**)&final_files_char, &final_size, &err); checkErr(err);
	SetupBuildProgramWithSource(program_cl, NULL, devices, include_path, addtion);

	SaveProgramToBinary(program_cl, "Kernel1_Opt.bin");
	*/

	program_cl = SetupBuildProgramWithBinary(context, devices, "Kernel1_Opt.bin");

	err = clUnloadPlatformCompiler(platforms[gpu_platform_id]); checkErr(err);
	//Set kernel arguments
	cl_kernel kernels[1];
	char kernel_name[][50] = { "kernel1" };
	SetupKernel(kernels, program_cl, 1, kernel_name);

/**************************************************************************/
/************************    Original Vina code    ************************/
/**************************************************************************/

	szv needed;                                     //unint����
	VINA_FOR_IN(i, atom_types_needed) {            //for(i=0;i<atom_types_needed.size;i++)
		sz t = atom_types_needed[i];               //��atom_types_needed�����е�Ԫ�����θ�ֵ��t
		if (!grids[t].initialized()) {              //grid�ṹ������άdouble����m_data�Ĳ���sz m_i, m_j, m_kͬʱ>0�򷵻�1  ����0
			needed.push_back(t);                   // needed����ĩβ����t
			grids[t].init(gd);                     //��grid�ṹ���е����vec�ṹ�帳ֵ
		}
	}

	if (needed.empty())                             //���neededΪ�ռ�����grid�ṹ������άdouble����m_data�Ĳ���sz m_i, m_j, m_kͬʱ>0
		return;                                    //�������ؿ�ֵ 
												   //����needed��Ϊ��
	flv affinities(needed.size());                 //double��������needed����Ԫ��һ����

	sz nat = num_atom_types(atu);				   //atu��atom_type�ṹ���е�ö���ͱ���t��ΪEL�ͷ���11��AD�ͷ���20��XS�ͷ���17��SY�ͷ���18

	grid& g = grids[needed.front()];               //needed.front()������needed��һ��Ԫ�ص�����

	const fl cutoff_sqr = p.cutoff_sqr();          //����precalculate�ṹ���е�m_cutoff_sqr

	grid_dims gd_reduced = szv_grid_dims(gd);      //gd����3��Ԫ�ص�grid_dim�ṹ�����飬����grid_dims��tmp[i]�Ĳ���begin   end    n��ֵ  Ȼ�󷵻�tmp
	szv_grid ig(m, gd_reduced, cutoff_sqr);        //��index[i]��ֵ    ����   m_data[index[0] + m_i*(index[1] + m_j*index[2])]      m_data��array3d<szv>����       

	//Check if ig.m_init and ig.m_range equal to g.m_init and g.m_range, 
	//if so, m_init and m_range can be used in kernel, otherwise, ig.m_init and ig.m_range need to be allocated 
	///���ｫ�ṹ��szv_grid�еĳ�Ա��Ϊ��public������ֱ�ӷ���
	for (int i = 0; i < 3; i++) {
		if (ig.m_init[i] != g.m_init[i]) {
			printf("m_init not equal!");
			exit(-1);
		}
		if (ig.m_range[i] != g.m_range[i]) {
			printf("m_range not equal!");
			exit(-1);
		}
	}

/**************************************************************************/
/************************    Allocate CPU memory    ***********************/
/**************************************************************************/
	grid_cl g_cl;
	for (int i = 0; i < 3; i++) {
		g_cl.m_init[i]				= g.m_init[i];
		g_cl.m_factor[i]			= g.m_factor[i];
		g_cl.m_dim_fl_minus_1[i]	= g.m_dim_fl_minus_1[i];
		g_cl.m_factor_inv[i]		= g.m_factor_inv[i];
	}
	if (g.m_data.dim0() != 0) {
		g_cl.m_i = g.m_data.dim0(); assert(GRID_MI == g_cl.m_i);
		g_cl.m_j = g.m_data.dim1(); assert(GRID_MJ == g_cl.m_j);
		g_cl.m_k = g.m_data.dim2(); assert(GRID_MK == g_cl.m_k);
		for (int j = 0; j < GRID_MI * GRID_MJ * GRID_MK; j++) {
			g_cl.m_data[j] = g.m_data.m_data[j];
		}
		//memcpy(g_cl.m_data, g.m_data.m_data.data(), GRID_MI * GRID_MJ * GRID_MK * sizeof(float));
	}
	else {
		g_cl.m_i = 0;
		g_cl.m_j = 0;
		g_cl.m_k = 0;
	}

	grid_atoms_cl ga_cl;
	for (int i = 0; i < m.grid_atoms.size(); i++) {
		ga_cl.atoms[i].types[0] = m.grid_atoms[i].el;
		ga_cl.atoms[i].types[1] = m.grid_atoms[i].ad;
		ga_cl.atoms[i].types[2] = m.grid_atoms[i].xs;
		ga_cl.atoms[i].types[3] = m.grid_atoms[i].sy;
		for (int j = 0; j < 3; j++)ga_cl.atoms[i].coords[j] = m.grid_atoms[i].coords.data[j];
	}

	// Preaparing p related data
	p_cl p_cl;
	p_cl.m_cutoff_sqr = p.cutoff_sqr();
	p_cl.factor = p.factor;
	p_cl.n = p.n;
	assert(MAX_P_DATA_M_DATA_SIZE > p.data.m_data.size());
	for (int i = 0; i < p.data.m_data.size(); i++) {
		p_cl.m_data[i].factor = p.data.m_data[i].factor;
		assert(FAST_SIZE == p.data.m_data[i].fast.size());
		assert(SMOOTH_SIZE == p.data.m_data[i].smooth.size());

		for (int j = 0; j < FAST_SIZE; j++) {
			p_cl.m_data[i].fast[j] = p.data.m_data[i].fast[j];
		}
		for (int j = 0; j < SMOOTH_SIZE; j++) {
			p_cl.m_data[i].smooth[j][0] = p.data.m_data[i].smooth[j].first;
			p_cl.m_data[i].smooth[j][1] = p.data.m_data[i].smooth[j].second;
		}

		//memcpy(p_cl.m_data[i].fast, p.data.m_data[i].fast.data(), FAST_SIZE * sizeof(double));
		//memcpy(p_cl.m_data[i].smooth, p.data.m_data[i].smooth.data(), 2 * SMOOTH_SIZE * sizeof(double));
	}

	// Preaparing ig related data
	int m_data_dims[3] = { ig.m_data.dim0(),ig.m_data.dim1(),ig.m_data.dim2() };
	float ig_m_init[3] = { ig.m_init.data[0],ig.m_init.data[1] ,ig.m_init.data[2] };
	float ig_m_range[3] = { ig.m_range.data[0],ig.m_range.data[1] ,ig.m_range.data[2] };
	// Others
	int needed_size = needed.size();
	float epsilon_fl_float = epsilon_fl;

/**************************************************************************/
/************************    Allocate GPU memory    ***********************/
/**************************************************************************/
	// Preparing g related data
	cl_mem g_cl_gpu;
	CreateDeviceBuffer(&g_cl_gpu, CL_MEM_READ_WRITE, sizeof(g_cl), context);
	err = clEnqueueWriteBuffer(queue, g_cl_gpu, false, 0, sizeof(g_cl), &g_cl, 0, NULL, NULL); checkErr(err);

	// Preparing ga_cl related data
	cl_mem ga_cl_gpu;
	CreateDeviceBuffer(&ga_cl_gpu, CL_MEM_READ_ONLY, (size_t)sizeof(ga_cl), context);
	err = clEnqueueWriteBuffer(queue, ga_cl_gpu, false, 0, (size_t)sizeof(ga_cl), &ga_cl, 0, NULL, NULL); checkErr(err);

	// Preparing p related data
	cl_mem p_cl_gpu;
	CreateDeviceBuffer(&p_cl_gpu, CL_MEM_READ_ONLY, (size_t)sizeof(p_cl), context);
	err = clEnqueueWriteBuffer(queue, p_cl_gpu, false, 0, (size_t)sizeof(p_cl), &p_cl, 0, NULL, NULL); checkErr(err);

	// Preparing affinities data
	cl_mem results;
	CreateDeviceBuffer(&results, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, (GRID_MI * GRID_MJ * GRID_MK) * sizeof(affinities_cl), context);
	// ���0
	float fill_data = 0;
	//clEnqueueFillBuffer(queue, results, false, 0, (GRID_MI* GRID_MJ* GRID_MK) * sizeof(affinities_cl), &fill_data, 0, NULL, NULL); checkErr(err);
	err = clEnqueueFillBuffer(queue, results, &fill_data, sizeof(float), 0, (GRID_MI* GRID_MJ* GRID_MK) * sizeof(affinities_cl), 0, NULL, NULL); checkErr(err);
	clFinish(queue);

	//Preparing ig related data
	cl_mem ig_m_data_gpu;
	if (m_data_check(ig) == true) {
		CreateDeviceBuffer(&ig_m_data_gpu, CL_MEM_READ_ONLY, MAX_NUM_OF_TOTAL_M_DATA * sizeof(int), context);
		//total_mem += MAX_NUM_OF_TOTAL_M_DATA * sizeof(int);
		int fill_zero = 0;
		//printf("\n infinity = %0.16f", fill_inf);
		err = clEnqueueFillBuffer(	queue, ig_m_data_gpu, &fill_zero, sizeof(fill_zero), 0,
									MAX_NUM_OF_TOTAL_M_DATA * sizeof(int), 0, NULL, NULL); checkErr(err);
		for (int i = 0; i < ig.m_data.m_data.size(); i++) {
			if (ig.m_data.m_data[i].size() == 0)continue;
			err = clEnqueueWriteBuffer(	queue, ig_m_data_gpu, false, i * MAX_NUM_OF_EVERY_M_DATA_ELEMENT * sizeof(int), 
										ig.m_data.m_data[i].size() * sizeof(int), ig.m_data.m_data[i].data(), 0, NULL, NULL); checkErr(err);
		}
	}

	cl_mem ig_m_init_gpu;
	CreateDeviceBuffer(&ig_m_init_gpu, CL_MEM_READ_WRITE, 3 * sizeof(float), context);
	err = clEnqueueWriteBuffer(queue, ig_m_init_gpu, false, 0, 3 * sizeof(float), ig_m_init, 0, NULL, NULL); checkErr(err);

	cl_mem ig_m_range_gpu;
	CreateDeviceBuffer(&ig_m_range_gpu, CL_MEM_READ_WRITE, 3 * sizeof(float), context);
	err = clEnqueueWriteBuffer(queue, ig_m_range_gpu, false, 0, 3 * sizeof(float), ig_m_range, 0, NULL, NULL); checkErr(err);

	cl_mem m_data_dims_gpu;
	CreateDeviceBuffer(&m_data_dims_gpu, CL_MEM_READ_ONLY, 3 * sizeof(int), context);
	err = clEnqueueWriteBuffer(queue, m_data_dims_gpu, false, 0, 3 * sizeof(int), m_data_dims, 0, NULL, NULL); checkErr(err);

	// Others
	cl_mem needed_gpu;
	CreateDeviceBuffer(&needed_gpu, CL_MEM_READ_ONLY, needed_size * sizeof(int), context);
	err = clEnqueueWriteBuffer(queue, needed_gpu, false, 0, needed_size * sizeof(int), needed.data(), 0, NULL, NULL); checkErr(err);


/**************************************************************************/
/************************   Set kernel arguments    ***********************/
/**************************************************************************/
	//int atu_cl = 2;
	SetKernelArg(kernels[0], 0, sizeof(cl_mem), &g_cl_gpu);
	SetKernelArg(kernels[0], 1, sizeof(cl_mem), &ga_cl_gpu);
	SetKernelArg(kernels[0], 2, sizeof(cl_mem), &p_cl_gpu);
	SetKernelArg(kernels[0], 3, sizeof(int),	&atu);
	SetKernelArg(kernels[0], 4, sizeof(int),	&nat);
	SetKernelArg(kernels[0], 5, sizeof(float),  &epsilon_fl_float);
	SetKernelArg(kernels[0], 6, sizeof(cl_mem), &ig_m_data_gpu);
	SetKernelArg(kernels[0], 7, sizeof(cl_mem), &ig_m_init_gpu);
	SetKernelArg(kernels[0], 8, sizeof(cl_mem), &ig_m_range_gpu);
	SetKernelArg(kernels[0], 9, sizeof(cl_mem), &m_data_dims_gpu);
	SetKernelArg(kernels[0], 10, sizeof(int),	&needed_size);
	SetKernelArg(kernels[0], 11, sizeof(int),	&needed_gpu);
	SetKernelArg(kernels[0], 12, sizeof(cl_mem), &results);

/**************************************************************************/
/****************************   Start kernel    ***************************/
/**************************************************************************/
	size_t global_size[3] = { 128, 128, 64 };
	size_t local_size[3] = { 4,4,4 };

	cl_event cache_cl;
	err = clEnqueueNDRangeKernel(queue, kernels[0], 3, 0, global_size, local_size, 0, NULL, &cache_cl); checkErr(err);

	clWaitForEvents(1, &cache_cl);

	affinities_cl* result_ptr = (affinities_cl*)clEnqueueMapBuffer(queue, results, CL_TRUE, CL_MAP_READ,
																	0, (GRID_MI * GRID_MJ * GRID_MK) * sizeof(affinities_cl),
																	0, NULL, NULL, &err); checkErr(err);

	VINA_FOR(x, g.m_data.dim0()) {                  
		VINA_FOR(y, g.m_data.dim1()) {              
			VINA_FOR(z, g.m_data.dim2()) {          
				VINA_FOR_IN(j, needed) {			
					sz t = needed[j];				
					assert(t < nat);
					grids[t].m_data(x, y, z) = result_ptr[z * (GRID_MJ * GRID_MI) + y * (GRID_MI) + x].data[j];   // m_data[x+ m_i*(y + m_j*z)]=affinities[j]
				}
			}
		}
	}

	err = clEnqueueUnmapMemObject(queue, results, result_ptr, 0, NULL, NULL); checkErr(err);

	// Memory objects release
	err = clReleaseMemObject(g_cl_gpu); checkErr(err);
	err = clReleaseMemObject(ga_cl_gpu); checkErr(err);
	err = clReleaseMemObject(p_cl_gpu); checkErr(err);
	err = clReleaseMemObject(results); checkErr(err);

#ifdef DISPLAT_RUNTIME
	// Output Analysis
	cl_ulong time_start, time_end;
	double total_time;
	err = clGetEventProfilingInfo(cache_cl, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL); checkErr(err);
	err = clGetEventProfilingInfo(cache_cl, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL); checkErr(err);
	total_time = time_end - time_start;
	printf("\n GPU cache runtime = %0.16f s\n", (total_time / 1000000000.0));
#endif
}
#endif