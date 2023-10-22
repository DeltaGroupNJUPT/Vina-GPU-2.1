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

#ifndef VINA_TREE_H
#define VINA_TREE_H

#include "conf.h"
#include "atom.h"



/*
	orientation_m——mat结构体；orientation_q——四元数
	origin=origin_——vec结构体别名，用来存放四元数三个虚部；
	初始化：
	orientation_q——四元数，(1, 0, 0, 0)；
	orientation_m——mat结构体（对角矩阵）
	data[0]=1 ；data[3]=0 ；data[6]=0 ;
	data[1]=0 ；data[4]=1 ；data[7]=0 ;
	data[2]=0 ；data[5]=0 ；data[8]=1 ;
*/

struct frame {
	frame(const vec& origin_) : origin(origin_), orientation_q(qt_identity), orientation_m(quaternion_to_r3(qt_identity)) {}
	/*
		input——local_coords：vec结构体;
		output——tmp：vec结构体;
		demo：
		origin（data0[0],data0[1],data0[2]）
		orientation_m（data1[0],data1[3],data1[6]
				       data1[1],data1[4],data1[7]
					   data1[2],data1[5],data1[8]）
		local_coords（data2[0],data2[1],data2[2]）

		tmp(data0[0]+data1[0]*data2[0] + data1[3]*data2[1] + data1[6]*data2[2],
		    data0[1]+data1[1]*data2[0] + data1[4]*data2[1] + data1[7]*data2[2],
		    data0[2]+data1[2]*data2[0] + data1[5]*data2[1] + data1[8]*data2[2])
	*/
	vec local_to_lab(const vec& local_coords) const {
		vec tmp;
		tmp = origin + orientation_m*local_coords; 
		return tmp;
	}

	/*
		input——local_direction：vec结构体;
		output——tmp：vec结构体;
		demo:
		orientation_m（data1[0],data1[3],data1[6]
					   data1[1],data1[4],data1[7]
					   data1[2],data1[5],data1[8]）
		local_direction（data2[0],data2[1],data2[2]）

		tmp(data1[0]*data2[0] + data1[3]*data2[1] + data1[6]*data2[2],
			data1[1]*data2[0] + data1[4]*data2[1] + data1[7]*data2[2],
			data1[2]*data2[0] + data1[5]*data2[1] + data1[8]*data2[2])
	
	*/
	vec local_to_lab_direction(const vec& local_direction) const {
		vec tmp;
		tmp = orientation_m * local_direction;
		return tmp;
	}


	const qt& orientation() const { return orientation_q; }//返回orientation_q四元数
	const vec& get_origin() const { return origin; }       //返回origin结构体
	const mat& get_orientation_m() const { return orientation_m; }
protected:
	vec origin;
	/*	设置四元数orientation_q=q和orientation_m结构体
		input——q：四元数
	*/
	void set_orientation(const qt& q) { // does not normalize the orientation
		orientation_q = q;
		orientation_m = quaternion_to_r3(orientation_q);
	}
	mat orientation_m;
	qt  orientation_q;
};

/*
	begin、end：unint
*/
struct atom_range {
    sz begin;
    sz end;
	atom_range(sz begin_, sz end_) : begin(begin_), end(end_) {}
	template<typename F>
	/*	转换
		input——f：任何类型
		diff：unint
	*/
	void transform(const F& f) {
		sz diff = end - begin;
		begin = f(begin);
		end   = begin + diff;
	}
};


/*	原子结构，继承atom_range、frame结构体

    frame结构体中：
	origin=origin_——vec结构体别名，用来存放四元数三个虚部；
	初始化：
	orientation_q——四元数，(1, 0, 0, 0)；
	orientation_m——mat结构体（对角矩阵）
	data[0]=1 ；data[3]=0 ；data[6]=0 ;
	data[1]=0 ；data[4]=1 ；data[7]=0 ;
	data[2]=0 ；data[5]=0 ；data[8]=1 ;
	
	atom_range结构体中：
	begin=begin_；end=end_；

*/
struct atom_frame : public frame, public atom_range {
	atom_frame(const vec& origin_, sz begin_, sz end_) : frame(origin_), atom_range(begin_, end_) {}
	/*
		input——atoms：atom结构体向量；coords：vec结构体向量
		demo：
		origin（data0[0],data0[1],data0[2]）
		orientation_m（data1[0],data1[3],data1[6]
					   data1[1],data1[4],data1[7]
					   data1[2],data1[5],data1[8]）
		coords （data2[0],data2[1],data2[2]）

		tmp(data0[0]+data1[0]*data2[0] + data1[3]*data2[1] + data1[6]*data2[2],
		    data0[1]+data1[1]*data2[0] + data1[4]*data2[1] + data1[7]*data2[2],
		    data0[2]+data1[2]*data2[0] + data1[5]*data2[1] + data1[8]*data2[2])
			这是coords结构体一个元素的内容，有（begin—end）个
	*/
	void set_coords(const atomv& atoms, vecv& coords) const {
		VINA_RANGE(i, begin, end)
			coords[i] = local_to_lab(atoms[i].coords);
	}

	/*
		input——coords、forces：vec结构体向量；
		output——tmp：<vec,vec>结构体向量；
		先初始化tmp结构体向量<(0,0,0),(0,0,0)>，tmp有两个成员变量first、second
		demo：
			forces{（0，0，0），（1，1，1），（2，2，2）}
			coords{（1，1，1），（1，1，1），（2，2，2）}
			origin（1,1,1）
			coords[i] - origin={(1-1,1-1,1-1),(1-1,1-1,1-1),(2-1,2-1,2-1)}={（0，0，0），（0，0，0），（1，1，1）}
			tmp<(0+0+1+2=3,0+0+1+2=3,0+0+1+2=3),(0+0+0+0,0+0+0+0,0+0+0+0)>=<(3,3,3),(0,0,0)>
	
	*/
	vecp sum_force_and_torque(const vecv& coords, const vecv& forces) const {
		vecp tmp;
		tmp.first.assign(0);
		tmp.second.assign(0);
		VINA_RANGE(i, begin, end) {
			tmp.first  += forces[i]; 
			tmp.second += cross_product(coords[i] - origin, forces[i]);
		}
		return tmp;
	}
};

/*  刚体，继承atom_frame结构体
	origin_：vec结构体，begin_、end_：unint；
	frame结构体中的origin=origin_；atom_range结构体中begin=begin_；end=end_
*/
struct rigid_body : public atom_frame {
	rigid_body(const vec& origin_, sz begin_, sz end_) : atom_frame(origin_, begin_, end_) {}

	/*配置设置
	input——atoms：atom结构体向量；coords：vec结构体向量；c：rigid_conf结构体
		·atom_frame.frame.origin=rigid_conf.position
		·设置四元数orientation_q=rigid_conf.orientation和orientation_m结构体
		·设置coords结构体向量
		     demo：
		     origin（data0[0],data0[1],data0[2]）
		     orientation_m（data1[0],data1[3],data1[6]
		     			   data1[1],data1[4],data1[7]
		     			   data1[2],data1[5],data1[8]）
		     coords （data2[0],data2[1],data2[2]）
		     
		     tmp(data0[0]+data1[0]*data2[0] + data1[3]*data2[1] + data1[6]*data2[2],
		    data0[1]+data1[1]*data2[0] + data1[4]*data2[1] + data1[7]*data2[2],
		    data0[2]+data1[2]*data2[0] + data1[5]*data2[1] + data1[8]*data2[2])
			这是coords结构体一个元素的内容，有（begin—end）个
	*/
	void set_conf(const atomv& atoms, vecv& coords, const rigid_conf& c) {
		origin = c.position;
		set_orientation(c.orientation);
		set_coords(atoms, coords);
	}

	void count_torsions(sz& s) const {} // do nothing这个函数没有什么作用

	/*设置派生物
		input——force_torque：<vec, vec>结构体向量有两个成员函数first、second，c:rigid_change结构体
	*/
	void set_derivative(const vecp& force_torque, rigid_change& c) const {
		c.position     = force_torque.first;
		c.orientation  = force_torque.second;
	}
};
/*坐标结构，继承atom_frame结构体
	origin_、axis_root、diff：vec结构体，begin_、end_：unint，nrm：double；
	·frame结构体中的origin=origin_；atom_range结构体中begin=begin_；end=end_
	·diff=origin结构体减去axis_root结构体（“-”重构）
	·nrm=diff结构体中data数组各元素平方和，开根号(空间向量模)
	·axis：vec 结构体
	demo：
		diff（2，3，4）；1/nrm=2;
		 axis=（4,6,8)
*/
struct axis_frame : public atom_frame {
	axis_frame(const vec& origin_, sz begin_, sz end_, const vec& axis_root) : atom_frame(origin_, begin_, end_) {
		vec diff; diff = origin - axis_root;
		fl nrm = diff.norm();
		VINA_CHECK(nrm >= epsilon_fl);     //要大于2.22045e-16
		axis = (1/nrm) * diff;
	}

	/* 设置派生物
		input——force_torque：<vec,vec>结构体向量；c：double
		axis：vec 结构体
		demo：
		force_torque.second（1，2，3），axis（2，3，4）
		c=1*2+2*3+3*4=20
	*/
	void set_derivative(const vecp& force_torque, fl& c) const {
		c = force_torque.second * axis;
	}
	vec get_axis() { return axis; } // 返回axis

protected:
	vec axis;
};
/*   继承axis_frame结构体
	origin_、axis_root、relative_axis、relative_origin：vec结构体，begin_、end_：unint，parent：frame结构体；
	·frame结构体中的origin=origin_；atom_range结构体中begin=begin_；end=end_
	·diff=origin结构体减去axis_root结构体（“-”重构）
	·axis_frame结构体中nrm=diff结构体中data数组各元素平方和，开根号(空间向量模)
	·axis_frame结构体中axis：vec 结构体
		demo：
			diff（2，3，4）；1/nrm=2;
			 axis=（4,6,8)
	·判断frame结构体中orientation_q是否等于(1, 0, 0, 0)
	·relative_axis=axis
	·relative_origin：
	  demo：origin（4，5，6），frame结构体中origin（1，1，1）
	  relative_origin=（3，4，5）
*/
struct segment : public axis_frame {
	segment(const vec& origin_, sz begin_, sz end_, const vec& axis_root, const frame& parent) : axis_frame(origin_, begin_, end_, axis_root) {
		VINA_CHECK(eq(parent.orientation(), qt_identity)); // the only initial parent orientation this c'tor supports
		relative_axis = axis;
		relative_origin = origin - parent.get_origin();
	}

	/*获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
		input——parent：frame结构体，atoms、coords：atom结构体向量，c：迭代器
		torsion：double
		tmp：角度四元数*frame结构体中orientation_q四元数
	
	*/
	void set_conf(const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
		const fl torsion = *c;  //指向的元素的引用
		++c;                    //指向下一个元素
		origin = parent.local_to_lab(relative_origin);          //获得坐标
		axis = parent.local_to_lab_direction(relative_axis);    //获得起点
		qt tmp = angle_to_quaternion(axis, torsion) * parent.orientation();//角度四元数*frame结构体中orientation_q四元数
		quaternion_normalize_approx(tmp); // normalization added in 1.1.2  判断四元数是否近似标准化
		//quaternion_normalize(tmp); // normalization added in 1.1.2
		set_orientation(tmp);
		set_coords(atoms, coords);
	}
	/*扭矩计数
	 input——s:unint
	*/
	void count_torsions(sz& s) const {
		++s;
	}
public:
	vec relative_axis;
	vec relative_origin;
};


/*继承axis_frame结构体

*/
struct first_segment : public axis_frame {
	first_segment(const segment& s) : axis_frame(s) {}
	first_segment(const vec& origin_, sz begin_, sz end_, const vec& axis_root) : axis_frame(origin_, begin_, end_, axis_root) {}
	/*设置角度四元数orientation_q=q和orientation_m结构体
		input——atoms：atom结构体，coords：vec结构体，torsion：double
	*/
	void set_conf(const atomv& atoms, vecv& coords, fl torsion) {
		set_orientation(angle_to_quaternion(axis, torsion));//角度四元数
  	}

	/*扭矩计数
	input——s:unint
	*/
	void count_torsions(sz& s) const {
		++s;
	}
};

/*
 获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
*/
template<typename T> // T == branch
void branches_set_conf(std::vector<T>& b, const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
	VINA_FOR_IN(i, b)                            //for(i=0;i<v.size;i++)
		b[i].set_conf(parent, atoms, coords, c);
}


/*需要根据b的类型调用函数
	input——b：任意向量，origin：vec结构体，forces、out、coords：vec结构体向量
	d：迭代器
*/
template<typename T> // T == branch
void branches_derivative(const std::vector<T>& b, const vec& origin, const vecv& coords, const vecv& forces, vecp& out, flv::iterator& d) { // adds to out
	VINA_FOR_IN(i, b) {											//for(i=0;i<v.size;i++)
		vecp force_torque = b[i].derivative(coords, forces, d);
		out.first  += force_torque.first;
		vec r; r = b[i].node.get_origin() - origin;
		out.second += cross_product(r, force_torque.first) + force_torque.second;
	}
}


/*
    node：任意类型变量 ，children：tree结构体向量
	node=node_
*/
template<typename T> // T == segment
struct tree {
	T node;
	std::vector< tree<T> > children;
	tree(const T& node_) : node(node_) {}
	/*	配置，根据node类型去调用相对于的函数
		input——parent：frame结构体，atoms：atom结构体向量，coords：vec结构体向量
		c：迭代器
	*/
	void set_conf(const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
		node.set_conf(parent, atoms, coords, c);
		branches_set_conf(children, node, atoms, coords, c);
	}

	/*	根据node类型去调用相对于的函数
		input——coords、forces：vec结构体向量，p：迭代器
		output——force_torque：<vec, vec>结构体向量
	*/
	vecp derivative(const vecv& coords, const vecv& forces, flv::iterator& p) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		fl& d = *p; // reference   取当前迭代器指向的元素
		++p;        //迭代器指向下一个元素
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, d);
		return force_torque;
	}
};

typedef tree<segment> branch;       //形参为segment结构体的tree结构体
typedef std::vector<branch> branches;  //branch结构体向量

template<typename Node> // Node == first_segment || rigid_body
struct heterotree {
	Node node;          //任意类型变量
	branches children;  //branch结构体向量
	heterotree(const Node& node_) : node(node_) {}
	/*
	    获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
		input——atoms：atom结构体向量，coords：vec结构体向量，c：ligand_conf结构体
	*/
	void set_conf(const atomv& atoms, vecv& coords, const ligand_conf& c) {
		node.set_conf(atoms, coords, c.rigid);
		flv::const_iterator p = c.torsions.begin();   //p迭代器指向torsions向量的第一个元素
		branches_set_conf(children, node, atoms, coords, p);
		assert(p == c.torsions.end());               //p迭代器指向torsions向量的最后一个元素
	}

	/* 获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
	 input——atoms：atom结构体向量，coords：vec结构体向量，c：residue_conf结构体
	*/
	void set_conf(const atomv& atoms, vecv& coords, const residue_conf& c) {
		flv::const_iterator p = c.torsions.begin();//p迭代器指向torsions向量的第一个元素
		node.set_conf(atoms, coords, *p);
		++p;
		branches_set_conf(children, node, atoms, coords, p);
		assert(p == c.torsions.end());               //p迭代器指向torsions向量的最后一个元素
	}

	/*
	input——coords、forces：vec结构体向量，c：ligand_change结构体
	*/
	void derivative(const vecv& coords, const vecv& forces, ligand_change& c) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		flv::iterator p = c.torsions.begin();             //p迭代器指向torsions向量的第一个元素
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, c.rigid);
		assert(p == c.torsions.end());                    //p迭代器指向torsions向量的最后一个元素
	}

	/*需要根据node类型调用相对于的函数
	input——atoms：atom结构体向量，coords：vec结构体向量，c：residue_conf结构体
	*/
	void derivative(const vecv& coords, const vecv& forces, residue_change& c) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		flv::iterator p = c.torsions.begin();
		fl& d = *p; // reference
		++p;
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, d);
		assert(p == c.torsions.end());
	}
};

/*扭矩计数
	input——t：任意类型，s：unint
*/
template<typename T> // T = main_branch, branch, flexible_body
void count_torsions(const T& t, sz& s) {
	t.node.count_torsions(s);
	VINA_FOR_IN(i, t.children)
		count_torsions(t.children[i], s);
}

typedef heterotree<rigid_body> flexible_body;   //形参为rigid_body结构体的heterotree结构体
typedef heterotree<first_segment> main_branch; //形参为first_segment结构体的heterotree结构体

/*
	vector_mutable结构体继承系统vector类

*/
template<typename T> // T == flexible_body || main_branch

struct vector_mutable : public std::vector<T> {
	template<typename C>
	/*获得坐标、起点以及设置四元数orientation_q=q、orientation_m结构体和coords
		input—— atoms：atom结构体向量；coords：vec结构体向量；c：c类型的向量
	*/
	void set_conf(const atomv& atoms, vecv& coords, const std::vector<C>& c) { // C == ligand_conf || residue_conf
		VINA_FOR_IN(i, (*this))                                //for(i=0;i<v.size;i++)
			(*this)[i].set_conf(atoms, coords, c[i]);
	}
	/*扭矩计数
		output——tmp：unint向量
	*/
	szv count_torsions() const {
		szv tmp(this->size(), 0);
		VINA_FOR_IN(i, (*this))
			::count_torsions((*this)[i], tmp[i]);
		return tmp;
	}

	/*需要根据c类型调用相对于的函数
		input——coords、forces：vec结构体向量，c：任意类型向量
	
	*/
	template<typename C>
	void derivative(const vecv& coords, const vecv& forces, std::vector<C>& c) const { // C == ligand_change || residue_change
		VINA_FOR_IN(i, (*this))
			(*this)[i].derivative(coords, forces, c[i]);
	}
};

/* 需要根据t、f类型调用函数         
	input——t、f：任意类型变量；
*/
template<typename T, typename F> // tree or heterotree - like structure
void transform_ranges(T& t, const F& f) {
	t.node.transform(f);
	VINA_FOR_IN(i, t.children)
		transform_ranges(t.children[i], f);
}

#endif
