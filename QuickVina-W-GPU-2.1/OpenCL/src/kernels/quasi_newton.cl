#define EL_TYPE_SIZE 11
#define AD_TYPE_SIZE 20
#define XS_TYPE_SIZE 17
#define SY_TYPE_SIZE 18

inline int num_atom_types(int atu) {
	switch (atu) {
	case 0: return EL_TYPE_SIZE;
	case 1: return AD_TYPE_SIZE;
	case 2: return XS_TYPE_SIZE;
	case 3: return SY_TYPE_SIZE;
	default: printf("Kernel1:num_atom_types() ERROR!"); return INFINITY;
	}
}

inline void elementwise_product(float* out, const float* a, const float* b) {
	out[0] = a[0] * b[0];
	out[1] = a[1] * b[1];
	out[2] = a[2] * b[2];
}

inline float elementwise_product_sum(const float* a, const float* b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline float access_m_data(__constant float* m_data, int m_i, int m_j, int i, int j, int k) {
	return m_data[i + m_i * (j + m_j * k)];
}

inline bool not_max(float x) {
	return (x < 0.1 * INFINITY);// Problem: replace max_fl with INFINITY?
}

inline void curl_with_deriv(float* e, float* deriv, float v, const float epsilon_fl) {
	if (*e > 0 && not_max(v)) {
		float tmp = (v < epsilon_fl) ? 0 : (v / (v + *e));
		*e *= tmp;
		for (int i = 0; i < 3; i++)deriv[i] *= pown(tmp, 2);
	}
}

inline void curl_without_deriv(float* e, float v, const float epsilon_fl) {
	if (*e > 0 && not_max(v)) {
		float tmp = (v < epsilon_fl) ? 0 : (v / (v + *e));
		*e *= tmp;
	}
}

float g_evaluate(			__constant	grid_cl*	g,
					const				float*		m_coords,			// double[3]
					const				float		slope,				// double
					const				float		v,					// double
										float*		deriv,				// double[3]
					const				float		epsilon_fl
) {
	int m_i = g->m_i;
	int m_j = g->m_j;
	int m_k = g->m_k;
	if(m_i * m_j * m_k == 0)printf("\nkernel2: g_evaluate ERROR!#1");
	float tmp_vec[3] = { m_coords[0] - g->m_init[0],m_coords[1] - g->m_init[1] ,m_coords[2] - g->m_init[2] };
	float tmp_vec2[3] = { g->m_factor[0],g->m_factor[1] ,g->m_factor[2] };
	float s[3];
	elementwise_product(s, tmp_vec, tmp_vec2); // 

	float miss[3] = { 0,0,0 };
	int region[3];
	int a[3];
	int m_data_dims[3] = { m_i,m_j,m_k };
	for (int i = 0; i < 3; i++){
		if (s[i] < 0) {
			miss[i] = -s[i];
			region[i] = -1;
			a[i] = 0;
			s[i] = 0;
		}
		else if (s[i] >= g->m_dim_fl_minus_1[i]) {
			miss[i] = s[i] - g->m_dim_fl_minus_1[i];
			region[i] = 1;
			if (m_data_dims[i] < 2)printf("\nKernel2: g_evaluate ERROR!#2");
			a[i] = m_data_dims[i] - 2;
			s[i] = 1;
		}
		else {
			region[i] = 0;
			a[i] = (int)s[i];
			s[i] -= a[i];
		}
		if (s[i] < 0)printf("\nKernel2: g_evaluate ERROR!#3");
		if (s[i] > 1)printf("\nKernel2: g_evaluate ERROR!#4");
		if (a[i] < 0)printf("\nKernel2: g_evaluate ERROR!#5");
		if (a[i] + 1 >= m_data_dims[i])printf("\nKernel2: g_evaluate ERROR!#5");
	}

	float tmp_m_factor_inv[3] = { g->m_factor_inv[0],g->m_factor_inv[1],g->m_factor_inv[2] };
	const float penalty = slope * elementwise_product_sum(miss, tmp_m_factor_inv);
	if (penalty <= -epsilon_fl)printf("\nKernel2: g_evaluate ERROR!#6");

	const int x0 = a[0];
	const int y0 = a[1];
	const int z0 = a[2];
		 
	const int x1 = x0 + 1;
	const int y1 = y0 + 1;
	const int z1 = z0 + 1;

	//const float f000 = access_m_data(g->m_data, m_i, m_j, x0, y0, z0);
	//const float f100 = access_m_data(g->m_data, m_i, m_j, x1, y0, z0);
	//const float f010 = access_m_data(g->m_data, m_i, m_j, x0, y1, z0);
	//const float f110 = access_m_data(g->m_data, m_i, m_j, x1, y1, z0);
	//const float f001 = access_m_data(g->m_data, m_i, m_j, x0, y0, z1);
	//const float f101 = access_m_data(g->m_data, m_i, m_j, x1, y0, z1);
	//const float f011 = access_m_data(g->m_data, m_i, m_j, x0, y1, z1);
	//const float f111 = access_m_data(g->m_data, m_i, m_j, x1, y1, z1);

	int base = (x0 + m_i * (y0 + m_j * z0)) * 8;
	__constant float* base_ptr = &g->m_data[base];

	const float f000 = *base_ptr;
	const float f100 = *(base_ptr + 1);
	const float f010 = *(base_ptr + 2);
	const float f110 = *(base_ptr + 3);
	const float f001 = *(base_ptr + 4);
	const float f101 = *(base_ptr + 5);
	const float f011 = *(base_ptr + 6);
	const float f111 = *(base_ptr + 7);

	const float x = s[0];
	const float y = s[1];
	const float z = s[2];
		  
	const float mx = 1 - x;
	const float my = 1 - y;
	const float mz = 1 - z;

	float f =
		f000 * mx * my * mz +
		f100 * x  * my * mz +
		f010 * mx * y  * mz +
		f110 * x  * y  * mz +
		f001 * mx * my * z	+
		f101 * x  * my * z	+
		f011 * mx * y  * z	+
		f111 * x  * y  * z  ;

	if (deriv) { // valid pointer
		const float x_g =
			f000 * (-1) * my * mz +
			f100 *   1  * my * mz +
			f010 * (-1) * y  * mz +
			f110 *	 1  * y  * mz +
			f001 * (-1) * my * z  +
			f101 *   1  * my * z  +
			f011 * (-1) * y  * z  +
			f111 *   1  * y  * z  ;


		const float y_g =
			f000 * mx * (-1) * mz +
			f100 * x  * (-1) * mz +
			f010 * mx *   1  * mz +
			f110 * x  *   1  * mz +
			f001 * mx * (-1) * z  +
			f101 * x  * (-1) * z  +
			f011 * mx *   1  * z  +
			f111 * x  *   1  * z  ;


		const float z_g =
			f000 * mx * my * (-1) +
			f100 * x  * my * (-1) +
			f010 * mx * y  * (-1) +
			f110 * x  * y  * (-1) +
			f001 * mx * my *   1  +
			f101 * x  * my *   1  +
			f011 * mx * y  *   1  +
			f111 * x  * y  *   1  ;

		float gradient[3] = { x_g, y_g, z_g };

		curl_with_deriv(&f, gradient, v, epsilon_fl);

		float gradient_everywhere[3];

		for (int i = 0; i < 3; i++) {
			gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0);
			deriv[i] = g->m_factor[i] * gradient_everywhere[i] + slope * region[i];
		}
		return f + penalty;
	}	
	else {  // none valid pointer
		printf("\nKernel2: g_evaluate ERROR!#7");
		curl_without_deriv(&f, v, epsilon_fl);
		return f + penalty;
	}
}


float ig_eval_deriv(						output_type_cl*		x,
											change_cl*			g, 
						const				float				v,
								__constant	grids_cl*			grids,
											m_cl*				m,
						const				float				epsilon_fl
) {
	float e = 0;
	int nat = num_atom_types(grids->atu);
	for (int i = 0; i < m->m_num_movable_atoms; i++) {
		int t = m->atoms[i].types[grids->atu];
		if (t >= nat) {
			for (int j = 0; j < 3; j++)m->minus_forces.coords[i][j] = 0;
			continue;
		}
		float deriv[3];

		e = e + g_evaluate(&grids->grids[t], m->m_coords.coords[i], grids->slope, v, deriv, epsilon_fl);

		for (int j = 0; j < 3; j++) m->minus_forces.coords[i][j] = deriv[j];
	}
	return e;
}

inline void quaternion_to_r3(const float* q, float* orientation_m) {
	// Omit assert(quaternion_is_normalized(q));
	const float a = q[0];
	const float b = q[1];
	const float c = q[2];
	const float d = q[3];

	const float aa = a * a;
	const float ab = a * b;
	const float ac = a * c;
	const float ad = a * d;
	const float bb = b * b;
	const float bc = b * c;
	const float bd = b * d;
	const float cc = c * c;
	const float cd = c * d;
	const float dd = d * d;

	//Omit assert(eq(aa + bb + cc + dd, 1));
	matrix tmp;
	mat_init(&tmp, 0); // matrix with fixed dimension 3(here we treate this as a regular matrix(not triangular matrix!))

	matrix_set_element(&tmp, 3, 0, 0,		(aa + bb - cc - dd)	);
	matrix_set_element(&tmp, 3, 0, 1, 2 *	(-ad + bc)			);
	matrix_set_element(&tmp, 3, 0, 2, 2 *	(ac + bd)			);
							 
	matrix_set_element(&tmp, 3, 1, 0, 2 *	(ad + bc)			);
	matrix_set_element(&tmp, 3, 1, 1,		(aa - bb + cc - dd)	);
	matrix_set_element(&tmp, 3, 1, 2, 2 *	(-ab + cd)			);
							 
	matrix_set_element(&tmp, 3, 2, 0, 2 *	(-ac + bd)			);
	matrix_set_element(&tmp, 3, 2, 1, 2 *	(ab + cd)			);
	matrix_set_element(&tmp, 3, 2, 2,		(aa - bb - cc + dd)	);

	for (int i = 0; i < 9; i++) orientation_m[i] = tmp.data[i];
}

inline void local_to_lab_direction(			float* out,
									const	float* local_direction,
									const	float* orientation_m
) {
	out[0] =	orientation_m[0] * local_direction[0] +
				orientation_m[3] * local_direction[1] +
				orientation_m[6] * local_direction[2];
	out[1] =	orientation_m[1] * local_direction[0] +
				orientation_m[4] * local_direction[1] +
				orientation_m[7] * local_direction[2];
	out[2] =	orientation_m[2] * local_direction[0] +
				orientation_m[5] * local_direction[1] +
				orientation_m[8] * local_direction[2];
}

inline void local_to_lab(						float*		out,
							const				float*		origin,
							const				float*		local_coords,
							const				float*		orientation_m
) {
	out[0] = origin[0] + (	orientation_m[0] * local_coords[0] +
							orientation_m[3] * local_coords[1] +
							orientation_m[6] * local_coords[2]
							);			 
	out[1] = origin[1] + (	orientation_m[1] * local_coords[0] +
							orientation_m[4] * local_coords[1] +
							orientation_m[7] * local_coords[2]
							);			 
	out[2] = origin[2] + (	orientation_m[2] * local_coords[0] +
							orientation_m[5] * local_coords[1] +
							orientation_m[8] * local_coords[2]
							);
}

inline void angle_to_quaternion2(				float*		out,
									const		float*		axis,
												float		angle
) {
	if (norm3(axis) - 1 >= 0.001)printf("\nkernel2: angle_to_quaternion() ERROR!"); // Replace assert(eq(axis.norm(), 1));
	normalize_angle(&angle);
	float c = cos(angle / 2);
	float s = sin(angle / 2);
	out[0] = c;
	out[1] = s * axis[0];
	out[2] = s * axis[1];
	out[3] = s * axis[2];
}

void set(	const				output_type_cl* x,
								rigid_cl*		lig_rigid_gpu,
								m_coords_cl*	m_coords_gpu,	
			const				atom_cl*		atoms,				
			const				int				m_num_movable_atoms,
			const				float			epsilon_fl
) {
	//************** (root) node.set_conf **************// (CHECKED)
	for (int i = 0; i < 3; i++) lig_rigid_gpu->origin[0][i] = x->position[i]; // set origin
	for (int i = 0; i < 4; i++) lig_rigid_gpu->orientation_q[0][i] = x->orientation[i]; // set orientation_q
	quaternion_to_r3(lig_rigid_gpu->orientation_q[0], lig_rigid_gpu->orientation_m[0]);// set orientation_m
	// set coords
	int begin = lig_rigid_gpu->atom_range[0][0];
	int end =	lig_rigid_gpu->atom_range[0][1];
	for (int i = begin; i < end; i++) {
		local_to_lab(m_coords_gpu->coords[i], lig_rigid_gpu->origin[0], &atoms[i].coords[0], lig_rigid_gpu->orientation_m[0]);
	}
	//************** end node.set_conf **************//

	//************** branches_set_conf **************//
	//update nodes in depth-first order
	for (int current = 1; current < lig_rigid_gpu->num_children + 1; current++) { // current starts from 1 (namely starts from first child node)
		int parent = lig_rigid_gpu->parent[current];
		float torsion = x->lig_torsion[current - 1]; // torsions are all related to child nodes
		local_to_lab(	lig_rigid_gpu->origin[current],
						lig_rigid_gpu->origin[parent],
						lig_rigid_gpu->relative_origin[current],
						lig_rigid_gpu->orientation_m[parent]
						); // set origin
		local_to_lab_direction(	lig_rigid_gpu->axis[current],
								lig_rigid_gpu->relative_axis[current],
								lig_rigid_gpu->orientation_m[parent]
								); // set axis
		float tmp[4];
		float parent_q[4] = {	lig_rigid_gpu->orientation_q[parent][0],
								lig_rigid_gpu->orientation_q[parent][1] ,
								lig_rigid_gpu->orientation_q[parent][2] ,
								lig_rigid_gpu->orientation_q[parent][3] };
		float current_axis[3] = {	lig_rigid_gpu->axis[current][0],
									lig_rigid_gpu->axis[current][1],
									lig_rigid_gpu->axis[current][2] };

		angle_to_quaternion2(tmp, current_axis, torsion);
		angle_to_quaternion_multi(tmp, parent_q);
		quaternion_normalize_approx(tmp, epsilon_fl);

		for (int i = 0; i < 4; i++) lig_rigid_gpu->orientation_q[current][i] = tmp[i]; // set orientation_q
		quaternion_to_r3(lig_rigid_gpu->orientation_q[current], lig_rigid_gpu->orientation_m[current]); // set orientation_m

		// set coords
		begin = lig_rigid_gpu->atom_range[current][0];
		end =	lig_rigid_gpu->atom_range[current][1];
		for (int i = begin; i < end; i++) {
			local_to_lab(m_coords_gpu->coords[i], lig_rigid_gpu->origin[current], &atoms[i].coords[0], lig_rigid_gpu->orientation_m[current]);
		}
	}
	//************** end branches_set_conf **************//
}

void p_eval_deriv(						float*		out,
										int			type_pair_index,
										float		r2,
							__constant	pre_cl*		pre,
					const				float		epsilon_fl
) {
	const float cutoff_sqr = pre->m_cutoff_sqr;
	if(r2 > cutoff_sqr) printf("\nkernel2: p_eval_deriv() ERROR!");
	__constant p_m_data_cl* tmp = &pre->m_data[type_pair_index];
	float r2_factored = tmp->factor * r2;
	if (r2_factored + 1 >= SMOOTH_SIZE) printf("\nkernel2: p_eval_deriv() ERROR!");
	int i1 = (int)(r2_factored);
	int i2 = i1 + 1;
	if (i1 >= SMOOTH_SIZE || i1 < 0)printf("\n kernel2: p_eval_deriv() ERROR!");
	if (i2 >= SMOOTH_SIZE || i2 < 0)printf("\n : p_eval_deriv() ERROR!");
	float rem = r2_factored - i1;
	if (rem < -epsilon_fl)printf("\nkernel2: p_eval_deriv() ERROR!");
	if (rem >= 1 + epsilon_fl)printf("\nkernel2: p_eval_deriv() ERROR!");
	float p1[2] = { tmp->smooth[i1][0], tmp->smooth[i1][1] };
	float p2[2] = { tmp->smooth[i2][0], tmp->smooth[i2][1] };
	float e = p1[0] + rem * (p2[0] - p1[0]);
	float dor = p1[1] + rem * (p2[1] - p1[1]);
	out[0] = e;
	out[1] = dor;
}

inline void curl(float* e, float* deriv, float v, const float epsilon_fl) {
	if (*e > 0 && not_max(v)) {
		float tmp = (v < epsilon_fl) ? 0 : (v / (v + *e));
		(*e) = tmp * (*e);
		for (int i = 0; i < 3; i++)deriv[i] = deriv[i] * (tmp * tmp);
	}
}

float eval_interacting_pairs_deriv(			__constant	pre_cl*			pre,
									const				float			v,
									const				lig_pairs_cl*   pairs,
									const			 	m_coords_cl*	m_coords,
														m_minus_forces* minus_forces,
									const				float			epsilon_fl
) {
	float e = 0;

	for (int i = 0; i < pairs->num_pairs; i++) {
		const int ip[3] = { pairs->type_pair_index[i], pairs->a[i] ,pairs->b[i] };
		float coords_b[3] = { m_coords->coords[ip[2]][0], m_coords->coords[ip[2]][1], m_coords->coords[ip[2]][2] };
		float coords_a[3] = { m_coords->coords[ip[1]][0], m_coords->coords[ip[1]][1], m_coords->coords[ip[1]][2] };
		float r[3] = { coords_b[0] - coords_a[0], coords_b[1] - coords_a[1] ,coords_b[2] - coords_a[2] };
		float r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
	
		if (r2 < pre->m_cutoff_sqr) {
			float tmp[2];
			p_eval_deriv(tmp, ip[0], r2, pre, epsilon_fl);
			float force[3] = { r[0] * tmp[1], r[1] * tmp[1] ,r[2] * tmp[1] };
			curl(&tmp[0], force, v, epsilon_fl);
			e += tmp[0];
			for (int j = 0; j < 3; j++)minus_forces->coords[ip[1]][j] -= force[j];
			for (int j = 0; j < 3; j++)minus_forces->coords[ip[2]][j] += force[j];
		}
	}
	return e;
}

inline void product(float* res, const float*a,const float*b) {
	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}

void POT_deriv(	const					m_minus_forces* minus_forces,
				const					rigid_cl*		lig_rigid_gpu,
				const					m_coords_cl*		m_coords,
										change_cl*		g
) {
	int num_torsion = lig_rigid_gpu->num_children;
	int num_rigid = num_torsion + 1;
	float position_derivative_tmp[MAX_NUM_OF_RIGID][3];
	float position_derivative[MAX_NUM_OF_RIGID][3];
	float orientation_derivative_tmp[MAX_NUM_OF_RIGID][3];
	float orientation_derivative[MAX_NUM_OF_RIGID][3];
	float torsion_derivative[MAX_NUM_OF_RIGID]; // torsion_derivative[0] has no meaning(root node has no torsion)

	for (int i = 0; i < num_rigid; i++) {
		int begin = lig_rigid_gpu->atom_range[i][0];
		int end = lig_rigid_gpu->atom_range[i][1];
		for (int k = 0; k < 3; k++)position_derivative_tmp[i][k] = 0; 
		for (int k = 0; k < 3; k++)orientation_derivative_tmp[i][k] = 0;
		for (int j = begin; j < end; j++) {
			for (int k = 0; k < 3; k++)position_derivative_tmp[i][k] += minus_forces->coords[j][k];

			float tmp1[3] = {	m_coords->coords[j][0] - lig_rigid_gpu->origin[i][0],
								m_coords->coords[j][1] - lig_rigid_gpu->origin[i][1],
								m_coords->coords[j][2] - lig_rigid_gpu->origin[i][2] };
			float tmp2[3] = {  minus_forces->coords[j][0],
								minus_forces->coords[j][1],
								minus_forces->coords[j][2] };
			float tmp3[3];
			product(tmp3, tmp1, tmp2);
			for (int k = 0; k < 3; k++)orientation_derivative_tmp[i][k] += tmp3[k];
		}
	}

	// position_derivative 
	for (int i = num_rigid - 1; i >= 0; i--) {// from bottom to top
		for (int k = 0; k < 3; k++)position_derivative[i][k] = position_derivative_tmp[i][k]; // self
		// looking for chidren node
		for (int j = 0; j < num_rigid; j++) {
			if (lig_rigid_gpu->children_map[i][j] == true) {
				for (int k = 0; k < 3; k++)position_derivative[i][k] += position_derivative[j][k]; // self+children node
			}
		}
	}

	// orientation_derivetive
	for (int i = num_rigid - 1; i >= 0; i--) { // from bottom to top
		for (int k = 0; k < 3; k++)orientation_derivative[i][k] = orientation_derivative_tmp[i][k]; // self
		// looking for chidren node
		for (int j = 0; j < num_rigid; j++) {
			if (lig_rigid_gpu->children_map[i][j] == true) { // self + children node + product
				for (int k = 0; k < 3; k++)orientation_derivative[i][k] += orientation_derivative[j][k];
				float product_out[3];
				float origin_temp[3] = {	lig_rigid_gpu->origin[j][0] - lig_rigid_gpu->origin[i][0],
											lig_rigid_gpu->origin[j][1] - lig_rigid_gpu->origin[i][1],
											lig_rigid_gpu->origin[j][2] - lig_rigid_gpu->origin[i][2] };
				product(product_out, origin_temp, position_derivative[j]);
				for (int k = 0; k < 3; k++)orientation_derivative[i][k] += product_out[k];
			}
		}
	}

	// torsion_derivative
	for (int i = num_rigid - 1; i >= 0; i--) { // from bottom to top
		float sum = 0;
		for (int j = 0; j < 3; j++) sum += orientation_derivative[i][j] * lig_rigid_gpu->axis[i][j];
		torsion_derivative[i] = sum;
	}

	for (int k = 0; k < 3; k++)	g->position[k] = position_derivative[0][k];
	for (int k = 0; k < 3; k++) g->orientation[k] = orientation_derivative[0][k];
	for (int k = 0; k < num_torsion; k++) g->lig_torsion[k] = torsion_derivative[k + 1];// no meaning for node 0
}


float m_eval_deriv(						output_type_cl*		c,
										change_cl*			g,
										m_cl*				m,
							__constant	pre_cl*				pre,
							__constant	grids_cl*			grids,
					const	__global	float*				v,
					const				float				epsilon_fl
) {
	set(c, &m->ligand.rigid, &m->m_coords, m->atoms, m->m_num_movable_atoms, epsilon_fl);

	float e = ig_eval_deriv(	c,
								g, 
								v[1],				
								grids,
								m,
								epsilon_fl							
							);
	
	e += eval_interacting_pairs_deriv(	pre,
										v[0],
										&m->ligand.pairs,
										&m->m_coords,
										&m->minus_forces,
										epsilon_fl
									);

	POT_deriv(&m->minus_forces, &m->ligand.rigid, &m->m_coords, g);

	return e;
}


// Only support one ligand, no flex !
inline float find_change_index_read(const change_cl* g, int index, int lig_torsion_size) {
	if (index < 3)return g->position[index];
	index -= 3;
	if (index < 3)return g->orientation[index];
	index -= 3;
	if (index < lig_torsion_size)return g->lig_torsion[index]; 
	printf("\nKernel2:find_change_index_read() ERROR!"); // Shouldn't be here
	return -1;
}

inline void find_change_index_write(change_cl* g, int index, float data, int lig_torsion_size) {
	if (index < 3) { g->position[index] = data; return; }
	index -= 3;
	if (index < 3) { g->orientation[index] = data; return; }
	index -= 3;
	if (index < lig_torsion_size) { g->lig_torsion[index] = data; return; } 
	printf("\nKernel2:find_change_index_write() ERROR!"); // Shouldn't be here
}

void minus_mat_vec_product(	const		matrix*		h,
							const		change_cl*	in,
										change_cl*  out,
							const		int			lig_torsion_size
) {
	int n = h->dim;
	for (int i = 0; i < n; i++) {
		float sum = 0;
		for (int j = 0; j < n; j++) {
			sum += h->data[index_permissive(h, i, j)] * find_change_index_read(in, j, lig_torsion_size);
		}
		find_change_index_write(out, i, -sum, lig_torsion_size);
	}
}


inline float scalar_product(	const	change_cl*	a,
								const	change_cl*	b,
								const	int			n,
								const	int			lig_torsion_size
) {
	float tmp = 0;
	for (int i = 0; i < n; i++) {
		tmp += find_change_index_read(a, i, lig_torsion_size) * find_change_index_read(b, i, lig_torsion_size);
	}
	return tmp;
}

float line_search(					 	m_cl*				m,
							__constant	pre_cl*				pre,
							__constant	grids_cl*			grids,
										int					n,
					const				output_type_cl*		x,
					const				change_cl*			g,
					const				float				f0,
					const				change_cl*			p,
										output_type_cl*		x_new,
										change_cl*			g_new,
										float*				f1,
					const				float				epsilon_fl,
					const	__global	float*				hunt_cap,
					const				int					lig_torsion_size
) {
	const float c0 = 0.0001;
	const int max_trials = 10;
	const float multiplier = 0.5;
	float alpha = 1;

	const float pg = scalar_product(p, g, n, lig_torsion_size);

	for (int trial = 0; trial < max_trials; trial++) {

		*x_new =  *x;

		output_type_cl_increment(x_new, p, alpha, epsilon_fl, lig_torsion_size);

		*f1 =  m_eval_deriv(x_new,
							g_new,
							m,
							pre,
							grids,
							hunt_cap,
							epsilon_fl
							);

		if (*f1 - f0 < c0 * alpha * pg)
			break;
		alpha *= multiplier;
	}
	return alpha;
}


bool bfgs_update(			matrix*			h,
					const	change_cl*		p,
					const	change_cl*		y,
					const	float			alpha,
					const	mis_cl*			mis,
					const	int				torsion_size
) {
	const float yp = scalar_product(y, p, h->dim, torsion_size);

	if (alpha * yp < mis->epsilon_fl) return false;

	change_cl minus_hy = *y;
	minus_mat_vec_product(h, y, &minus_hy, torsion_size);
	const float yhy = -scalar_product(y, &minus_hy, h->dim, torsion_size);
	const float r = 1 / (alpha * yp);
	const int n = 6 + torsion_size;
	int s = torsion_size;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			float tmp = alpha * r * (find_change_index_read(&minus_hy, i, s) * find_change_index_read(p, j,s)
									+ find_change_index_read(&minus_hy, j, s) * find_change_index_read(p, i,s)) +
									+alpha * alpha * (r * r * yhy + r) * find_change_index_read(p, i,s) * find_change_index_read(p, j, s);

			h->data[i + j * (j + 1) / 2] += tmp;
		}
	}

	return true;
}


void bfgs(					output_type_cl*			x,
								change_cl*			g,
								m_cl*				m,
					__constant	pre_cl*				pre,
					__constant	grids_cl*			grids,
			const	__global	mis_cl*				mis,
			const				int					torsion_size,
			const				int					max_bfgs_steps,
							bool					global_flag,
							global_container*		g_container,
				__global   individual_container*	circularvisited,
							int						search_depth,
							vec3_cl*				origin,
							vec3_cl*				box_size,
				__global    ele_cl*					global_ptr,
				__global    int*					count_id
) 
{
	int lig_torsion_size = torsion_size;
	int n = 3 + 3 + lig_torsion_size; // the dimensions of matirx
	int thread = mis->thread;

	matrix h;
	matrix_init(&h, n, 0);
	matrix_set_diagonal(&h, 1);

	change_cl g_new = *g;

	output_type_cl x_new = *x ;
	float f0 = m_eval_deriv(	x,
								g,
								m,
								pre,
								grids,
								mis->hunt_cap,
								mis->epsilon_fl
							);
	int ret = 0;
	ret = global_interesting_cl(g_container, x, f0, g, 0, thread, search_depth, origin, box_size, global_ptr, count_id, torsion_size);
	
	if (ret >= 0 && individual_interesting_cl(circularvisited, x, f0, g, ret, torsion_size) >= 0) {
		x->e = f0;
	}
	else
	{
		add_to_individual_buffer(circularvisited, x, f0, g, torsion_size);
		float f_orig = f0;

		change_cl g_orig = *g;
		output_type_cl x_orig = *x;

		change_cl p = *g;

		float f_values[MAX_NUM_OF_BFGS_STEPS + 1];
		f_values[0] = f0;

		for (int step = 0; step < max_bfgs_steps; step++) {

			minus_mat_vec_product(&h, g, &p, torsion_size);
			float f1 = 0;

			const float alpha = line_search(	m,
												pre,
												grids,
												n,
												x,
												g,
												f0,
												&p,
												&x_new,
												&g_new,
												&f1,
												mis->epsilon_fl,
												mis->hunt_cap,
												torsion_size
											);
			
			change_cl y = g_new;
			// subtract_change
			for (int i = 0; i < n; i++) {
				float tmp = find_change_index_read(&y, i, torsion_size) - find_change_index_read(g, i, torsion_size);
				find_change_index_write(&y, i, tmp, torsion_size);
			}

			f_values[step + 1] = f1;

			f0 = f1;
			*x = x_new;
			if (!(sqrt(scalar_product(g, g, n, torsion_size)) >= 1e-5))break;
			*g = g_new;

			if (step == 0) {
				float yy = scalar_product(&y, &y, n, torsion_size);
				if (fabs(yy) > mis->epsilon_fl) {
					matrix_set_diagonal(&h, alpha * scalar_product(&y, &p, n, torsion_size) / yy);
				}
			}
			bool h_updated = bfgs_update(&h, &p, &y, alpha, mis, torsion_size);
			add_to_individual_buffer(circularvisited, x, f0, g, torsion_size);
		}

		if (!global_flag) {
			add_to_global_buffer(g_container, thread, search_depth, x, f0, g, origin, box_size, global_ptr, count_id, torsion_size);
		}

		if (!(f0 <= f_orig)) {
			f0 = f_orig;
			*x = x_orig;
			*g = g_orig;
		}

		// write output_type_cl energy
		x->e = f0;
	}
	
}

// add zxy start//
inline void to_minus(change_cl* out, const int n)
{
	// out = -out
	for (int i = 0; i < 3; i++)out->position[i] = -out->position[i];
	for (int i = 0; i < 3; i++)out->orientation[i] = -out->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)out->lig_torsion[i] = -out->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)out->flex_torsion[i] = -out->flex_torsion[i];
}
inline void get_to_minus(change_cl* a, change_cl* b, const int n)
{
	for (int i = 0; i < 3; i++)a->position[i] = -b->position[i];
	for (int i = 0; i < 3; i++)a->orientation[i] = -b->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)a->lig_torsion[i] = -b->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)a->flex_torsion[i] = -b->flex_torsion[i];
}
inline float to_norm(const change_cl* in, const int n, const int torsion_size)
{
	float d_test = 0;
	for (int i = 0; i < n; i++)
	{
		d_test += find_change_index_read(in, i, torsion_size) * find_change_index_read(in, i, torsion_size);
	}
	d_test = sqrt(d_test);
	return d_test;
}


inline int line_search_lewisoverton(m_cl* m_cl_gpu,
	__constant pre_cl* p_cl_gpu,
	__constant grids_cl* ig_cl_gpu,
	int n,
	float* stp,
	output_type_cl* x,
	float* f,
	change_cl* g,
	const change_cl* d,
	const output_type_cl* xp,
	const change_cl* gp,
	const float epsilon_fl,
	const __global float* hunt_cap,
	const int torsion_size

	)
{
	int count = 0;
	bool brackt = false;
	float finit, dginit, dgtest, dstest;
	float mu = 0.0, nu = 1.0e+20;

	dginit = scalar_product(gp, d, n, torsion_size);

	/* Make sure that s points to a descent direction. */
	if (0.0 < dginit)
	{
		return -1;
	}
	/* The initial value of the cost function. */
	finit = *f;
	//F_DEC_COEFF
	dgtest = 1.0e-4 * dginit;
	// S_CURV_COEFF 
	dstest = 0.1 * dginit;

	while (true)
	{
		//output_type_cl_init_with_output(x, xp);
		*x = *xp;
		output_type_cl_increment(x, d, *stp, epsilon_fl, torsion_size);
		/* Evaluate the function and gradient values. */
		*f = m_eval_deriv(x,
			g,
			m_cl_gpu,
			p_cl_gpu,
			ig_cl_gpu,
			hunt_cap,
			epsilon_fl
		);
		++count;
		/* Check the Armijo condition. */
		if (*f > finit + *stp * dgtest)
		{
			nu = *stp;
			brackt = true;
		}
		else
		{
			if (scalar_product(g, d, n, torsion_size) < dstest)
			{
				mu = *stp;
			}
			else
			{
				return count;
			}
		}
		// MAX_LINESEARCH
		if (10 <= count)
		{
			if (*f > finit)
				return -1;
			else
				return count;
		}
		if (brackt)
		{
			(*stp) = 0.5 * (mu + nu);
		}
		else
		{
			(*stp) *= 2.0;
		}
	}
}
void rilc_bfgs(					output_type_cl*			x,
								change_cl*			g,
								m_cl*				m,
					__constant	pre_cl*				pre,
					__constant	grids_cl*			grids,
			const	__global	mis_cl*				mis,
			const				int					torsion_size,
			const				int					max_bfgs_steps,
							bool					global_flag,
							global_container*		g_container,
				__global   individual_container*	circularvisited,
							int						search_depth,
							vec3_cl*				origin,
							vec3_cl*				box_size,
				__global    ele_cl*					global_ptr,
				__global    int*					count_id
) 
{

	int thread = mis->thread;

	output_type_cl xp = *x;
	change_cl gp = *g;
	float fx = m_eval_deriv(	x,
								g,
								m,
								pre,
								grids,
								mis->hunt_cap,
								mis->epsilon_fl
							);
	int ret = 0;
	ret = global_interesting_cl(g_container, x, fx, g, 0, thread, search_depth, origin, box_size, global_ptr, count_id, torsion_size);
	
	if (ret >= 0 && individual_interesting_cl(circularvisited, x, fx, g, ret, torsion_size) >= 0) {
		x->e = fx;
	}
	else
	{
		add_to_individual_buffer(circularvisited, x, fx, g, torsion_size);
		int n = 3 + 3 + torsion_size; // the dimensions of matirx
		// rilc_bfgs init
		int k, ls;
		float step, ys, yy, step_min, step_max;
		float beta = 0, cau;
		float lm_s_dot_d, cau_t;
		float d_update_tmp = 0;

		float lm_alpha = 0;
		float lm_s[MAX_HESSIAN_MATRIX_SIZE] = { 0 };
		float lm_y[MAX_HESSIAN_MATRIX_SIZE] = { 0 };
		float lm_ys = 0;
		
		float fxp = fx;
		float fx_orig = fx;
		change_cl g_orig = *g;
		output_type_cl x_orig = *x;

		change_cl d = *g;

		//d = -g
		get_to_minus(&d, g, n);
		if (!(sqrt(scalar_product(g, g, n, torsion_size)) >= 1e-5))
		{
			x->e = fx;
			return;
		}
		else
		{
			step = 1.0 / to_norm(&d, n, torsion_size);
			k = 1;

			while (true) {

				xp = *x;
				gp = *g;

				fxp = fx;


				ls = line_search_lewisoverton(m,
					pre,
					grids,
					n,
					&step,
					x,
					&fx,
					g,
					&d,
					&xp,
					&gp,
					mis->epsilon_fl,
					mis->hunt_cap,
					torsion_size
				);
				if (ls < 0) {

					*x = xp;
					*g = gp;

					fx = fxp;
					add_to_individual_buffer(circularvisited, x, fx, g, torsion_size);
					break;
				}
				if (!(sqrt(scalar_product(g, g, n, torsion_size)) >= 1e-5)) {
					x->e = fx;
					return;
				}
				if (max_bfgs_steps != 0 && max_bfgs_steps <= k) {
					add_to_individual_buffer(circularvisited, x, fx, g, torsion_size);
					break;
				}
				/* Count the iteration number. */
				++k;

				for (int i = 0; i < n; i++) {
					lm_s[i] = step * find_change_index_read(&d, i, torsion_size);
				}
				for (int i = 0; i < n; i++) {
					lm_y[i] = find_change_index_read(g, i, torsion_size) - find_change_index_read(&gp, i, torsion_size);
				}

				ys = 0;
				for (int i = 0; i < n; i++) {
					ys += lm_y[i] * lm_s[i];
				}

				yy = 0;
				for (int i = 0; i < n; i++) {
					yy += lm_y[i] * lm_y[i];
				}

				lm_ys = ys;

				/* Compute the negative of gradients. */
				// d = -g
				get_to_minus(&d, g, n);
				cau = 0;
				for (int i = 0; i < n; i++) {
					cau += lm_s[i] * lm_s[i];
				}
				cau_t = 0;
				for (int i = 0; i < n; i++)
				{
					cau_t = cau_t + (find_change_index_read(&gp, i, torsion_size) * find_change_index_read(&gp, i, torsion_size));
				}
				cau = cau * sqrt(cau_t);
				// CAUTIOUS_FACTOR
				cau = cau * 1.0e-6;

				if (ys > cau) {
					/* \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}. */
					lm_s_dot_d = 0;
					for (int a = 0; a < n; a++)
					{
						lm_s_dot_d = lm_s_dot_d + lm_s[a] * find_change_index_read(&d, a, torsion_size);
					}
					lm_alpha = lm_s_dot_d / lm_ys;
					for (int b = 0; b < n; b++)
					{
						float mimus_lm_alpha_mulit_lm_y = (-lm_alpha * lm_y[b]);
						find_change_index_write(&d, b, find_change_index_read(&d, b, torsion_size) + lm_alpha * lm_y[b], torsion_size);
					}


					// d *= ys / yy;
					for (int i = 0; i < n; i++)
					{
						find_change_index_write(&d, i, find_change_index_read(&d, i, torsion_size) * (ys / yy), torsion_size);
					}

					beta = 0;
					for (int a = 0; a < n; a++)
					{
						beta += (lm_y[a] * find_change_index_read(&d, a, torsion_size));
					}
					beta /= lm_ys;


					for (int i = 0; i < n; i++)
					{
						find_change_index_write(&d, i, find_change_index_read(&d, i, torsion_size) + ((lm_alpha - beta) * lm_s[i]), torsion_size);
					}

				}
				step = 1.0;
				add_to_individual_buffer(circularvisited, x, fx, g, torsion_size);
			}

		}
		
		if (!global_flag) {
			add_to_global_buffer(g_container, thread, search_depth, x, fx, g, origin, box_size, global_ptr, count_id, torsion_size);
		}

		if (!(fx <= fx_orig)) {
			fx = fx_orig;

			*x = x_orig;

			*g = g_orig;

		}

		// write output_type_cl energy
		x->e = fx;
	}
	
}