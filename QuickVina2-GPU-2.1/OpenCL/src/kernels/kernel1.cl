int fl_to_sz(float x, float max_sz) {
	if (x <= 0) return 0;
	if (x >= max_sz) return max_sz;
	return (int)x;
}
int num_atom_types(int atu) {
	switch (atu) {
	case 0: return EL_TYPE_SIZE;
	case 1: return AD_TYPE_SIZE;
	case 2: return XS_TYPE_SIZE;
	case 3: return SY_TYPE_SIZE;
	default: printf("\nkernel1:num_atom_types ERROR!"); return INFINITY;
	}
}

const float vec_distance_sqr(const float* a, const float* b) {
	return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
}

const int triangular_matrix_index(int n, int i, int j) {
	if (j >= n) printf("\nkernel1:triangular_matrix_index ERROR!");
	if (i > j) printf("\nkernel1:triangular_matrix_index ERROR!");
	return i + j * (j + 1) / 2;
}

const int triangular_matrix_index_permissive(int n, int i, int j) {
	return (i <= j) ? triangular_matrix_index(n, i, j) : triangular_matrix_index(n, j, i);
}

float eval_fast(int type_pair_index, float r2, float cutoff_sqr, const __global pre_cl* pre) {
	if (r2 > cutoff_sqr) printf("\nkernel1:eval_fast ERROR!");
	if (r2 * pre->factor >= FAST_SIZE)printf("\nkernel1:eval_fast ERROR!");
	int i = (int)(pre->factor * r2);
	if (i >= FAST_SIZE)printf("\nkernel1:eval_fast ERROR!");
	float res = pre->m_data[type_pair_index].fast[i];
	return res;
}

const __global int* possibilities(						float*		coords,
												const	__global	ar_cl*		ar,
												const				float		epsilon_fl, 
												const	__global	gb_cl*		gb,
																	int*		relation_count,
												const	__global	mis_cl*		mis
) {
	int index[3];
	float temp_array[3];
	int m_data_dims[3] = { mis->ar_mi, mis->ar_mj, mis->ar_mk };
	for (int i = 0; i < 3; i++) {
		if (coords[i] + epsilon_fl < gb->init[i]) printf("\nkernel1:possibilities ERROR!1 coords = %f,  init = %f", coords[i], gb->init[i]);//replace assert()
		if (coords[i] > gb->init[i] + gb->range[i] + epsilon_fl) printf("\nkernel1:possibilities ERROR!2 coords = %f,  init = %f, range =%f", coords[i], gb->init[i], gb->range[i]);//replace assert()
		const float tmp = (coords[i] - gb->init[i]) * m_data_dims[i] / gb->range[i];
		temp_array[i] = tmp;
		index[i] = fl_to_sz(tmp, m_data_dims[i] - 1);
	}
	int temp = index[0] + m_data_dims[0] * (index[1] + m_data_dims[1] * index[2]);
	const __global int* address;
	address = &(ar->relation[temp][0]);
	*relation_count = ar->relation_size[temp];
	return address;
}


__kernel
void kernel1(
				const		__global		pre_cl*		pre,
				const		__global		pa_cl*		pa,
				const		__global		gb_cl*		gb,
				const		__global		ar_cl*		ar,
							__global		grids_cl*	grids,
				const		__global		mis_cl*		mis,
				const		__global		float*		needed,
				const						int			atu,
				const						int			nat
) {
	int x = get_global_id(0);
	int y = get_global_id(1);
	int z = get_global_id(2);
	int grids_front = mis->grids_front;

	//printf("here!");
	if (x >=grids->grids[grids_front].m_i)return;
	if (y >=grids->grids[grids_front].m_j)return;
	if (z >=grids->grids[grids_front].m_k)return;
	
	float affinities[17]; if (mis->needed_size > 17) printf("\nkernel1: ERROR!");
	for (int i = 0; i < mis->needed_size; i++) affinities[i] = 0;

	float probe_coords[3] = {	grids->grids[grids_front].m_init[0] + grids->grids[grids_front].m_factor_inv[0] * x,
								grids->grids[grids_front].m_init[1] + grids->grids[grids_front].m_factor_inv[1] * y,
								grids->grids[grids_front].m_init[2] + grids->grids[grids_front].m_factor_inv[2] * z };
	int relation_count;

	const __global int* possibilities_ptr = possibilities(probe_coords, ar, mis->epsilon_fl, gb, &relation_count, mis);

	for (int possibilities_i = 0; possibilities_i < relation_count; possibilities_i++) {
		const int i = possibilities_ptr[possibilities_i];
		const atom_cl* a = &pa->atoms[i];
		const int t1 = a->types[atu];
		if (t1 >= nat) continue;
		const float r2 = vec_distance_sqr(a->coords, probe_coords);

		if (r2 <= mis->cutoff_sqr) {
			for (int j = 0; j < mis->needed_size; j++) {
				const int t2 = needed[j]; if (t2 >= nat) printf("\nkernel1: ERROR!1 t2 = %d,nat=%d", t2, nat);
				const int type_pair_index = triangular_matrix_index_permissive(num_atom_types(atu), t1, t2);
				affinities[j] += eval_fast(type_pair_index, r2, mis->cutoff_sqr, pre);
			}
		}
	}


	for (int j = 0; j < mis->needed_size; j++) {
		int t = needed[j]; if (t >= nat)printf("\nkernel1: ERROR!2 t = %d,nat=%d",t,nat);
		int mi = grids->grids[t].m_i;
		int mj = grids->grids[t].m_j;

		int addr_base =  x		+ mi * ( y      + mj *  z     );
		int addr_1	  = (x - 1) + mi * ( y      + mj *  z     );
		int addr_2    =  x      + mi * ((y - 1) + mj *  z     );
		int addr_3    = (x - 1) + mi * ((y - 1) + mj *  z     );
		int addr_4    =  x      + mi * ( y      + mj * (z - 1));
		int addr_5    = (x - 1) + mi * ( y      + mj * (z - 1));
		int addr_6    =  x      + mi * ((y - 1) + mj * (z - 1));
		int addr_7    = (x - 1) + mi * ((y - 1) + mj * (z - 1));

								grids->grids[t].m_data[addr_base * 8    ] = affinities[j];
		if (x != 0)				grids->grids[t].m_data[addr_1    * 8 + 1] = affinities[j];
		if (y != 0)				grids->grids[t].m_data[addr_2    * 8 + 2] = affinities[j];
		if (x != 0 && y != 0)	grids->grids[t].m_data[addr_3    * 8 + 3] = affinities[j];
		if (z != 0)				grids->grids[t].m_data[addr_4    * 8 + 4] = affinities[j];
		if (x != 0 && z != 0)	grids->grids[t].m_data[addr_5    * 8 + 5] = affinities[j];
		if (y != 0 && z != 0)	grids->grids[t].m_data[addr_6    * 8 + 6] = affinities[j];
		if (x * y * z != 0)		grids->grids[t].m_data[addr_7    * 8 + 7] = affinities[j];

	}
}