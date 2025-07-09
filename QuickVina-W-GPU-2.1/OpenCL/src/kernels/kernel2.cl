void get_heavy_atom_movable_coords(output_type_cl* tmp, const m_cl* m, ligand_atom_coords_cl* coords) {
	int counter = 0;
	for (int i = 0; i < m->m_num_movable_atoms; i++) {
		if (m->atoms[i].types[0] != EL_TYPE_H) {
			for (int j = 0; j < 3; j++)coords->coords[counter][j] = m->m_coords.coords[i][j];
			counter++;
		}
		else {
			//printf("\n kernel2: removed H atom coords in get_heavy_atom_movable_coords()!");
		}
	}
	//assign 0 for others
	for (int i = counter; i < MAX_NUM_OF_ATOMS; i++) {
		for (int j = 0; j < 3; j++)coords->coords[i][j] = 0;
	}
}


//Generate a random number according to step
float generate_n(__constant float* pi_map, const int step) {
	return fabs(pi_map[step]) / M_PI;
}

bool metropolis_accept(float old_f, float new_f, float temperature, float n) {
	if (new_f < old_f)return true;
	const float acceptance_probability = exp((old_f - new_f) / temperature);
	bool res = n < acceptance_probability;
	return n < acceptance_probability;
}

__kernel
void kernel2(
				const		__global		output_type_cl* ric,
						__global		m_cl* mg,
						__constant		pre_cl* pre,
						__constant		grids_cl* grids,
						__constant		random_maps* random_maps,
						__global		ligand_atom_coords_cl* coords,
						__global		output_type_cl* results,
				const	__global		mis_cl* mis,
				const				int torsion_size,
				const				int search_depth,
				const				int max_bfgs_steps,
									float center_x,
									float center_y,
									float center_z,
									float size_x,
									float size_y,
									float size_z,
						__global    	ele_cl* global_ptr,
						__global    	int* count_id,
				const 				int rilc_bfgs_enable,
				__global    	individual_container* all_lists
) {
	int gx = get_global_id(0);
	int gy = get_global_id(1);
	int gs = get_global_size(0);
	int gl = get_global_linear_id();

	float best_e = INFINITY;
	output_type_cl best_out;
	ligand_atom_coords_cl best_coords;

	vec3_cl origin = origin_init(center_x, center_y, center_z);
	vec3_cl boxsize = boxsize_init(size_x, size_y, size_z);

	for (int gll = gl;
			 gll < mis->thread;
			 gll += mis->total_wi
		)
	{
		m_cl m = *mg;

		output_type_cl tmp = ric[gll];

		change_cl g;

		output_type_cl candidate;

		// Use global memory for list
		__global individual_container* list = &all_lists[gl];
		circularvisited_init_cl(list);

		global_container g_container;
		global_init_cl(&g_container);

		for (int step = 0; step < search_depth; step++) {
			candidate = tmp;

			int map_index = (step + gll * search_depth) % MAX_NUM_OF_RANDOM_MAP;
			mutate_conf_cl(	map_index,
							&candidate,
							random_maps->int_map,
							random_maps->sphere_map,
							random_maps->pi_map,
							m.ligand.begin,
							m.ligand.end,
							m.atoms,
							&m.m_coords,
							m.ligand.rigid.origin[0],
							mis->epsilon_fl,
							mis->mutation_amplitude,
							torsion_size
			);
			if (rilc_bfgs_enable==1){
				rilc_bfgs(	&candidate,
						&g,
						&m,
						pre,
						grids,
						mis,
						torsion_size,
						max_bfgs_steps,
						true,
						&g_container,
						list,
						search_depth,
						&origin,
						&boxsize,
						global_ptr,
						count_id
				);
			}else{
				bfgs(	&candidate,
						&g,
						&m,
						pre,
						grids,
						mis,
						torsion_size,
						max_bfgs_steps,
						true,
						&g_container,
						list,
						search_depth,
						&origin,
						&boxsize,
						global_ptr,
						count_id
				);
			}
			
			float n = generate_n(random_maps->pi_map, map_index);

			if (step == 0 || metropolis_accept(tmp.e, candidate.e, 1.2, n)) {

				tmp = candidate;

				set(&tmp, &m.ligand.rigid, &m.m_coords,
					m.atoms, m.m_num_movable_atoms, mis->epsilon_fl);

				if (tmp.e < best_e) {
					if (rilc_bfgs_enable==1){
						rilc_bfgs(	&candidate,
								&g,
								&m,
								pre,
								grids,
								mis,
								torsion_size,
								max_bfgs_steps,
								false,
								&g_container,
								list,
								search_depth,
								&origin,
								&boxsize,
								global_ptr,
								count_id
						);
					}else{
						bfgs(	&candidate,
								&g,
								&m,
								pre,
								grids,
								mis,
								torsion_size,
								max_bfgs_steps,
								false,
								&g_container,
								list,
								search_depth,
								&origin,
								&boxsize,
								global_ptr,
								count_id
						);
					}

					// set
					if (tmp.e < best_e) {
						set(&tmp, &m.ligand.rigid, &m.m_coords,
							m.atoms, m.m_num_movable_atoms, mis->epsilon_fl);

						best_out = tmp;
						get_heavy_atom_movable_coords(&best_out, &m, &best_coords); // get coords
						best_e = tmp.e;
					}

				}
			}
		}

		// write the best conformation back to CPU
		results[gll] = best_out;
		coords[gll] = best_coords;
	}
}